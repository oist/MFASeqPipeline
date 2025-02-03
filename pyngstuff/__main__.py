import argparse
import copy
import pysam
import sys
from collections import defaultdict
from pyngstuff import Umi, Strand, UmiTagSource, UmiTagSourceAutodetect, \
    skip_read, write_read, read_transformer, read_pair_ordered_iterator, \
    umi_group_raw_iterator, umi_group_ec_iterator, progress_reporter

# Parse arguments
ap = argparse.ArgumentParser(description="remove PCR duplicates stemming from the same original molecule")
ap.add_argument('input_bam')
ap.add_argument('--output-bam', metavar="FILE", help="file to write flagged or filtered reads to")
ap.add_argument('--output-readdist', metavar="FILE", help="fiile to write per-UMI read count distribution to")
#ap.add_argument('--output-umis', metavar="FILE", help="file to a table of UMIs and read counts to")
ap.add_argument('--include-duplicates', help="include duplicats (in the output (flagged as duplicates)", action="store_true")
ap.add_argument('--include-ambiguous', help="include reads without unambiguous original molecule (i.e. UMI)  ", action="store_true")
ap.add_argument('--paired', help="input file containes paired-end data", action="store_true")
ap.add_argument('--umitag-source', help="location of each read's UMI tag", default="auto", type=UmiTagSource.parse)
ap.add_argument('--ignore-cellbc', help="allow duplicates to stem from different cells", action="store_true")
ap.add_argument('--ignore-strand', help="allow duplicates to map to different strands", action="store_true")
ap.add_argument('--ignore-position', help="allow duplicates to map to different positions (set max-map-dist to infinity)", action="store_true")
ap.add_argument('--ignore-umitag', help="allow duplicates to have distinct UMI tags (i.e., ignore UMI tags)", action="store_true")
ap.add_argument('--max-map-dist', metavar="DIST", help="maximum difference between mapping distances of duplicates",
                type=int, default=1)
ap.add_argument('--max-umitag-dist', metavar="DIST", help="maximum number of differences between UMI tags of duplicates",
                type=int, default=1)
#ap.add_argument('--min-overlap', metavar="DIST", help="minimal overlap between duplicates", default=None)
ap.add_argument('--progress', help="show progress", action="store_true")
ap.add_argument('--progress-bar', help="show progress bar", action="store_true")
args = ap.parse_args()

# Handle argument dependencies
if args.ignore_umitag and args.ignore_cellbc:
    # Use a dummy source that always returns None if we're ignoring the tag anyway
    args.umitag_source = lambda read: (None, None)
if args.ignore_cellbc and isinstance(args.umitag_source, UmiTagSourceAutodetect):
    # If we're auto-detecting the UMI tag source, let the auto-detector know if we're ignoring cell barcodes anyway
    args.umitag_source.cellbarcode = False

# Auto-detect settings
input_file = pysam.AlignmentFile(args.input_bam, mode="rb")
if isinstance(args.umitag_source, UmiTagSourceAutodetect):
    read = next(input_file, None)
    if read is not None:
        args.umitag_source(read)
    else:
        raise RuntimeError("input file is empty, cannot auto-detect UMI tag source")
input_file.close()

# Re-open input file, enable progress output if desired
input_file = pysam.AlignmentFile(args.input_bam, mode="rb")
if args.progress or args.progress_bar:
    def statusline(_, total):
        processed, usable = read_transformer.stats_total, read_transformer.stats_usable
        output, duplicate = read_transformer.stats_category("unique"), read_transformer.stats_category("duplicate")
        reads_per_umi = float(output + duplicate) / float(output) if output > 0 else float('nan')
        return "%s%d%s reads, %d (%.1f%%) useable, %.2f reads/UMI" % (
            "[%.1f%%] " % (100.0*processed/total) if total is not None and total > 0 else "",
            processed, "/%d" % total if total is not None else "",
            usable, 100.0*usable/processed if processed > 0 else float('nan'), reads_per_umi
        )

    total_reads = None
    try:
        total_reads = input_file.mapped + input_file.unmapped
    except:
        print("*** BAM file not indexed, will not show progress bar", file=sys.stdout, flush=True)

    read_iterator = progress_reporter(input_file, total=total_reads, every=10000, statusline=statusline,
                                      target=sys.stdout, show_progressbar=args.progress_bar)
else:
    read_iterator = input_file

# Open output file
if args.output_bam is not None:
    output_file = pysam.AlignmentFile(args.output_bam, mode="wb", header=copy.deepcopy(dict(input_file.header)))
else:
    output_file = None

# First step in the read processing pipeline is the read processor which allows
# (modified) reads to be written to the output file while keeping the file ordered,
# and tracks statistics about how many reads were encountered in different categories
usable_categories = {"unique", "duplicate"}
output_categories = set({})
if args.output_bam:
    output_categories |= {"unique"}
    if args.include_duplicates:
        output_categories |= {"duplicate"}
    if args.include_ambiguous:
        output_categories |= {"ambiguous umi"}
if args.include_ambiguous:
    usable_categories |= {"ambiguous umi"}
read_transformer = read_transformer(read_iterator, output_file, usable_categories=usable_categories, output_categories=output_categories)

# Setup read processing pipeline that returns groups of identical UMIs extracted from either paired or single-end reads
if args.paired:
    def order_key(pair):
        return pair.reverse.reference_end

    def umi_extractor(pair):
        assert pair.forward.query_name == pair.reverse.query_name
        (umitag, cellbc) = args.umitag_source(pair.left)
        strand = Strand.PLUS if not pair.read1.is_reverse else Strand.MINUS
        return Umi(sample=cellbc if not args.ignore_cellbc else None,
                   seq=pair.reference_name,
                   strand=strand if not args.ignore_strand else None,
                   pos1=pair.reverse.reference_end if not args.ignore_position else None,
                   pos2=pair.forward.reference_start if not args.ignore_position else None,
                   tag=umitag if not args.ignore_umitag else None)

    # Iterate over read pairs in the order defined by order_key (increasing pair.reverse.reference_end),
    # and group by UMI as defined by the UMI extractor (reference_name, strand, forward.reference_start, reverse.reference.end, umitag)
    ordered_pair_i = read_pair_ordered_iterator(read_transformer, skip_unpaired=True)
    raw_umigroup_i = umi_group_raw_iterator(ordered_pair_i, umi_extractor=umi_extractor, order_key=order_key)

else:
    def order_key(read):
        return read.reference_start

    def umi_extractor(read):
        (umitag, cellbc) = args.umitag_source(read)
        strand = Strand.PLUS if not read.is_reverse else Strand.MINUS
        return Umi(sample=cellbc if not args.ignore_cellbc else None,
                   seq=read.reference_name,
                   strand=strand if not args.ignore_strand else None,
                   pos1=read.reference_start if not args.ignore_position else None,
                   tag=umitag if not args.ignore_umitag else None)

    # Iterate over reads in the order defined by order_key (read.reference_start),
    # and group by UMI as defined by the UMI extractor (reference_name, strand, reference_start, umitag)
    raw_umigroup_i = umi_group_raw_iterator(read_transformer, umi_extractor=umi_extractor, order_key=order_key)

# Add error-correction step to the processing pipeline
if (not args.ignore_umitag and args.max_umitag_dist > 0) or (not args.ignore_position and args.max_map_dist > 0):
    ec_umigroup_i = umi_group_ec_iterator(raw_umigroup_i,
                                          max_hamming_distance=args.max_umitag_dist if not args.ignore_umitag else 0,
                                          max_1st_mapping_distance=args.max_map_dist if not args.ignore_position else 0,
                                          max_2nd_mapping_distance=args.max_map_dist if not args.ignore_position else 0)
else:
    ec_umigroup_i = raw_umigroup_i

# Collect per-UMI read-count distribution
if args.output_readdist is not None:
    readdist = defaultdict(lambda: 0)

# Iterate through reads grouped by UMIs, flags all except one pair per read as PCR duplicate
with read_iterator as _:
    for umi, groups in ec_umigroup_i:
        # Iterate over the individual reads with the same (error-corrected) UMI, add tags and output
        # We do this even if we don't output a BAM file to ensure accurate read statistics
        for rawumi, reads_or_pairs in groups.items():
            for i, read_or_pair in enumerate(reads_or_pairs):
                is_principal_pair = (i == 0 and rawumi == umi)
                reads = (read_or_pair.left, read_or_pair.right) if args.paired else (read_or_pair,)
                for read in reads:
                    if umi is not None:
                        # Mark a single read (pair) as not being a duplicate,
                        # tag all reads with the error-corrected full UMI (tag MI),
                        # and with the number of reads that observe this UMI (tag XD).
                        # We also set the CB tag to the cell barcode if there is one.
                        if args.output_bam:
                            read.is_duplicate = not is_principal_pair
                            if umi.sample is not None:
                                read.set_tag(tag="CB", value_type="Z", value=umi.sample)
                            if rawumi.tag is not None:
                                read.set_tag(tag="OX", value_type="Z", value=rawumi.tag)
                            if umi.tag is not None:
                                read.set_tag(tag="RX", value_type="Z", value=umi.tag)
                            read.set_tag(tag="MI", value_type="Z", value=str(umi))
                            read.set_tag(tag="XD", value_type="i", value=len(reads_or_pairs))
                        if is_principal_pair:
                            write_read(read, "unique")
                        elif args.include_duplicates:
                            write_read(read, "duplicate")
                        else:
                            skip_read(read, "duplicate")
                    else:
                        # Reads that couldn't unambiguously be error-corrected
                        if args.output_bam:
                            read.is_duplicate = True
                            if umi.sample is not None:
                                read.set_tag(tag="CB", value_type="Z", value=umi.sample)
                            if rawumi.tag is not None:
                                read.set_tag(tag="OX", value_type="Z", value=rawumi.tag)
                        if args.include_ambiguous:
                            write_read(read, "ambiguous umi")
                        else:
                            skip_read(read, "ambiguous umi")

        # Collect distribution of reads per UMI
        if umi is not None and args.output_readdist:
            readcount = sum(len(reads_or_pairs) for reads_or_pairs in groups.values())
            readdist[readcount] += 1

# We're done, close self.output_file
read_transformer.close()
input_file.close()

# Output  per-UMI read-count distribution
if args.output_readdist is not None:
    with open(args.output_readdist, 'w') as file:
        file.write("reads\tnumis\n")
        for rc in sorted(readdist.keys()):
            file.write("%d\t%d\n" % (rc, readdist[rc]))
