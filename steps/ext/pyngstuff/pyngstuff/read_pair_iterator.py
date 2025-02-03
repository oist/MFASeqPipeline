import sys
from collections import deque
from collections import defaultdict
from pyngstuff import postpone_read, skip_read, write_read

class read_pair_iterator:
    """
    Iterate over the read in a BAM file as mate pairs

    Reads are *roughly* reported in the order they appear in in the input file.
    Non-overlapping pairs are output ordered by ascending mapping coordinates.
    For overlapping pairs no particular order is guaranteed. To have overlapping
    pairs output in a deterministic order, use `read_pair_ordered_iterator`.

    Mate pairs where only one mate is mapped are reported as singleton reads if
    skip_unpaired is set to false, or skipped if set to true.

    Mate pairs where both mates are mapped, but in a discordant fasion (i.e. not
    flagged as proper in the BAM file, meaning their distance was too large, they
    ended up on different chromosomes, or had the wrong orientations) are skipped.

    For each pair an instance of `read_pair` is returned, which allows accessing
    the constituent reads using various access like `read1`, `read2`, `left`, `right`,
    and similar. See `read_pair` for details.
    """

    # Queue entries
    class read_pair:
        def __init__(self, read_left, read_right):
            self.left = read_left
            self.right = read_right

        @property
        def reference_name(self):
            return self.left.reference_name

        @property
        def reference_id(self):
            return self.left.reference_id

        @property
        def read1(self):
            if self.left is not None and self.left.is_read1:
                return self.left
            elif self.right is not None and self.right.is_read1:
                return self.right
            else:
                return None

        @read1.setter
        def read1(self, read):
            if not read.is_read1:
                raise RuntimeError("replacement read for read1 must also be a read1")
            if self.left is not None and self.left.is_read1:
                self.left = read
            elif self.right is not None and self.right.is_read1:
                self.right = read
            else:
                raise RuntimeError("pair contains no read1 to replace")

        @property
        def read2(self):
            if self.left is not None and not self.left.is_read1:
                return self.left
            elif self.right is not None and not self.right.is_read1:
                return self.right
            else:
                return None

        @read2.setter
        def read2(self, read):
            if read.is_read1:
                raise RuntimeError("replacement read for read2 must also be a read2")
            if self.left is not None and not self.left.is_read1:
                self.left = read
            elif self.right is not None and not self.right.is_read1:
                self.right = read
            else:
                raise RuntimeError("pair contains no read2 to replace")

        @property
        def forward(self):
            if self.left is not None and not self.left.is_reverse:
                return self.left
            elif self.right is not None and not self.right.is_reverse:
                return self.right
            else:
                return None

        @forward.setter
        def forward(self, read):
            if not read.is_read1:
                raise RuntimeError("replacement read for forward read also be a forward read")
            if self.left is not None and not self.left.is_reverse:
                self.left = read
            elif self.right is not None and not self.right.is_reverse:
                self.right = read
            else:
                raise RuntimeError("pair contains no forward read to replace")

        @property
        def reverse(self):
            if self.left is not None and self.left.is_reverse:
                return self.left
            elif self.right is not None and self.right.is_reverse:
                return self.right
            else:
                return None

        @reverse.setter
        def reverse(self, read):
            if not read.is_read1:
                raise RuntimeError("replacement read for reverse read also be a reverse read")
            if self.left is not None and self.left.is_reverse:
                self.left = read
            elif self.right is not None and self.right.is_reverse:
                self.right = read
            else:
                raise RuntimeError("pair contains no forward read to replace")

    def __init__(self, read_iterator, skip_unpaired=True, skip_secondary=True, skip_supplementary=True):
        # Settings
        self.skip_unpaired = skip_unpaired
        self.skip_secondary = skip_secondary
        self.skip_supplementary = skip_supplementary
        # Input iterator
        self.iterator = read_iterator
        # Incomplete pairs
        self.current_reference_id = None
        self.current_reference_name = None
        self.incomplete_pairs = dict()

    def __iter__(self):
        return self

    def __next__(self):
        # Repeat until a valid pair is found
        pair = None
        while True:

            # 1. Fetch read
            #
            try:
                r = next(self.iterator)
            except StopIteration:
                # Reference id changes to None
                self.on_sequence_will_change()
                self.current_reference_id = None
                self.current_reference_name = None
                raise
            # Check if sequence name changes, run checks if it does
            if self.current_reference_id is not None and self.current_reference_id != r.reference_id:
                self.on_sequence_will_change()
            self.current_reference_id = r.reference_id
            self.current_reference_name = r.reference_name

            # 2. Attempt to form a pair
            #
            # (A) If both mates are unmapped, the pair is ignored.
            # (B) If exactly one mate is mapped, the pair is accepted, and the other mate
            # set to "None". If skip_unpaired was set, the pair is dropped instead.
            # (C) Supplementary alignments (i.e., if a single alignment is split over
            # multiple alignments) are skipped unless skip_supplementary is False. Its is
            # unclear how pairs are identified in the presence of supplementary alignments
            # (D) Secondary alignments are skipped unless skip_secondary is set to False.
            # (E) If both mates are mapped but discordantely, the pair is ignored.
            # (F) If both are mapped concordantly,
            # the pair is (naturally) accepted.
            if r.is_supplementary and self.skip_supplementary:
                # C in the skip_supplementary case. Ignore.
                skip_read(r, "supplementary alignment")
                continue
            elif r.is_unmapped:
                # A, B for unmapped mate. Ignore.
                skip_read(r, "read unmapped")
                continue
            elif r.is_secondary and self.skip_secondary:
                # D in the skip_secondary case. Ignore.
                skip_read(r, "secondary alignment")
                continue
            elif r.mate_is_unmapped:
                if not self.skip_unpaired:
                    # B for mapped mate. Create pair entry with missing mate, enqueue.
                    # Since we dont return the pair right away, postpone it
                    postpone_read(r)
                    pair = read_pair_iterator.read_pair(r, None)
                    continue
                else:
                    # B for mapped mate. Ignore because skip_unpaired is set
                    skip_read(r, "mate unmapped")
                    continue
            elif not r.is_proper_pair or r.reference_name != r.next_reference_name:
                # C. Ignore.
                skip_read(r, "discordant pair")
                continue
            else:
                # D. Form unique key that we use to identify read pairs
                pair_key = (r.query_name, r.reference_id, r.reference_start)
                if pair_key not in self.incomplete_pairs:
                    # Mate was not yet seen, create incomplete_pairs entry
                    postpone_read(r)
                    mate_key = (r.query_name, r.next_reference_id, r.next_reference_start)
                    self.incomplete_pairs[mate_key] = read_pair_iterator.read_pair(r, None)
                    continue
                else:
                    # Mate was already seen, found a pair
                    pair = self.incomplete_pairs.pop(pair_key)
                    pair.right = r
                    break

            # 3. Validate pair and return
            #
            # Unless both reads map to the same sequence and one of the reads maps to the forward strand, and the other
            # to the reverse strand, the pair is skipped. Similarly, if the forward read lies fully to the right of the
            # reverse read (i.e., if they look like this: <--| |-->), the pair is skipped.
            # Note that since we clear the read queue when the sequence name changes, the two reads must necessarily
            # map to the same sequence for us to get here
            assert pair.left.reference_id == pair.right.reference_id
            if ((pair.left is not None) and (pair.right is not None) and
                    (pair.left.is_reverse == pair.right.is_reverse)):
                # Concordant pairs expected to have orientation FR or RF
                skip_read_or_pair(pair, 'pair orientation')
                continue
            if ((pair.left is not None) and (pair.right is not None) and
                    (pair.left.reference_end <= pair.right.reference_start) and pair.left.is_reverse):
                # Concordant pairs expected to point towards one another
                skip_read_or_pair(pair, 'pair direction')
                continue
            # Found a valid pair
            break
        return pair

    # Iterator reached the next sequence.
    # No proper pairs should extend over multiple sequences, so warn if that is the case
    def on_sequence_will_change(self):
        if self.incomplete_pairs:
            print("WARNING: %d incomplete concordant mate pairs left af the end of %s, skipping" %
                  (len(self.incomplete_pairs), self.current_reference_name),
                  file=sys.stderr)
            for pair in self.incomplete_pairs.values():
                skip_read_or_pair(pair, "incomplete pair")
        self.incomplete_pairs.clear()

def postpone_read_or_pair(read):
    if isinstance(read, read_pair_iterator.read_pair):
        if read.left is not None:
            postpone_read(read.left)
        if read.right is not None:
            postpone_read(read.right)
    else:
        postpone_read(read)

def skip_read_or_pair(read, reason=None):
    if isinstance(read, read_pair_iterator.read_pair):
        if read.left is not None:
            skip_read(read.left, reason)
        if read.right is not None:
            skip_read(read.right, reason)
    else:
        skip_read(read, reason)

def write_read_or_pair(read):
    if isinstance(read, read_pair_iterator.read_pair):
        if read.left is not None:
            write_read(read.left)
        if read.right is not None:
            write_read(read.right)
    else:
        write_read(read)
