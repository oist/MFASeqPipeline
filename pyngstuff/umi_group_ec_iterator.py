from collections import defaultdict
from pyngstuff import Umi, UmiFixed, UmiFuzzy, Strand, UmiIncrementalErrorCorrector

def umi_group_ec_iterator_return_clusters(fixedumi, cluster_iterator):
    for cluster in cluster_iterator:
        # Return unambiguously assigned pairs
        unambiguous_pairs = { Umi.combine(fixedumi, fuzzyumi): pairs for fuzzyumi, weight, pairs in cluster.unambiguous()}
        yield Umi.combine(fixedumi, cluster.principal_umi), unambiguous_pairs,
        # Return ambiguous pairs, with "None" is the clusters principal UMI
        ambiguous_pairs = { Umi.combine(fixedumi, fuzzyumi): pairs for fuzzyumi, weight, pairs in cluster.ambiguous()}
        yield None, ambiguous_pairs

def umi_group_ec_iterator(umi_group_iterator, max_weight_relative=1.0, min_weight_difference=0,
                          max_1st_mapping_distance=0, max_2nd_mapping_distance=0, max_hamming_distance=0):
    error_correctors = defaultdict(lambda: UmiIncrementalErrorCorrector(
        max_weight_relative=max_weight_relative, min_weight_difference=min_weight_difference,
        max_1st_mapping_distance=max_1st_mapping_distance, max_2nd_mapping_distance=max_2nd_mapping_distance,
        max_hamming_distance=max_hamming_distance))

    current_sequence = None
    current_umi_pos1 = None
    while True:
        # Fetch next UMI from underlying iterator, if the iterator is exhaused set next_read to None
        stop_iteration = False
        try:
            next_umi, next_pairs = next(umi_group_iterator)
        except StopIteration:
            stop_iteration = True

        # If the iteration was exhausted or the current sequence changed, flush error-correctors
        next_sequence = next_umi.seq if not stop_iteration else None
        next_umi1_pos = next_umi.fuzzy.pos1 if not stop_iteration else None
        if stop_iteration or next_sequence != current_sequence:
            # Flush error-correctors
            for fixedumi, ec in error_correctors.items():
                yield from umi_group_ec_iterator_return_clusters(fixedumi, ec.complete_clusters(flush=True))
            # Check if the iteration has ended, if not update current_reference_name
            if stop_iteration:
                current_sequence = None
                current_umi_pos1 = None
                return
            else:
                current_sequence = next_sequence
                current_umi_pos1 = -1

        # Verify that UMIs are reported in order of ascending pos1
        if next_umi1_pos < current_umi_pos1:
            raise RuntimeError("umi_group_ec_iterator expects UMIs in order of ascending pos1")
        current_umi_pos1 = next_umi1_pos

        # Add fetched UMI to the error corrector
        next_fixedumi = next_umi.fixed
        ec = error_correctors[next_fixedumi]
        weight = sum(len(pairs) for _, pairs in next_pairs.items())
        data = [pair for rawumi, pairs in next_pairs.items() for pair in pairs]
        assert len(data) == weight
        ec.add_umi(next_umi.fuzzy, weight=weight, data=data)

        # Output any complete clusters
        yield from umi_group_ec_iterator_return_clusters(next_fixedumi, ec.complete_clusters())
