from pyngstuff import postpone_read_or_pair

class umi_group_raw_iterator:
    def __init__(self, pair_iterator, umi_extractor, order_key):
        self.pair_iterator = pair_iterator
        self.umi_extractor = umi_extractor
        self.order_key = order_key
        self.umis = list()
        self.reads_or_pairs = dict()
        self.next_read_or_pair = None
        self.current_position = None

    def __iter__(self):
        return self

    def __next__(self):
        while True:
            # Fetch next read from underlying iterator, if the iterator is exhaused set next_read to None
            stop_iteration = False
            if self.next_read_or_pair is None:
                try:
                    self.next_read_or_pair = next(self.pair_iterator)
                    postpone_read_or_pair(self.next_read_or_pair)
                except StopIteration:
                    self.next_read_or_pair = None
                    stop_iteration = True

            # If the iteration was exhausted or the current position changed, return next not-yet returned group
            if not stop_iteration:
                next_pos = (self.next_read_or_pair.reference_id, self.order_key(self.next_read_or_pair))
            else:
                next_pos = None
            if stop_iteration or next_pos != self.current_position:
                if self.umis:
                    # Return next not-yet-returned group
                    umi = self.umis.pop(0)
                    return (umi, {umi: self.reads_or_pairs.pop(umi)})
                # Check if the iteration has ended, if not update current_position
                if stop_iteration:
                    self.current_position = None
                    raise StopIteration()
                else:
                    # Verify that UMIs appear in order of ascending pos1
                    if (self.current_position is not None and
                            next_pos[0] == self.current_position[0] and
                            next_pos[1] < self.current_position[1]):
                        raise RuntimeError(
                            "umi_group_raw_iterator expects the value returned by order_key to increase "
                            "monotonically within one reference sequence")
                    self.current_position = next_pos

            # Inserted next_pair into the collection of non-returned pairs and continue
            umi = self.umi_extractor(self.next_read_or_pair)
            if umi not in self.reads_or_pairs:
                self.reads_or_pairs[umi] = list()
                self.umis.append(umi)
            self.reads_or_pairs[umi].append(self.next_read_or_pair)
            self.next_read_or_pair = None
