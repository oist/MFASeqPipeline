from sortedcontainers import SortedDict
from pyngstuff import read_pair_iterator, postpone_read_or_pair

class read_pair_ordered_iterator:
    def __init__(self, read_iterator, skip_unpaired=True):
        self.pair_iterator = read_pair_iterator(read_iterator, skip_unpaired=skip_unpaired)
        self.ordered_pairs = SortedDict()
        self.next_pair = None
        self.current_reference_id = None

    def __iter__(self):
        return self

    def __next__(self):
        while True:
            # Fetch next read from underlying iterator, if the iterator is exhaused set next_pair to None
            stop_iteration = False
            if self.next_pair is None:
                try:
                    self.next_pair = next(self.pair_iterator)
                    # Since we dont return the pair right away, postpone it
                    postpone_read_or_pair(self.next_pair)
                except StopIteration:
                    self.next_pair = None
                    stop_iteration = True

            # If the iteration was exhausted or the reference_id changed, return all not-yet returned pairs in order
            if stop_iteration or (self.next_pair.reference_id != self.current_reference_id):
                if self.ordered_pairs:
                    # Return the first pair of the pairs with the lower order key, remove returned pair, and remove
                    # empty ordered_pairs entries.
                    order_key, pairs = self.ordered_pairs.peekitem(0)
                    pair = pairs.pop(0)
                    if not pairs:
                        del self.ordered_pairs[order_key]
                    return pair
                # Check if the iteration has ended, if not update current_reference_name
                if stop_iteration:
                    self.current_reference_id = None
                    raise StopIteration()
                else:
                    self.current_reference_id = self.next_pair.reference_id

            # All pairs that lie fully to the left of next_pair can be returned now, return them in order
            if self.ordered_pairs:
                order_key, pairs = self.ordered_pairs.peekitem(0)
                if order_key <= self.next_pair.forward.reference_start:
                    # Return the first pair of the pairs with the lower order key, remove returned pair, and remove
                    # empty ordered_pairs entries.
                    pair = pairs.pop(0)
                    if not pairs:
                        del self.ordered_pairs[order_key]
                    return pair

            # Inserted next_pair into the collection of non-returned pairs and continue
            order_key = self.next_pair.reverse.reference_end
            self.ordered_pairs.setdefault(order_key, list())
            self.ordered_pairs[order_key].append(self.next_pair)
            self.next_pair = None
