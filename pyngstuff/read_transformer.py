import sys
import traceback
from enum import Enum
from collections import defaultdict
from recordclass import recordclass
from wrapt import ObjectProxy

class read_transformer:
    class read_output_state(Enum):
        CURRENT = 0
        POSTPONED = 1
        SKIPPED = 2
        WRITTEN = 3

    class output_buffer:
        output_slot = recordclass('slot', ['output_buffer', 'read', 'position'])
        vacant = object()
        skip = object()

        def __init__(self, output_file):
            self.output_file = output_file
            self.current_position = None
            self.postponed = list()

        def has_vacant_slots(self):
            return self.postponed == True

        def append_read(self, read):
            # Check ordering constraints
            read_pos = (read.reference_id, read.reference_start)
            self.current_position = read_pos

            # If there are postponed objects, enqueue new object at the end
            # Otherwise, output immediately.
            if self.postponed:
                self.postponed.append(read)
            else:
                self.output_file.write(read)

        def append_slot(self, reference_id, reference_start):
            # Check ordering constraints
            read_pos = (reference_id, reference_start)
            self.current_position = read_pos

            # Create placeholder object and put onto postponed list
            slot = read_transformer.output_buffer.output_slot(output_buffer=self, read=read_transformer.output_buffer.vacant,
                                                              position=(reference_id, reference_start))
            self.postponed.append(slot)
            return slot

        def fill_slot(self, slot, read):
            # Validate placeholder
            if not isinstance(slot, read_transformer.output_buffer.output_slot):
                raise TypeError('slot must be of type read_transformer.output_buffer.output_slot')
            if id(slot.output_buffer) != id(self):
                raise ValueError('slot belongs to different output_buffer')
            if id(slot.read) != id(read_transformer.output_buffer.vacant):
                raise ValueError('slot was filled in or skipped already')
            # Validate ordering constraints
            read_pos = (read.reference_id, read.reference_start)
            if slot.position != read_pos:
                raise RuntimeError("slot at position %s cannot be filled by read at position %s" %
                                   (self.format_position(slot.position), self.format_position(read_pos)))
            # Fill slot
            slot.read = read
            # Write out filled-in slots as far as possible
            self.process_postponed()

        def skip_slot(self, slot):
            # Validate placeholder
            if not isinstance(slot, read_transformer.output_buffer.output_slot):
                raise TypeError('slot must be of type read_transformer.output_buffer.output_slot')
            if id(slot.output_buffer) != id(self):
                raise ValueError('slot belongs to different output_buffer')
            if id(slot.read) != id(read_transformer.output_buffer.vacant):
                raise ValueError('slot was filled in or skipped already')
            # Set slot to special value skip
            slot.read = read_transformer.output_buffer.skip
            # Write out filled-in slots as far as possible
            self.process_postponed()

        def process_postponed(self):
            # Write postponed objects up to the first slot that is neither filled nor skipped
            while self.postponed:
                head = self.postponed[0]
                if not isinstance(head, read_transformer.output_buffer.output_slot):
                    # Not a slot, write it out
                    self.postponed.pop(0)
                    self.output_file.write(head)
                elif id(head.read) == id(read_transformer.output_buffer.vacant):
                    # A vacant slot, stop
                    break
                elif id(head.read) == id(read_transformer.output_buffer.skip):
                    # A skipped slot, just remove it and continue
                    self.postponed.pop(0)
                else:
                    # A filled slot, write the object it contains
                    self.postponed.pop(0)
                    self.output_file.write(head.read)
            # If we cleared the postponed list, we're ready to switch to a new reference sequence
            if not self.postponed:
                self.current_position = None

        def close(self):
            if self.postponed:
                raise RuntimeError('output_buffer contains unfilled slots, cannot close()')
            self.output_file.close()

        def format_position(self, position):
            reference_id, reference_start = position
            if reference_id >= 0:
                return "%s:%d" % (self.output_file.get_reference_name(reference_id), reference_start)
            else:
                assert reference_start == -1
                return "*"

    class read_proxy(ObjectProxy):
        def __init__(self, actual_read, read_transformer, debug_orphaned_slots=False):
            super(read_transformer.read_proxy, self).__init__(actual_read)
            self._self_actual_read = actual_read
            self._self_state = read_transformer.read_output_state.CURRENT
            self._self_read_transformer = read_transformer
            self._self_output_slot = None
            self._self_debug_orphaned_slots_postponed_from = None
            self._self_debug_orphaned_slots = debug_orphaned_slots

        def __del__(self):
            if self._self_debug_orphaned_slots and self._self_state == read_transformer.read_output_state.POSTPONED:
                print("WARNING: orphaned output slot at %s from" %
                      self._self_read_transformer.output_buffer.format_position(self._self_output_slot.position),
                      file=sys.stderr)
                print("".join(traceback.format_stack()[:-1]), file=sys.stderr)
                print("slot originally created from", file=sys.stderr)
                print("".join(self._self_debug_orphaned_slots_postponed_from), file=sys.stderr)

        @property
        def output_state(self):
            return self._self_state

        def output_skip_read(self, category=None):
            if self._self_state == read_transformer.read_output_state.CURRENT:
                # No slot allocated yet, nothing to do
                self._self_state = read_transformer.read_output_state.SKIPPED
                self._self_read_transformer.stats_categories[category] += 1
            elif self._self_state == read_transformer.read_output_state.POSTPONED:
                if self._self_read_transformer.output_buffer is not None:
                    # Remove already allocated slot from the output stream
                    self._self_read_transformer.output_buffer.skip_slot(self._self_output_slot)
                    self._self_output_slot = None
                self._self_state = read_transformer.read_output_state.SKIPPED
                self._self_read_transformer.stats_categories[category] += 1
            elif self._self_state == read_transformer.read_output_state.SKIPPED:
                # Already skipped, nothing to do
                pass
            elif self._self_state == read_transformer.read_output_state.WRITTEN:
                # Cant go back to skipped
                raise RuntimeError("cannot writen an already skipped read")

        def output_postpone_read(self):
            if self._self_state == read_transformer.read_output_state.CURRENT:
                if self._self_read_transformer.output_buffer is not None:
                    # Insert a slot into the output stream that can be filled later
                    if self._self_debug_orphaned_slots:
                        # Rememver were the read was postponed, this is later reported should the slot be orphaned
                        self._self_debug_orphaned_slots_postponed_from = traceback.format_stack()[:-2]
                    self._self_output_slot = self._self_read_transformer.output_buffer.append_slot(
                        self._self_actual_read.reference_id,
                        self._self_actual_read.reference_start)
                self._self_state = read_transformer.read_output_state.POSTPONED
            elif self._self_state == read_transformer.read_output_state.POSTPONED:
                # Already postponed, do nothing
                pass
            elif self._self_state == read_transformer.read_output_state.SKIPPED:
                # Cant go back to postponed
                raise RuntimeError("cannot postpone an already skipped read")
            elif self._self_state == read_transformer.read_output_state.WRITTEN:
                # Cant go back to postponed
                raise RuntimeError("cannot postpone an already written read")

        def output_write_read(self, category="output"):
            if self._self_state == read_transformer.read_output_state.CURRENT:
                # Read is still the last thing read from the input, just write to output
                self._self_read_transformer.output_buffer.append_read(self.__self_actual_read)
                self._self_state = read_transformer.read_output_state.WRITTEN
                self._self_read_transformer.stats_category[category] += 1
            elif self._self_state == read_transformer.read_output_state.POSTPONED:
                # Replace previously inserted output slot with the read to be written
                if self._self_read_transformer.output_buffer is not None:
                    self._self_read_transformer.output_buffer.fill_slot(self._self_output_slot, self._self_actual_read)
                    self._self_output_slot = None
                self._self_state = read_transformer.read_output_state.WRITTEN
                self._self_read_transformer.stats_categories[category] += 1
            elif self._self_state == read_transformer.read_output_state.SKIPPED:
                # Cant go back to postponed
                raise RuntimeError("cannot write an already skipped read")
            elif self._self_state == read_transformer.read_output_state.WRITTEN:
                # Cant go back to postponed
                raise RuntimeError("cannot write an already written read")

    def __init__(self, input_iterator, output_file, usable_categories={"output"}, output_categories={"output"},
                 debug_orphaned_slots=False):
        self.input_iterator = input_iterator
        self.output_buffer = read_transformer.output_buffer(output_file) if output_file else None
        self.current_sequence = None
        self.current_read = None
        # Statistics
        self.stats_read = 0
        self.stats_categories = defaultdict(lambda: 0)
        self.usable_categories = frozenset(usable_categories)
        self.output_categories = frozenset(output_categories)
        # Debugging
        self.debug_orphaned_slots = debug_orphaned_slots

    def __iter__(self):
        return self

    def __next__(self):
        # If the read wasn't already skipped, postponed or written, skip now
        if self.current_read is not None and self.current_read.output_state == read_transformer.read_output_state.CURRENT:
            self.current_read.output_skip_read()
            self.current_read = None

        # Fetch read, wrap inside a read_proxy
        try:
            self.current_read = read_transformer.read_proxy(next(self.input_iterator), self,
                                                            debug_orphaned_slots=self.debug_orphaned_slots)
            self.stats_read += 1
        except StopIteration:
            self.current_read = None
            raise
        return self.current_read

    @property
    def stats_total(self):
        return sum(n for c, n in self.stats_categories.items())

    @property
    def stats_output(self):
        return sum(n for c, n in self.stats_categories.items() if c in self.output_categories)

    @property
    def stats_usable(self):
        return sum(n for c, n in self.stats_categories.items() if c in self.usable_categories)

    @property
    def stats_ignored(self):
        return sum(n for c, n in self.stats_categories.items() if c not in self.usable_categories)

    def stats_category(self, category):
        return self.stats_categories[category]

    def print_statistics(self, file=sys.stdout):
        read, output, usable, ignored = self.stats_read, self.stats_output, self.stats_usable, self.stats_ignored
        print("*** Read statistics")
        print("Input contained %d reads, of these" % read, file=file)
        def print_line(k, c, n):
            print("%s%d (%.0f%%) were %s" % (" " * k, n, 100.0 * n / read if read > 0 else float('nan'), c), file = file)
        def print_lines(k, generator):
            for c, n in generator:
                if n > 0:
                    print_line(k, c, n)

        if output > 0 and not self.output_categories.issubset(self.usable_categories):
            # Output categories are not a subset of usable categories, report separately
            print_line(2, "output", output)
            print_lines(2, ((c, n) for c, n in self.stats_categories.items() if c in self.output_categories and c != "output"))
        if usable > 0:
            # Nest output categories inside usable
            print_line(2, "usable", usable)
            if output > 0 and self.output_categories.issubset(self.usable_categories):
                print_line(4, "output", output)
                print_lines(6, ((c, n) for c, n in self.stats_categories.items() if c in self.output_categories and c != "output"))
            print_lines(4, ((c, n) for c, n in self.stats_categories.items() \
                            if c in self.usable_categories and c not in self.output_categories))
        if ignored > 0:
            print_line(2, "ignored", ignored)
            print_lines(4, ((c, n) for c, n in self.stats_categories.items() if c not in self.usable_categories))
        file.flush()

    def close(self, statistics_file=sys.stdout):
        self.input_iterator = None
        if self.output_buffer is not None:
            self.output_buffer.close()
        if statistics_file is not None:
            self.print_statistics(file=statistics_file)

def postpone_read(read):
    if isinstance(read, read_transformer.read_proxy):
        read.output_postpone_read()

def skip_read(read, category=None):
    if isinstance(read, read_transformer.read_proxy):
        read.output_skip_read(category)

def write_read(read, category="output"):
   if isinstance(read, read_transformer.read_proxy):
       read.output_write_read(category)
