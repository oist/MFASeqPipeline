from collections import namedtuple
from sortedcontainers import SortedDict
from math import floor, inf
from pyngstuff import hamming_distance


class UmiCluster:
    """ A cluster of similar UMIs, some of which may be ambiguous

    A cluster consists of a set of UMIs (plus weights) which are all
    considered to be errorneous versions of the cluster's ``principal_umi``.
    Some of the cluster members may be ambiguous, i.e. may also have
    successor node in other clusters.

    The member ``total_weight`` contains the sum of weights of all
    unambiguous cluster members.

    The member function ``unambiguous()`` returns an iterator over
    all pairs ``(umi, weight)`` for the unambiguous UMIs in the cluster

    The member function ``ambiguous()`` returns an iterator over
    all pairs ``(umi, weight)`` for the ambiguous UMIs in the cluster

    The principal UMI's node in the ``UmiDag``instace before the cluster was
    extracted from the DAG can be access via the ``terminal`` member.
    Note that edges towards nodes not in the cluster are removed when
    a cluster is extracted from an ``UmiDag`` instance!
    """

    def __init__(self, terminal):
        self.principal_umi = terminal.umi
        self.principal_weight = terminal.weight
        self.total_weight = None
        self.terminal = terminal
        # Filled in by ``extract_from_dag()``
        self.unambiguous_uminodes = None
        self.ambiguous_uminodes = None
        # Temporaries used by ``extract_from_dag()``
        self.dag = None

    def unambiguous(self):
        """ Iterates of ``(UmiFuzzy, weight, data)`` pairs for all unambiguous UMIs in the cluster """
        for uminode in self.unambiguous_uminodes:
            yield (uminode.umi, uminode.weight, uminode.data)

    def ambiguous(self):
        """ Iterates of ``(UmiFuzzy, weight, data)`` pairs for all ambiguous UMIs in the cluster,
            i.e. UMIs that were similar to UMIs within as well as outside of the cluster """
        for uminode in self.ambiguous_uminodes:
            yield (uminode.umi, uminode.weight, uminode.data)

    def extract_from_dag(self, dag):
        """ Extracts the cluster from its ``UmiDag``. Normally by ``UmiDag`` after construction """
        try:
            # Visitors need to access the DAG
            self.dag = dag
            # Scan predecessors of terminal uminode, determine unambiguous and ambiguous
            # cluster members (the latter have successors without as well as outside of
            # the cluster), and remove all uminodes we encounter from the DAG
            self.unambiguous_uminodes = set()
            self.ambiguous_uminodes = set()
            UmiDag.visit_predecessors(self.terminal, self.extraction_visitor)
            # The first scan finds ambiguous nodes, but does not mark their predecessors
            # as ambiguous as well. Do that now, and also remove all edges that leave the
            # current cluster
            for uminode in frozenset(self.ambiguous_uminodes):
                # Remove the outgoing (i.e. successor) edges that are not
                # part of the cluster, so that the cluster is self-contained
                assert len(uminode.successors) > 1
                for s in frozenset(uminode.successors):
                    if s not in self.ambiguous_uminodes and s not in self.unambiguous_uminodes:
                        uminode.successors.remove(s)
                        s.predecessors.remove(uminode)
                # Add predecessors to the set of ambiguous uminodes
                UmiDag.visit_predecessors(uminode, self.indirect_ambiguity_visitor, include_self=False)
            # Compute total weight
            self.total_weight = sum(w for _, w, _ in self.unambiguous())
        finally:
            self.dag = None

    def extraction_visitor(self, uminode):
        # Once a node is marked unambiguous it remains unambiguous. Since we descend
        # immediately after marking a node as unambiguous, if we re-enter the node
        # we don't have to descend again, so we may exit early.
        if uminode in self.unambiguous_uminodes:
            return False
        # Only if *all* the uminode's successors are unambiguous members of the cluster,
        # the uminode itself is too (in other words, no successor is member of another
        # cluster or itself ambiguous)
        if len(uminode.successors) > 1:
            if all(s in self.unambiguous_uminodes for s in uminode.successors):
                # Actually unambiguous, remove previously set ambiguity mark
                self.unambiguous_uminodes.add(uminode)
                if uminode in self.ambiguous_uminodes:
                    # We reached the node before and tentatively marked in ambiguous, undo that
                    self.ambiguous_uminodes.remove(uminode)
                else:
                    # First the we reach the node, remove from DAG
                    self.dag.remove_umi(uminode)
                # Our predecessors may switch from ambiguous to unambiguous, so descend in any case
                return True
            elif uminode in self.ambiguous_uminodes:
                # Previously visited this node, so don't descend to predecessors again
                return False
            else:
                # Tentatively ambiguous, removed again later if the other successors
                # are found to lie within the cluster as well (see case above)
                self.ambiguous_uminodes.add(uminode)
                # First time we reached this node, so remove from DAG and descend to predecessors
                self.dag.remove_umi(uminode)
                return True
        else:
            # Single-successor uminodes are always unambiguous
            self.unambiguous_uminodes.add(uminode)
            # Single-succesoor uminodes are visited only once, since we take care here
            # to not descend twice into the same node. So this must be the first visit,
            # otherwise we'd have bailed out at the top already. So descend!
            self.dag.remove_umi(uminode)
            return True

    def indirect_ambiguity_visitor(self, uminode):
        if uminode in self.ambiguous_uminodes:
            # Node was visited already, don't descend into predecessors
            return False
        assert uminode in self.unambiguous_uminodes
        self.unambiguous_uminodes.remove(uminode)
        self.ambiguous_uminodes.add(uminode)


class UmiDag:
    """A *directed acyclic graph* (DAG) to represent UMI similarity
`
    The ``UmiDag`` class represents a mutable graph whose nodes represent UMIs,
    and edges connect UMIs that are similar enough to be probably stem from the
    same original molecule. The class is designed to support *incremental* merging
    of similar UMIs, i.e. UMIs are added to the graph in order of ascending genomic
    coordinates, and once all UMIs in a cluster of similar UMIs are far enough
    "behind" the current position , it is output and its nodes removed from the DAG.

    **Edge Direction**

    Edges are oriented, and point from a predecessor node with lower weight towards
    a successor node with larger weight. Conceptually, edges thus point from an
    errorneous version of a UMI towards the error free (or a less errorneos) version
    into which it is supposed to be merged. The requirement that the predecessor's
    weight must not exceed the successor's weight ensures that the graph is acyclic.
    To allow edges between UMIs with equal weights while still ensuring acyclicity,
    we allow for such pairs of UMIs edges only from the smaller to larger UMIs in the
    lexicographic order of UMIs as triples (pos1, pos2, tag).

    **Terminal Nodes**

    Nodes without successors represent (putative) error-free UMIs, and are the target
    into which their surrounding cluster of similar UMIs (i.e. direct or indirect
    predecessors) are merged. To facilitate fast detection of clusters of similar UMIs
    that can no longer grow and may thus be output, ``UmiDag`` tracks the current set
    of terminal nodes, ordered by the rightmost position of any of their predecessors
    (``pos1max``).

    **Similarity Index***

    In addition to the ``nodes`` dictionary that maps UMIs to their corresponding
    ``UmiNode`` instance, ``UmiDag`` also maintains an index that allows fast enumration
    of all UMIs that are similar to a given UMI, meaning that the differences between
    their positions and tags are do not exceed a set of given thresholds (distances
    between tags are number of differing characters, i.e. hamming distance)

    **Usage**::

        dag = UmiDag(...)
        for umi, weight in umis_sorted_by_pos1:
            node = dag.add_umi(umi, weight, additional_data)
            for other_node in dag.find_potentially_similar_umis(...):
                if umis_are_really_similar:
                    dag.add_edge(node, other_node)
            for cluster in dag.find_and_remove_clusters_before(umi.pos1 - slack):
                output(cluster)
    """

    class UmiNode(namedtuple('UmiNode', ['umi', 'weight', 'pos1max', 'predecessors', 'successors', 'data'])):
        """ A node in the DAG """

        def __new__(cls, umi, weight, data=None):
            return super().__new__(cls, umi, weight, UmiDag.MutableValue(umi.pos1), set(), set(), data)

        def __hash__(self):
            return id(self)

        def __eq__(self, other):
            return id(self) == id(other)

        def __str__(self):
            return "UmiNode(%s, w=%s) [pos1max=%d, predecessors={%s}, successors={%s}]" % (
                self.umi, self.weight, self.pos1max.value,
                ",".join(str(p.umi) for p in self.predecessors),
                ",".join(str(p.umi) for p in self.successors),
            )

        def __repr__(self):
            return str(self)

    class MutableValue:
        """ Wrapper around a single value, used for single mutable values in ``UmiNode`` """

        def __init__(self, value):
            self.value = value

        def __str__(self):
            return str(self.value)

    @staticmethod
    def visit_successors(uminode, visitor, include_self=True):
        """ Run the ``visitor`` function for all direct or indirect successors of ``uminode``

        If ``include_self`` is ``True``, the visitor is first called for ``uminode`` itself
        """
        c = visitor(uminode) if include_self else None
        if type(c) != bool or c:
            # Ascend if visitor returned anything other than ``False``.
            for s in uminode.successors:
                UmiDag.visit_successors(s, visitor, include_self=True)

    @staticmethod
    def visit_predecessors(uminode, visitor, include_self=True):
        """ Run the ``visitor`` function for all direct or indirect predecessors of ``uminode``

        If ``include_self`` is ``True``, the visitor is first called for ``uminode`` itself
        """
        c = visitor(uminode) if include_self else None
        if type(c) != bool or c:
            # Descend if visitor returned anything other than ``False``.
            for p in uminode.predecessors:
                UmiDag.visit_predecessors(p, visitor, include_self=True)

    @staticmethod
    def is_terminal(uminode):
        """ Return ``True`` if ``uminode`` has no successors """
        return not uminode.successors

    @staticmethod
    def is_umi_lessthan(umi1, umi2):
        if umi1.pos1 != umi2.pos1:
            # pos1 inverted on purpose to make the earlier UMI the successor in case of equal weights,
            return umi1.pos1 > umi2.pos1
        elif umi1.pos2 != umi2.pos2:
            # pos2 inverted on purpose to be consistent with pos1
            return umi1.pos2 > umi2.pos2
        else:
            # Tag is ordered ascending
            return umi1.tag < umi2.tag

    @staticmethod
    def is_possible_successor(uminode1, uminode2):
        """ Returns True if ``add_successor(uminode1, uminode2)`` would succeed """
        if not (uminode1.weight == uminode2.weight):
            return uminode1.weight < uminode2.weight
        else:
            # Use the total order of UMIs to disambiguate
            return UmiDag.is_umi_lessthan(uminode1.umi, uminode2.umi)

    def __init__(self, similarity_threshold_hamming_distance):
        """ Initialize the DAG for a particular similarity threshold """
        # Parameters
        self.similarity_threshold_hamming_distance = similarity_threshold_hamming_distance
        # Start and end indices of max_hamming_distance+1 tag parts used for indexing
        # To support varying tag lengths, these are computed on-demand for every length we encounter, see slice_tag
        self.tag_slices = dict()
        # All UmiNodes in the DAG, indexed by their UmiFuzzy
        self.nodes = dict()
        # umis[ position1 ][ position2 ][ tag length ][ (i, ith_part_of_tag ) ] contains all umis
        # with position = (position1, position2), and with the i-th part of their tag equal to ith_part_of_tag.
        self.index = SortedDict()
        self.index_pos1_horizon = -inf
        # Set of terminal UmiNodes, i.e. UMIs without a successor, indexed by
        # their pos1max value.
        self.terminals = SortedDict()

    def move_similarity_horizon(self, pos1_horizon):
        """ Moves the horizon that determines how far back ``find_similar_umis'' can look back."""
        if self.index_pos1_horizon > pos1_horizon:
            raise ValueError("horizon cannot move backwards from %d to %d" %
                             (self.index_pos1_horizon, pos1_horizon))
        # Remove index entries with a pos1 left of the horizon
        while self.index and self.index.peekitem(0)[0] < pos1_horizon:
            self.index.popitem(0)
        self.index_pos1_horizon = pos1_horizon

    def find_similar_umis(self, uminode, pos1_deltas=(0, 0), pos2_deltas=(0, 0)):
        """ Returns all UMIs with one or more identical tag slices, and that lie within the specified distances """
        range_pos1 = tuple(uminode.umi.pos1 + d for d in pos1_deltas)
        range_pos2 = tuple(uminode.umi.pos2 + d for d in pos2_deltas)
        if range_pos1[0] < self.index_pos1_horizon:
            raise ValueError("pos1 range [%d, %d] lies (partially) before horizon (%d)" %
                             (range_pos1[0], range_pos1[1], self.index_pos1_horizon))
        already_considered = set()
        for j1 in self.index.irange(minimum=range_pos1[0], maximum=range_pos1[1], inclusive=(True, True)):
            i1 = self.index[j1]
            for j2 in i1.irange(minimum=range_pos2[0], maximum=range_pos2[1], inclusive=(True, True)):
                i2 = i1[j2]
                # The tag length must match perfectly for UMIs to be similar, so we use the slice
                # index corresponding to the UMI tag's length.
                i3 = i2[len(uminode.umi.tag)]
                # Scan all UMIs where at least one slice agrees perfectly, skip nodes we already considered
                for p in self.slice_tag(uminode.umi.tag):
                    if p in i3:
                        # We have an index slot for the current slice p, scan it
                        i4 = i3[p]
                        for uminode2 in i4.difference(already_considered):
                            if uminode2 in already_considered:
                                continue
                            # Must re-check, since the slice-based index is only a preselection
                            if self.similarity_threshold_hamming_distance == 0:
                                # No need to re-check, in that case the slice-based index is precise
                                yield uminode2
                            elif hamming_distance(uminode.umi.tag, uminode2.umi.tag) <= self.similarity_threshold_hamming_distance:
                                # Rechecked, and UMI is actually similar
                                yield uminode2
                        # Other slices may contain the same UMIs again, skip them from now on
                        already_considered.update(i4)

    def add_umi(self, umi, weight, data=None):
        """ Adds the ``UmiFuzzy`` instance with a given ``weight`` to the DAG """
        # Create UmiNode instance, there must be only one for a given UMI
        if umi in self.nodes:
            raise KeyError("UMI %s is already part of the DAG" % (umi,))
        uminode = UmiDag.UmiNode(umi, weight, data)
        self.nodes[umi] = uminode
        # Add to index & terminals (newly created nodes are always terminal)
        self.index_add(uminode)
        self.terminals_add(uminode)
        return uminode

    def add_successor(self, uminode, uminode_succ):
        """ Adds an edge from ``uminode`` towards ``uminode_succ`` """
        # Weights must at least be half-ordered. We then allow edges only
        # The lexicographic order on (weight, umi.pos1, umi.pos2, umi.tag) ensure an acyclic graph
        if ((not uminode.weight < uminode_succ.weight) and
                (not uminode.weight == uminode_succ.weight or not
                self.is_umi_lessthan(uminode.umi, uminode_succ.umi))):
            raise ValueError("edges a -> b only allowed if a.w < b.w or a.w == b.w and w.umi < b.umi")
        # UmiNode will have a successor now, will no longer be terminal
        if self.is_terminal(uminode):
            self.terminals_remove(uminode)
        # Add edge to graph
        uminode.successors.add(uminode_succ)
        uminode_succ.predecessors.add(uminode)

        # Track maximum extent of pos1 within each UMI's set of predecessors
        # by merging UMI's extent into it's new successor's extents.
        def pos1max_merge(s):
            if s.pos1max.value >= uminode.umi.pos1:
                # Stop (and skip all successors of um) if it already encompasses our pos1
                return False
            # Must remove and re-add to terminals if pos1max changes
            if self.is_terminal(s):
                self.terminals_remove(s)
            s.pos1max.value = uminode.umi.pos1
            if self.is_terminal(s):
                self.terminals_add(s)

        self.visit_successors(uminode, pos1max_merge, include_self=False)

    def remove_umi(self, uminode, if_exists=False):
        """ Remove the node from the DAG """
        if (if_exists and not uminode.umi in self.nodes):
            return
        # Remove UmiNode if its indeed part of our DAG
        assert uminode.umi in self.nodes
        del self.nodes[uminode.umi]
        # Remove from terminals list and index
        # Note: The caller might have already pruned the list of
        # terminals and/or the index, must deal gracefully with that
        if self.terminals and self.terminals.peekitem(0)[0] <= uminode.pos1max.value:
            # terminals list not yet pruned
            if self.is_terminal(uminode):
                self.terminals_remove(uminode)
        if self.index and self.index.peekitem(0)[0] <= uminode.umi.pos1:
            # index not yet pruned
            self.index_remove(uminode)

    def find_and_remove_clusters_before(self, pos1_cutoff):
        while self.terminals:
            if self.terminals.peekitem(0)[0] >= pos1_cutoff:
                # Not to the left of the cutoff point, stop
                break
            # To the left of the cutoff point, extract clusters
            # Note that we already prune the terminals list here,
            # remove_umi can deal with that and it's faster than
            # removing everything individually. Note that we *DO*
            # *NOT* prune the index here, but instead prune if
            # ``move_similarity_horizon'' is called. This allows
            # callers to prune the similarity index more aggresively.
            pos1max, terminals = self.terminals.popitem(0)
            for terminal in terminals:
                cluster = UmiCluster(terminal)
                cluster.extract_from_dag(self)
                yield cluster

    def terminals_add(self, uminode):
        """ Must be called if a new terminal node is created, or after changing a node's ``pos1max``"""
        assert self.is_terminal(uminode)
        t = self.terminals.setdefault(uminode.pos1max.value, set())
        t.add(uminode)

    def terminals_remove(self, uminode):
        """ Must be called when removing a node, or before changing a node's ``pos1max``"""
        assert self.is_terminal(uminode)
        t = self.terminals[uminode.pos1max.value]
        assert uminode in t
        t.remove(uminode)
        if not t:
            del self.terminals[uminode.pos1max.value]

    def index_add(self, uminode):
        """ Must be called after creating a new node to add it to the index """
        i1 = self.index.setdefault(uminode.umi.pos1, SortedDict())
        i2 = i1.setdefault(uminode.umi.pos2, dict())
        i3 = i2.setdefault(len(uminode.umi.tag), dict())
        for p in self.slice_tag(uminode.umi.tag):
            i4 = i3.setdefault(p, set())
            i4.add(uminode)

    def index_remove(self, uminode):
        """ Must be called after removing a node to remove it from the index """
        i1 = self.index[uminode.umi.pos1]
        i2 = i1[uminode.umi.pos2]
        i3 = i2[len(uminode.umi.tag)]
        for p in self.slice_tag(uminode.umi.tag):
            i4 = i3[p]
            assert uminode in i4
            i4.remove(uminode)
            # Remove empty set at leaf level
            if not i4:
                del i3[p]
        # Remove empty SortedDicts in intermediate levels
        if not i3:
            del i2[len(uminode.umi.tag)]
        if not i2:
            del i1[uminode.umi.pos2]
        if not i1:
            del self.index[uminode.umi.pos1]

    def slice_tag(self, tag):
        """ Splits the given UMI tag into ``similarity_threshold_hamming_distance+1`` slices """
        tag_length = len(tag)
        if self.similarity_threshold_hamming_distance == 0:
            # Just use the whole tag as the index key
            return [(0, tag)]
        elif tag_length < self.similarity_threshold_hamming_distance + 1:
            # Cannot split tags shorted than hamming_distance + 1 slices characters into hamming_distance + 1 slices,
            # so just use a constant index key, thus effectively disabling the index alltogether.
            return [(0, "")]
        if tag_length not in self.tag_slices:
            # Find a slicing scheme that slices tags of length tag_length into hamming_distance + 1 slices
            tag_slice_nchars = int(floor(tag_length / (self.similarity_threshold_hamming_distance + 1)))
            assert tag_slice_nchars > 0
            self.tag_slices[tag_length] = [
                (i,
                 i * tag_slice_nchars,
                 (i + 1) * tag_slice_nchars if i < self.similarity_threshold_hamming_distance else tag_length)
                for i in range(self.similarity_threshold_hamming_distance + 1)
            ]
        # Split tag
        return [(i, tag[start:end]) for (i, start, end) in self.tag_slices[tag_length]]


class UmiIncrementalErrorCorrector:
    """ Incremental Graph-Based UMI Error Correction

    Finds clusters of similar UMIs within a stream of UMIs sorted  by ``pos1``.
    UMIs and their weights are added incrementally with ``add_umi``. At any time,
    ``complete_clusters`` can be called, and will return an iterator that lists
    complete clusters that haven't yet been reported. A cluster is considered to
    be complete if its right-most member (according to ``pos1``) is too far left
    to possibly be considered similar to any UMI not yet added.

    **Usage**::

        ec = UmiIncrementalErrorCorrector(...)
        for umi, weight in umis_sorted_by_pos1:
            node = ec.add_umi(umi, weight, additional_data)
            for cluster in ec.complete_clusters():
                output(cluster)
        for cluster in ec.complete_clusters(flush=True):
            output(cluster)
    """

    def __init__(self, max_weight_relative=1.0, min_weight_difference=0,
                 max_1st_mapping_distance=0, max_2nd_mapping_distance=0, max_hamming_distance=0):
        # Parameters
        self.dag = UmiDag(max_hamming_distance)
        self.max_1st_mapping_distance = max_1st_mapping_distance
        self.max_2nd_mapping_distance = max_2nd_mapping_distance
        self.max_hamming_distance = max_hamming_distance
        self.max_weight_relative = max_weight_relative
        self.min_weight_difference = min_weight_difference
        # Track current position, must always increase
        self.position = None

    def add_umi(self, umi, weight, data=None):
        # Track the current position, which must never decrease
        assert (self.position is None) or (umi.pos1 >= self.position)
        self.position = umi.pos1
        # Inform the DAG that how far back we still need to search for similar UMIs
        self.dag.move_similarity_horizon(umi.pos1 - self.max_1st_mapping_distance)
        # Add UMI to graph, add edges to existing similar UMIs, process & remove complete clusters
        uminode = self.dag.add_umi(umi, weight, data)
        self.add_edges_to_similar_existing_umis(uminode)

    def complete_clusters(self, flush=False):
        # Now extract any clusters, i.e. groups of UMIs that are predecessors of
        # a single node (the cluster center) and all far enough to the left of
        # the current position that they cannot be "similar" to any UMI yet to come
        if flush:
            cutoff = self.position + 1
        else:
            cutoff = self.position - self.max_1st_mapping_distance
        for cluster in self.dag.find_and_remove_clusters_before(cutoff):
            yield cluster

    def is_mergable(self, uminode, uminode_succ):
        # Should only be guaranteed by the pre-filtering in UmiDag's find_similar_umis
        assert abs(uminode.umi.pos1 - uminode_succ.umi.pos1) <= self.max_1st_mapping_distance
        assert abs(uminode.umi.pos2 - uminode_succ.umi.pos2) <= self.max_2nd_mapping_distance
        assert hamming_distance(uminode.umi.tag, uminode_succ.umi.tag) <= self.max_hamming_distance
        return uminode.weight <= uminode_succ.weight * self.max_weight_relative and \
               uminode.weight <= uminode_succ.weight - self.min_weight_difference

    def add_edges_to_similar_existing_umis(self, uminode):
        # Add similar umis as successors (if they have a higher weight), other as a predecessor
        # Since we're processing UMIs in ascending order of position1, no UMI can have a larger pos1
        pos1_delta = (-self.max_1st_mapping_distance, 0)
        # For the position2 range, however, we must use a symmetric interval
        pos2_delta = (-self.max_2nd_mapping_distance, self.max_2nd_mapping_distance)
        for uminode2 in self.dag.find_similar_umis(uminode, pos1_delta, pos2_delta):
            if uminode2 is uminode:
                continue
            if self.dag.is_possible_successor(uminode, uminode2):
                uminode_pred, uminode_succ = uminode, uminode2
            else:
                uminode_pred, uminode_succ = uminode2, uminode
            if self.is_mergable(uminode_pred, uminode_succ):
                self.dag.add_successor(uminode_pred, uminode_succ)
