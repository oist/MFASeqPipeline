from unittest import TestCase
from pyngstuff.umi_error_correction_incremental import UmiDag
from pyngstuff import UmiIncrementalErrorCorrector, UmiCluster, UmiFuzzy

class TestUmiDag(TestCase):
    def test_create_dag(self):
        dag = UmiDag(similarity_threshold_hamming_distance=1)
        self.assertDagEmpty(dag)
        return dag

    def test_add_umi(self, dag=None, umi=UmiFuzzy(pos1=1, tag="ACT"), weight=1):
        if dag is None:
            dag = self.test_create_dag()
        node = dag.add_umi(umi, weight=weight)
        self.assertEqual(node.umi, umi)
        self.assertEqual(node.weight, weight)
        self.assertEqual(node.pos1max.value, umi.pos1)
        self.assertFalse(node.successors)
        self.assertFalse(node.predecessors)
        self.assertTrue(dag.is_terminal(node))
        self.assertNodeInDag(node, dag)
        return node

    def test_add_umi_error_dupkey(self):
        dag = self.test_create_dag()
        self.test_add_umi(dag=dag)
        self.assertRaises(KeyError, lambda: self.test_add_umi(dag=dag))

    def test_remove_umi(self):
        dag = self.test_create_dag()
        uminode1 = self.test_add_umi(dag)
        uminode2 = self.test_add_umi(dag, umi=UmiFuzzy(pos1=1, tag="ATT"))
        dag.remove_umi(uminode1)
        self.assertNodeInDag(uminode2, dag)
        dag.remove_umi(uminode2)
        self.assertDagEmpty(dag)

    def test_add_successor(self):
        # ACT  <-----|            w=2
        #            a
        #            |
        # AAT --b-> ATT           w=1
        #  ^
        #  |
        #  |---c--- AAT           w=0
        dag = self.test_create_dag()
        uminode1 = self.test_add_umi(dag, umi=UmiFuzzy(pos1=1, tag="ACT"), weight=2)
        uminode2 = self.test_add_umi(dag, umi=UmiFuzzy(pos1=2, tag="ATT"), weight=1)
        dag.add_successor(uminode2, uminode1) # (a)
        self.assertTrue(dag.is_terminal(uminode1))
        self.assertFalse(dag.is_terminal(uminode2))
        self.assertEqual(uminode1.pos1max.value, 2)
        self.assertIn(uminode1, dag.terminals[uminode1.pos1max.value])
        self.assertNotIn(1, dag.terminals)
        uminode3 = self.test_add_umi(dag, umi=UmiFuzzy(pos1=2, tag="AAT"), weight=1)
        dag.add_successor(uminode3, uminode2) # (b)
        self.assertTrue(dag.is_terminal(uminode1))
        self.assertFalse(dag.is_terminal(uminode2))
        self.assertFalse(dag.is_terminal(uminode3))
        self.assertEqual(uminode1.pos1max.value, 2)
        self.assertIn(uminode1, dag.terminals[uminode1.pos1max.value])
        self.assertNotIn(1, dag.terminals)
        uminode4 = self.test_add_umi(dag, umi=UmiFuzzy(pos1=3, tag="AAT"), weight=0)
        dag.add_successor(uminode4, uminode3) # (c)
        self.assertTrue(dag.is_terminal(uminode1))
        self.assertFalse(dag.is_terminal(uminode2))
        self.assertFalse(dag.is_terminal(uminode3))
        self.assertFalse(dag.is_terminal(uminode4))
        self.assertEqual(uminode1.pos1max.value, 3)
        self.assertIn(uminode1, dag.terminals[uminode1.pos1max.value])
        self.assertNotIn(1, dag.terminals)
        self.assertNotIn(2, dag.terminals)

    def test_add_successor_error_wrongdir(self):
        dag = self.test_create_dag()
        uminode1 = self.test_add_umi(dag, umi=UmiFuzzy(pos1=1, pos2=0, tag="ACT"), weight=2)
        uminode2 = self.test_add_umi(dag, umi=UmiFuzzy(pos1=2, pos2=0, tag="ATT"), weight=1)
        self.assertRaises(ValueError, lambda: dag.add_successor(uminode1, uminode2))

    def test_find_similar_umis(self):
        dag = self.test_create_dag()
        uminode1 = self.test_add_umi(dag, umi=UmiFuzzy(pos1=1, pos2=0, tag="ACT"), weight=1)
        uminode2 = self.test_add_umi(dag, umi=UmiFuzzy(pos1=2, pos2=0, tag="ATT"), weight=2)
        uminode3 = self.test_add_umi(dag, umi=UmiFuzzy(pos1=3, pos2=0, tag="GTT"), weight=1)
        uminode4 = self.test_add_umi(dag, umi=UmiFuzzy(pos1=2, pos2=0, tag="TTA"), weight=3)
        uminode5 = self.test_add_umi(dag, umi=UmiFuzzy(pos1=2, pos2=1, tag="TTA"), weight=3)
        self.assertSetEqual(set(dag.find_similar_umis(uminode1, (-1,1))),
                            {uminode1, uminode2})
        self.assertSetEqual(set(dag.find_similar_umis(uminode1, (-1,0))),
                            {uminode1})
        self.assertSetEqual(set(dag.find_similar_umis(uminode2, (-1,1))),
                            {uminode1, uminode2, uminode3})
        self.assertSetEqual(set(dag.find_similar_umis(uminode2, (-1,0))),
                            {uminode1, uminode2})
        self.assertSetEqual(set(dag.find_similar_umis(uminode2, ( 0,1))),
                            {uminode2, uminode3})
        self.assertSetEqual(set(dag.find_similar_umis(uminode3, (-1,1))),
                            {uminode2, uminode3})
        self.assertSetEqual(set(dag.find_similar_umis(uminode4, (-1,1))),
                            {uminode4})
        self.assertSetEqual(set(dag.find_similar_umis(uminode4, (-1,1), (0,1))),
                            {uminode4, uminode5})

    def test_find_and_remove_clusters_before(self):
        dag = self.test_create_dag()
        umi1 = (UmiFuzzy(pos1=1, pos2=0, tag="ACT"), 1)
        uminode1 = self.test_add_umi(dag, umi=umi1[0], weight=umi1[1])
        umi2 = (UmiFuzzy(pos1=2, pos2=0, tag="ATT"), 2)
        uminode2 = self.test_add_umi(dag, umi=umi2[0], weight=umi2[1])
        dag.add_successor(uminode1, uminode2)
        umi3 = (UmiFuzzy(pos1=3, pos2=0, tag="GTT"), 1)
        uminode3 = self.test_add_umi(dag, umi=umi3[0], weight=umi3[1])
        dag.add_successor(uminode3, uminode2)
        umi4 = (UmiFuzzy(pos1=2, pos2=0, tag="TTA"), 3)
        uminode4 = self.test_add_umi(dag, umi=umi4[0], weight=umi4[1])
        umi5 = (UmiFuzzy(pos1=2, pos2=1, tag="TTA"), 2)
        uminode5 = self.test_add_umi(dag, umi=umi5[0], weight=umi5[1])
        dag.add_successor(uminode5, uminode4)
        umi6 = (UmiFuzzy(pos1=0, pos2=0, tag="NNN"), 0)
        uminode6 = self.test_add_umi(dag, umi=umi6[0], weight=umi6[1])
        dag.add_successor(uminode6, uminode1)
        dag.add_successor(uminode6, uminode5)
        umi7 = (UmiFuzzy(pos1=1, pos2=0, tag="NNN"), 0)
        uminode7 = self.test_add_umi(dag, umi=umi7[0], weight=umi7[1])
        dag.add_successor(uminode7, uminode6)
        dag.add_successor(uminode7, uminode6)
        umi8 = (UmiFuzzy(pos1=2, pos2=0, tag="NNN"), 0)
        uminode8 = self.test_add_umi(dag, umi=umi8[0], weight=umi8[1])
        dag.add_successor(uminode8, uminode7)
        dag.add_successor(uminode8, uminode3)
        clusters_before_2 = set(dag.find_and_remove_clusters_before(2))
        clusters_before_3 = set(dag.find_and_remove_clusters_before(3))
        clusters_before_4 = set(dag.find_and_remove_clusters_before(4))
        def cluster_to_triple(cluster):
            return frozenset((c.principal_umi,
                              frozenset((u,w) for u,w,_ in c.unambiguous()),
                              frozenset((u,w) for u,w,_ in c.ambiguous()))
                             for c in cluster)
        self.assertSetEqual(cluster_to_triple(clusters_before_2),
                            set())
        self.assertSetEqual(cluster_to_triple(clusters_before_3),
                            { ( umi4[0], frozenset({umi4,umi5}), frozenset({umi6, umi7, umi8}) ) })
        self.assertSetEqual(cluster_to_triple(clusters_before_4),
                            { ( umi2[0], frozenset({umi1,umi2,umi3}), frozenset({}) ) })
        self.assertDagEmpty(dag)


    def test_example(self):
        #      A
        #    / | \
        #  C   B  \
        #       \  \
        #        -- D
        #           |
        #           E
        a = (UmiFuzzy(pos1=1168, pos2=-1, tag='TTGGGGTCTG'), 8)
        b = (UmiFuzzy(pos1=1168, pos2=-1, tag='TTCGGGTCTG'), 8)
        c = (UmiFuzzy(pos1=1168, pos2=-1, tag='TTGGTGTCTG'), 7)
        d = (UmiFuzzy(pos1=1169, pos2=-1, tag='TTTGGGTCTG'), 4)
        e = (UmiFuzzy(pos1=1168, pos2=-1, tag='TATGGGTCTG'), 1)
        dag = UmiDag(similarity_threshold_hamming_distance=1)
        # Add UMIs and similarity edges
        self.assertDagEmpty(dag)
        self.test_add_umi(dag, a[0], a[1])
        self.test_add_umi(dag, b[0], b[1])
        dag.add_successor(dag.nodes[b[0]], dag.nodes[a[0]])
        self.test_add_umi(dag, c[0], c[1])
        dag.add_successor(dag.nodes[c[0]], dag.nodes[a[0]])
        self.test_add_umi(dag, d[0], d[1])
        dag.add_successor(dag.nodes[d[0]], dag.nodes[a[0]])
        dag.add_successor(dag.nodes[d[0]], dag.nodes[b[0]])
        self.test_add_umi(dag, e[0], e[1])
        dag.add_successor(dag.nodes[e[0]], dag.nodes[d[0]])
        self.assertEqual(dag.nodes[a[0]].pos1max.value, 1169)
        self.assertSetEqual(dag.terminals[1169], {dag.nodes[a[0]]})
        self.assertSetEqual(dag.nodes[a[0]].successors, set())
        self.assertSetEqual(dag.nodes[b[0]].successors, {dag.nodes[a[0]]})
        self.assertSetEqual(dag.nodes[c[0]].successors, {dag.nodes[a[0]]})
        self.assertSetEqual(dag.nodes[d[0]].successors, {dag.nodes[a[0]], dag.nodes[b[0]]})
        self.assertSetEqual(dag.nodes[e[0]].successors, {dag.nodes[d[0]]})
        # Now extract cluster
        cluster = UmiCluster(dag.nodes[a[0]])
        cluster.extract_from_dag(dag)
        self.assertDictEqual(dag.nodes, {})


    def assertDagEmpty(self, dag):
        self.assertFalse(dag.nodes)
        self.assertFalse(dag.index)
        self.assertFalse(dag.terminals)

    def assertNodeInDag(self, node, dag):
        self.assertIs(dag.nodes[node.umi], node)
        self.assertIn(node, dag.terminals[node.umi.pos1])


class TestUmiIncrementalErrorCorrector(TestCase):
    def test_ambiguity(self):
        #    a                 | w=5
        # b     c        i     | w=4
        # d        f        k  | w=3
        #    e     g     m     | w=2
        #       h     n        | w=1
        # ---------------------|
        # 0  1  2  3  4  5  6
        #
        # Nodes are designed to be similar to the node directly
        # above them (however far), and their immediate upper-
        # diagonal neighbours, except for n which is *not* similar
        # to g, only to m
        #
        # Weights corresponds to the tier of the node.
        #
        # The expected clusters are
        #    Center a, Members a,b,c,d,e
        #    Center i, Members i,k,m,n
        #
        # Nodes f,g,h are ambiguous, could be member of either cluster
        umis = {
            'a': (UmiFuzzy(pos1=1, pos2=1, tag='ACGT'), 5),
            'b': (UmiFuzzy(pos1=0, pos2=0, tag='ACGA'), 4),
            'c': (UmiFuzzy(pos1=2, pos2=2, tag='ACAT'), 4),
            'd': (UmiFuzzy(pos1=0, pos2=1, tag='AAGA'), 3),
            'e': (UmiFuzzy(pos1=1, pos2=2, tag='AAGT'), 2),
            'f': (UmiFuzzy(pos1=3, pos2=3, tag='CCAT'), 3),
            'g': (UmiFuzzy(pos1=3, pos2=3, tag='CCGT'), 2),
            'h': (UmiFuzzy(pos1=2, pos2=3, tag='ACGT'), 1),
            'i': (UmiFuzzy(pos1=4, pos2=4, tag='CGAT'), 4),
            'k': (UmiFuzzy(pos1=5, pos2=5, tag='CGAG'), 3),
            'm': (UmiFuzzy(pos1=4, pos2=4, tag='CGGG'), 2),
            'n': (UmiFuzzy(pos1=3, pos2=4, tag='GGGG'), 1)
        }
        successor_labels = {
            'a': {},
            'b': {'a'},
            'c': {'a'},
            'd': {'b'},
            'e': {'a', 'd'},
            'f': {'c', 'i'},
            'g': {'f'},
            'h': {'c', 'e', 'g'}, # should not be there!
            'i': {},
            'k': {'i'},
            'm': {'k'},
            'n': {'m'}
        }
        member_labels = {
            'a': {'a', 'b', 'c', 'd', 'e'},
            'i': {'i', 'k', 'm', 'n'}
        }
        member_ambiguous_labels = {
            'a': {'f', 'g', 'h'},
            'i': {}
        }

        # Translate from labels to UMIs
        umi_weight = { u:w for u,w in umis.values() }
        successor_umis = { umis[lp][0]: set(umis[ls][0] for ls in s)
                       for lp, s in successor_labels.items() }
        member_umis = { umis[lc][0]: set(umis[lm][0] for lm in m)
                       for lc, m in member_labels.items() }
        member_ambiguous_umis = { umis[lc][0]: set(umis[lm][0] for lm in m)
                                  for lc, m in member_ambiguous_labels.items() }

        # Create error corrector & add UMIs
        ec = UmiIncrementalErrorCorrector(max_hamming_distance=1,
                                          max_1st_mapping_distance=2, max_2nd_mapping_distance=1)
        for label, (umi, weight) in sorted(umis.items(), key=lambda i: i[1][0].pos1):
            ec.add_umi(umi, weight)

        # Check that we ended up with the expected DAG structure
        # TODO: Actually, the edge h -> a should not be there according to the diagram above!
        for label, (umi, weight) in umis.items():
            self.assertSetEqual(ec.dag.nodes[umi].successors,
                                set(ec.dag.nodes[s] for s in successor_umis[umi]),
                                "Successors of node '%s' (%s)" % (label, umi))

        # Check that the error corrector produced the expected clusters
        complete_clusters = set(ec.complete_clusters(flush=True))
        self.assertEqual(len(complete_clusters), len(member_umis))
        for cluster in complete_clusters:
            if cluster.principal_umi not in member_umis:
                self.assertTrue(False, 'Unexpected terminal UMI %s' % (cluster.principal_umi,))
            else:
                self.assertSetEqual(set(umi for umi, _, _ in cluster.unambiguous()),
                                    member_umis[cluster.principal_umi])
                self.assertSetEqual(set(umi for umi, _, _ in cluster.ambiguous()),
                                    member_ambiguous_umis[cluster.principal_umi])
                self.assertEqual(cluster.total_weight,
                                 sum(umi_weight[u] for u in member_umis[cluster.principal_umi]))


    def test_cluster_completion(self):
        # Create error corrector
        ec = UmiIncrementalErrorCorrector(max_hamming_distance=1,
                                          max_1st_mapping_distance=2)

        # Function to extract complete clusters and store in umis_corrected, umis_versions
        umis_corrected = set()
        umis_versions = dict()
        def handle_complete_clusters(flush=False):
            for cluster in ec.complete_clusters(flush):
                umis_corrected.add(cluster.principal_umi)
                umis_versions[cluster.principal_umi] = set((u,w) for u,w,_ in cluster.unambiguous())

        # Add nodes, check that complete clusters are produced at the expected time points
        ec.add_umi(UmiFuzzy(pos1=1, tag='AAA'), weight=5)
        handle_complete_clusters()
        self.assertSetEqual(umis_corrected, set())
        ec.add_umi(UmiFuzzy(pos1=1, tag='ACA'), weight=3)
        handle_complete_clusters()
        self.assertSetEqual(umis_corrected, set())
        ec.add_umi(UmiFuzzy(pos1=1, tag='ATT'), weight=1)
        handle_complete_clusters()
        self.assertSetEqual(umis_corrected, set())
        ec.add_umi(UmiFuzzy(pos1=2, tag='ACT'), weight=2)
        handle_complete_clusters()
        self.assertSetEqual(umis_corrected, set())
        ec.add_umi(UmiFuzzy(pos1=4, tag='AAT'), weight=2)
        handle_complete_clusters()
        self.assertSetEqual(umis_corrected, set())
        ec.add_umi(UmiFuzzy(pos1=6, tag='CCC'), weight=1)
        handle_complete_clusters()
        self.assertSetEqual(umis_corrected, set())
        ec.add_umi(UmiFuzzy(pos1=7, tag='CCC'), weight=1)
        handle_complete_clusters()
        self.assertSetEqual(umis_corrected, {UmiFuzzy(pos1=1, tag='AAA')})
        self.assertSetEqual(umis_versions[UmiFuzzy(pos1=1, tag='AAA')],
                            {(UmiFuzzy(pos1=1, tag='AAA'), 5), (UmiFuzzy(pos1=1, tag='ACA'), 3),
                             (UmiFuzzy(pos1=1, tag='ATT'), 1), (UmiFuzzy(pos1=2, tag='ACT'), 2),
                             (UmiFuzzy(pos1=4, tag='AAT'), 2)})
        handle_complete_clusters(flush=True)
        self.assertSetEqual(umis_corrected, {UmiFuzzy(pos1=1, tag='AAA'), UmiFuzzy(pos1=6, tag='CCC')})
        self.assertSetEqual(umis_versions[UmiFuzzy(pos1=6, tag='CCC')],
                            {(UmiFuzzy(pos1=6, tag='CCC'), 1), (UmiFuzzy(pos1=7, tag='CCC'), 1)})

