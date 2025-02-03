from unittest import TestCase
from collections import namedtuple
from pyngstuff import UmiTagSource, UmiTagSourceReadname, UmiTagSourceBAMTags

class ReadMockup(namedtuple("ReadMockup", ["query_name", "tags"])):
    def get_tag(self, tag):
        return self.tags[tag]

class TestUmiTagSource(TestCase):
    def test_readname_nocellbc(self):
        s = UmiTagSource.parse("readname:U")
        self.assertIsInstance(s, UmiTagSourceReadname)
        self.assertEqual("readname:U", str(s))
        self.assertEqual(True, s.has_umitag)
        self.assertEqual(False, s.has_cellbc)
        tags = s(ReadMockup("myread:AGACTG", {}))
        self.assertEqual("AGACTG", tags[0])
        self.assertIsNone(tags[1])
        with self.assertRaises(RuntimeError) as cm:
            s(ReadMockup("myread", {}))

        s = UmiTagSource.parse("readname_U")
        self.assertIsInstance(s, UmiTagSourceReadname)
        self.assertEqual("readname_U", str(s))
        self.assertEqual(True, s.has_umitag)
        self.assertEqual(False, s.has_cellbc)
        tags = s(ReadMockup("myread_GGAGCTA", {}))
        self.assertEqual("GGAGCTA", tags[0])
        self.assertIsNone(tags[1])
        with self.assertRaises(RuntimeError) as cm:
            s(ReadMockup("myread:GGAGCTA", {}))

        s = UmiTagSource.parse("readname-U")
        self.assertIsInstance(s, UmiTagSourceReadname)
        self.assertEqual("readname-U", str(s))
        self.assertEqual(True, s.has_umitag)
        self.assertEqual(False, s.has_cellbc)
        tags = s(ReadMockup("myread-GGAGCTA", {}))
        self.assertEqual("GGAGCTA", tags[0])
        self.assertIsNone(tags[1])
        tags = s(ReadMockup("myread-GGAGCTA-AGACGT", {}))
        self.assertEqual("GGAGCTA-AGACGT", tags[0])
        self.assertIsNone(tags[1])
        with self.assertRaises(RuntimeError) as cm:
            s(ReadMockup("myread:GGAGCTA", {}))

    def test_readname_cellbc(self):
        s = UmiTagSource.parse("readname:U:C")
        self.assertIsInstance(s, UmiTagSourceReadname)
        self.assertEqual("readname:U:C", str(s))
        self.assertEqual(True, s.has_umitag)
        self.assertEqual(True, s.has_cellbc)
        tags = s(ReadMockup("myread:AGACTG:GAGTAC", {}))
        self.assertEqual("AGACTG", tags[0])
        self.assertEqual("GAGTAC", tags[1])
        with self.assertRaises(RuntimeError) as cm:
            s(ReadMockup("myread", {}))

        s = UmiTagSource.parse("readname_C:U")
        self.assertIsInstance(s, UmiTagSourceReadname)
        self.assertEqual("readname_C:U", str(s))
        self.assertEqual(True, s.has_umitag)
        self.assertEqual(True, s.has_cellbc)
        tags = s(ReadMockup("myread_GGAGCTA:ATACAGT", {}))
        self.assertEqual("ATACAGT", tags[0])
        self.assertEqual("GGAGCTA", tags[1])
        with self.assertRaises(RuntimeError) as cm:
            s(ReadMockup("myread_GGAGCTA", {}))

    def test_bamtag_nocellbc(self):
        s = UmiTagSource.parse("U=RX")
        self.assertIsInstance(s, UmiTagSourceBAMTags)
        self.assertEqual("RX", s.bamtag_umitag)
        self.assertEqual(None, s.bamtag_cellbc)
        tags = s(ReadMockup("myread", {"RX": "GGAGCTA"}))
        self.assertEqual("GGAGCTA", tags[0])
        self.assertIsNone(tags[1])

        s = UmiTagSource.parse("U=OX")
        self.assertIsInstance(s, UmiTagSourceBAMTags)
        self.assertEqual("OX", s.bamtag_umitag)
        self.assertEqual(None, s.bamtag_cellbc)
        tags = s(ReadMockup("myread", {"OX": "GATACCT"}))
        self.assertEqual("GATACCT", tags[0])
        self.assertIsNone(tags[1])

    def test_bamtag_cellbc(self):
        s = UmiTagSource.parse("U=RX,C=BC")
        self.assertIsInstance(s, UmiTagSourceBAMTags)
        self.assertEqual("RX", s.bamtag_umitag)
        self.assertEqual("BC", s.bamtag_cellbc)
        tags = s(ReadMockup("myread", {"RX": "GGAGCTA", "BC": "ATACCGT"}))
        self.assertEqual("GGAGCTA", tags[0])
        self.assertEqual("ATACCGT", tags[1])

        s = UmiTagSource.parse("C=CR,U=OX")
        self.assertIsInstance(s, UmiTagSourceBAMTags)
        self.assertEqual("OX", s.bamtag_umitag)
        self.assertEqual("CR", s.bamtag_cellbc)
        tags = s(ReadMockup("myread", {"OX": "GATACCT", "CR": "GGATACG"}))
        self.assertEqual("GATACCT", tags[0])
        self.assertEqual("GGATACG", tags[1])

    def test_autodetect_nocellbc(self):
        r = ReadMockup("ILLUMINA_0179:6:1101:10714:18881:ACTGCAAAGTCTAGCAGTTG-TTGGAAACTACGTTCTCGAC", {})
        s = UmiTagSource.autodetect(r, cellbarcode=False)
        self.assertIsInstance(s, UmiTagSourceReadname)
        self.assertEqual("readname:U", str(s))
        tags = s(r)
        self.assertEqual("ACTGCAAAGTCTAGCAGTTG-TTGGAAACTACGTTCTCGAC", tags[0])
        self.assertIsNone(tags[1])

        r = ReadMockup("ILLUMINA_0179:6:1101:10714:18881_ACTGCAAAGTCTAGCAGTTG-TTGGAAACTACGTTCTCGAC", {})
        s = UmiTagSource.autodetect(r, cellbarcode=False)
        self.assertIsInstance(s, UmiTagSourceReadname)
        self.assertEqual("readname_U", str(s))
        tags = s(r)
        self.assertEqual("ACTGCAAAGTCTAGCAGTTG-TTGGAAACTACGTTCTCGAC", tags[0])
        self.assertIsNone(tags[1])

        r = ReadMockup("ILLUMINA_0179:6:1101:10714:18881_ACTGCAAAGTCTAGCAGTTG-TTGGAAACTACGTTCTCGAC", {"RX": "ACTGCAAAGTCTAGCAGTTG-TTGGAAACTACGTTCTCGAC"})
        s = UmiTagSource.autodetect(r, cellbarcode=False)
        self.assertIsInstance(s, UmiTagSourceBAMTags)
        self.assertEqual("U=RX", str(s))
        tags = s(r)
        self.assertEqual("ACTGCAAAGTCTAGCAGTTG-TTGGAAACTACGTTCTCGAC", tags[0])
        self.assertIsNone(tags[1])

    def test_autodetect_cellbc(self):
        r = ReadMockup("ILLUMINA_0179:6:1101:10714:18881:AACGTACGTA:ACTGCAAAGTCTAGCAGTTG-TTGGAAACTACGTTCTCGAC", {})
        s = UmiTagSource.autodetect(r, cellbarcode=True)
        self.assertIsInstance(s, UmiTagSourceReadname)
        self.assertEqual("readname:C:U", str(s))
        tags = s(r)
        self.assertEqual("ACTGCAAAGTCTAGCAGTTG-TTGGAAACTACGTTCTCGAC", tags[0])
        self.assertEqual("AACGTACGTA", tags[1])

        r = ReadMockup("ILLUMINA_0179:6:1101:10714:18881_AACGTACGTA_ACTGCAAAGTCTAGCAGTTG-TTGGAAACTACGTTCTCGAC", {})
        s = UmiTagSource.autodetect(r, cellbarcode=True)
        self.assertIsInstance(s, UmiTagSourceReadname)
        self.assertEqual("readname_C_U", str(s))
        tags = s(r)
        self.assertEqual("ACTGCAAAGTCTAGCAGTTG-TTGGAAACTACGTTCTCGAC", tags[0])
        self.assertEqual("AACGTACGTA", tags[1])

        r = ReadMockup("ILLUMINA_0179:6:1101:10714:18881_ACTGCAAAGTCTAGCAGTTG-TTGGAAACTACGTTCTCGAC",
                       {"RX": "ACTGCAAAGTCTAGCAGTTG-TTGGAAACTACGTTCTCGAC", "CB": "AGCCGATA"})
        s = UmiTagSource.autodetect(r, cellbarcode=True)
        self.assertIsInstance(s, UmiTagSourceBAMTags)
        self.assertEqual("U=RX,C=CB", str(s))
        tags = s(r)
        self.assertEqual("ACTGCAAAGTCTAGCAGTTG-TTGGAAACTACGTTCTCGAC", tags[0])
        self.assertEqual("AGCCGATA", tags[1])
