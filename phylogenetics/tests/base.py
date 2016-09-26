import unittest


from phylogenetics import sequences

class BaseTestCase(unittest.TestCase):

    def setUp(self):
        """Construct a SequenceList object for use elsewhere.
        """
        h1 = sequences.Sequence("XX00000000", species="human", sequence="AAAAAVVAXXAAAAMMAAAAWWWMMVVVAAAAAA")
        h2 = sequences.Sequence("XX00000001", species="human", sequence="VVVVAAAAAAWWWMWMWMWWAVVVVV")
        h3 = sequences.Sequence("XX00000002", species="human", sequence="VVAAAAAAXXPPMMMAAAAAAAAVVVVAV")
        h4 = sequences.Sequence("XX00000003", species="human", sequence="VVAAVPPPMAAAPPMMMAAAAAAAAVVVVAV")
        h5 = sequences.Sequence("XX00000004", species="human", sequence="VVAAVPPPMAAAAMMAAAAAVVVVAV")
        h6 = sequences.Sequence("XX00000005", species="human", sequence="AAVVVAAAAPMMMMMMAAAPPMMMAAAA")
        h7 = sequences.Sequence("XX00000006", species="human", sequence="AAVVVAXXXXXXAAAAAPPMMMAAAA")
        print("Setting up HomologSet object.")
        self.SequenceList = sequences.SequenceList([h1,h2,h3,h4,h5,h6,h7])

    def tearDown(self):
        """Dispose of the SequenceList object.
        """
        print("Tearing down SequenceList object.")
        delattr(self, "SequenceList")
