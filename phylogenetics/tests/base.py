import unittest


from phylogenetics import homologs

class BaseTestCase(unittest.TestCase):

    def setUp(self):
        """Construct a HomologSet object for use elsewhere.
        """
        h1 = homologs.Homolog("XX00000000", species="human", sequence="AAAAAVVAXXAAAAMMAAAAWWWMMVVVAAAAAA")
        h2 = homologs.Homolog("XX00000001", species="human", sequence="VVVVAAAAAAWWWMWMWMWWAVVVVV")
        h3 = homologs.Homolog("XX00000002", species="human", sequence="VVAAAAAAXXPPMMMAAAAAAAAVVVVAV")
        h4 = homologs.Homolog("XX00000003", species="human", sequence="VVAAVPPPMAAAPPMMMAAAAAAAAVVVVAV")
        h5 = homologs.Homolog("XX00000004", species="human", sequence="VVAAVPPPMAAAAMMAAAAAVVVVAV")
        h6 = homologs.Homolog("XX00000005", species="human", sequence="AAVVVAAAAPMMMMMMAAAPPMMMAAAA")
        h7 = homologs.Homolog("XX00000006", species="human", sequence="AAVVVAXXXXXXAAAAAPPMMMAAAA")
        print("Setting up HomologSet object.")
        self.HomologSet = homologs.HomologSet([h1,h2,h3,h4,h5,h6,h7])

    def tearDown(self):
        """Dispose of the homolog-set object.
        """
        print("Tearing down HomologSet object.")
        delattr(self, "HomologSet")
