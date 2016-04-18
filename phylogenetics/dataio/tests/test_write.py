from nose.tools import assert_equal

from ..homologio import Write
from phylogenetics.homologs import Homolog

class TestWrite:

    def setup(self):
        self.Homolog = Homolog("XX00000001",
            sequence="ASDGASHAERRHDGASDASE",
            organism="homosapien",
            defline="human"
        )

        self.Write(self.Homolog)

    def teardown(self):
        self.Homolog = None
