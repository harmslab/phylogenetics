import os
from nose.tools import assert_equal

from .base import RWTestCase
from phylogenetics.homologs import Homolog

class TestWrite(RWTestCase):
    """Test the Write object of HomologSet object

    Anytime a file is written in the test, make sure to add its path to
    the `fnames` attribute of the test class!

    """
    def test_fasta(self):
        """ Test fasta writing.

            Not sure this is the best way to test this method.
        """
        # Check that at least this one sequence is in the full fasta
        known = """>XX00000001\nVVVVAAAAAAWWWMWMWMWWAVVVVV"""
        self.assertIn(known,self.HomologSet.Write.fasta())

    def test_fasta_fname(self):
        """ test that writing to a file works. """
        fasta_fname = "test_fasta_fname.fasta"
        fasta_fname_path = os.path.join(os.getcwd(), fasta_fname)
        self.fnames.append(fasta_fname_path)

        fasta_string = self.HomologSet.Write.fasta()
        self.HomologSet.Write.fasta(fname=fasta_fname_path)

        # Check that the output file was created.
        check = os.path.isfile(fasta_fname_path)
        self.assertTrue(check)

        # If the file exists, check that input
        if check:
            with open(fasta_fname_path, "r") as f:
                file_output = f.read()
            self.assertEqual(fasta_string, file_output)


class TestRead(RWTestCase):
    """Test the Read object of HomologSet object

    Anytime a file is written in the test, make sure to add its path to
    the `fnames` attribute of the test class!

    """
    def test_fasta(self):
        """ read a fasta file. """
        fasta_string = """>XX00000001\nVVVVWWWMWMMMMMMMWMWWAVVVVV"""
        sequence = "VVVVWWWMWMMMMMMMWMWWAVVVVV"
        self.HomologSet.Read.fasta(fasta_string)
        self.assertEqual(sequence, self.HomologSet.XX00000001.sequence)
