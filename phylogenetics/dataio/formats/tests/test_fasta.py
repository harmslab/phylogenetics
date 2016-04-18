# Test the fasta input/output modules

# Test imports
from nose import with_setup
from nose.tools import assert_equal

# Module imports
from ..fasta import read, write

# Test dataset 1
data1 = """>XX00000001\nASEFASVSADFAWEFSDFASDFAWEVGSDHDTWAETASFSFRAWEF
>XX00000002\nASFASFAEHSDHJADFAEGEFGASDFAWEFAFAWEF"""
sequences1 = [
    (("XX00000001",), "ASEFASVSADFAWEFSDFASDFAWEVGSDHDTWAETASFSFRAWEF"),
    (("XX00000002",), "ASFASFAEHSDHJADFAEGEFGASDFAWEFAFAWEF"),
]

# Test dataset 2
data2 = """>XX00000001|this is cool|WOO\nASDFASFAWEFASGSDAGASFASDFASFASGASDG
>XX00000002|testing|AWESOME\nASDFASFAWEFASGSDAGASFASDFASFASGASDG"""

sequences2 = [
    (("XX00000001","this is cool", "WOO",), "ASDFASFAWEFASGSDAGASFASDFASFASGASDG"),
    (("XX00000002","testing","AWESOME",), "ASDFASFAWEFASGSDAGASFASDFASFASGASDG"),
]

# Test reading only single sequence fastas and sequence data
data3 = """>XX00000003\nAGAEWGAHDSGAWSCDCADSGADAG\
"""

sequences3 = (("XX00000003",), "AGAEWGAHDSGAWSCDCADSGADAG")

def test_read():
    """ Test read """
    assert_equal(read(data1), sequences1)
    assert_equal(read(data2), sequences2)
    assert_equal(read(data3), sequences3)

#@with_setup(setup_func, teardown_func)
def test_write():
    """ Test read """
    assert_equal(write(sequences1), data1)
    assert_equal(write(sequences2), data2)
    assert_equal(write(sequences3), data3)
