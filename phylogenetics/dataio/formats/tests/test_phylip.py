from nose.tools import assert_equal

from ..phylip import read, write

data1 = """\
2    45
XX00000
ASGASHAHAHAE---------ASDHASDF---ASD-ADFASHASS
XX00000001
SDHASASHADSGASDFASDF----------ASDGHASSFASGASD\
"""

phylip_data1 = [
    ("XX00000", "ASGASHAHAHAE---------ASDHASDF---ASD-ADFASHASS"),
    ("XX00000001", "SDHASASHADSGASDFASDF----------ASDGHASSFASGASD")
]

data2 = """\
2    45
XX00000   ASGASHAHAHAE---------ASDHASDF---ASD-ADFASHASS
XX00000001SDHASASHADSGASDFASDF----------ASDGHASSFASGASD\
"""

phylip_data2 = [
    ("XX00000", "ASGASHAHAHAE---------ASDHASDF---ASD-ADFASHASS"),
    ("XX00000001", "SDHASASHADSGASDFASDF----------ASDGHASSFASGASD")
]


def test_read():
    """ Test read """
    assert_equal(read(data1), phylip_data1)
    assert_equal(read(data2), phylip_data2)

def test_write():
    assert_equal(write(phylip_data1), data1)
