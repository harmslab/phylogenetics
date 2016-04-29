import os

from phylogenetics.tests.base import BaseTestCase

class RWTestCase(BaseTestCase):

    def setUp(self):
        """ """
        super(RWTestCase, self).setUp()
        self.fnames = []


    def tearDown(self):
        """ Remove any files created in tests. """
        for f in self.fnames:
            os.remove(f)
