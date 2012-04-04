#!/usr/bin/python
import optparse
import os
import sys

# sys.path.insert(0, "/usr/local/google_appengine/")
sys.path.insert(0, "/Users/pamad/dev/PerseQ/")


# Install the Python unittest2 package before you run this script.
import unittest2
# import main.DemoTestCase

USAGE = """%prog SDK_PATH TEST_PATH
Run unit tests for App Engine apps.

SDK_PATH    Path to the SDK installation
TEST_PATH   Path to package containing test modules"""

from google.appengine.ext import db
from google.appengine.ext import testbed

from models.snp import snp

# test get_or_insert
class DemoTestCase(unittest2.TestCase):
    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, "/usr/local/google_appengine/")
        dev_appserver.fix_sys_path()

    def setUp(self):
        # First, create an instance of the Testbed class.
        self.testbed = testbed.Testbed()
        # Then activate the testbed, which prepares the service stubs for use.
        self.testbed.activate()
        # Next, declare which service stubs you want to use.
        self.testbed.init_datastore_v3_stub()
        # self.testbed.init_memcache_stub()

    def tearDown(self):
        self.testbed.deactivate()

    def testGetOrInsert(self):
        id1 = snp(key_name="test1", snpid="1805007", magnitude=0.2).put()
        act_id1 = snp.get_or_insert("test1")
        self.assertEqual(act_id1.magnitude, 0.2)
        self.assertEqual(act_id1.snpid, "1805007")
        self.assertEqual(act_id1.key(), id1)

    def test2(self):
        self.assertEqual(1,1)

def main(test_path):
    # sys.path.insert(0, sdk_path)
    import dev_appserver
    dev_appserver.fix_sys_path()
    suite = unittest2.TestLoader().loadTestsFromTestCase(DemoTestCase)
    unittest2.TextTestRunner(verbosity=2).run(suite)


if __name__ == '__main__':
    # parser = optparse.OptionParser(USAGE)
    # options, args = parser.parse_args()
    # if len(args) != 2:
    #     print 'Error: Exactly 2 arguments required.'
    #     parser.print_help()
    #     sys.exit(1)
    # print os.path.dirname(__file__)
    TEST_PATH = os.path.dirname(__file__)
    main(TEST_PATH)