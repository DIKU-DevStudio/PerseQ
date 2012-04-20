#!/usr/bin/python
import optparse
import os
import sys
# Install the Python unittest2 package before you run this script.
import unittest2

USAGE = """%prog SDK_PATH TEST_PATH
Run unit tests for App Engine apps.

SDK_PATH    Path to the SDK installation
TEST_PATH   Path to package containing test modules"""


def main(test_path):
    sys.path.insert(0, "/usr/local/google_appengine")
    import dev_appserver
    dev_appserver.fix_sys_path()
    suite = unittest2.loader.TestLoader().discover(test_path)
    unittest2.TextTestRunner(verbosity=2).run(suite)


if __name__ == '__main__':
    main(test_path="test")