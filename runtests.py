#!/usr/bin/python
import optparse
import os
import sys
import platform
# Install the Python unittest2 package before you run this script.
import unittest2

USAGE = """%prog SDK_PATH TEST_PATH
Run unit tests for App Engine apps.

SDK_PATH    Path to the SDK installation
TEST_PATH   Path to package containing test modules"""


def main(test_path, sdk_path="/usr/local/google_appengine"):
    sys.path.insert(0, sdk_path)
    import dev_appserver
    dev_appserver.fix_sys_path()
    suite = unittest2.loader.TestLoader().discover(test_path)
    unittest2.TextTestRunner(verbosity=2).run(suite)

def get_sdk_path():
    """
    get_sdk_path: Searches the system path for appengine.
    If no mention of appengine is in the PATH variable, 
    we return the OS default.
    """
    SDK_PATH = ""
    print 'Checking OS PATH string to determine SDK location'
    for path in os.environ['path'].split(';'):
        if path.find("appengine") != -1:
            SDK_PATH = path
            break
    if SDK_PATH == "":
        if platform.system() == "Windows":
            SDK_PATH = "C:\Program Files (x86)\Google\google_appengine"
        else:
            SDK_PATH = "/usr/local/google_appengine"
    print "Determined path to be " + SDK_PATH
    return SDK_PATH

def print_line():
    """
    Adds a nice dotted line to STDOUT
    """
    print '-' * 70

if __name__ == '__main__':
    parser = optparse.OptionParser(USAGE)
    options, args = parser.parse_args()
    
    if len(args) < 2:
        print_line()
        SDK_PATH = get_sdk_path()
        print_line()
        TEST_PATH = "test"
    
    if len(args) >= 1:
        TEST_PATH = args[0]
    
    if len(args) == 2:
        SDK_PATH = args[1]

    print "Running tests in following subdirectory: " + TEST_PATH
    print_line()
    main(TEST_PATH,SDK_PATH)