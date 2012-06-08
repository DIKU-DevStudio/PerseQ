import webapp2
import unittest2
import controllers
import webtest # easy_install webtest
from google.appengine.ext import testbed

class TestCase(unittest2.TestCase):
    def setUp(self):
        # Create a WSGI application.
        app = webapp2.WSGIApplication(controllers.__routes__)
        # Wrap the app with WebTest's TestApp.
        self.testapp = webtest.TestApp(app)

        self.testbed = testbed.Testbed()
        # Then activate the testbed, which prepares the service stubs for use.
        self.testbed.activate()
        # Next, declare which service stubs you want to use.
        self.testbed.init_datastore_v3_stub()
        # self.testbed.init_user_stub()

        # FUCK ME DET TOG LANG TID AT FINDE DET HER!!!
        self.testbed.setup_env(
            USER_EMAIL = 'test@example.com',
            USER_ID = '123',
            # USER_IS_ADMIN = '1',
            overwrite = True)

    def tearDown(self):
        self.testbed.deactivate()
