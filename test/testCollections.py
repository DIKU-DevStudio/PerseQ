import unittest2
# import sys
import os
import webtest # easy_install webtest
import webapp2
from google.appengine.ext import testbed
# webapp2.WSGIApplication
import controllers
from controllers.collectionController import Collection
from models.collection import SNPCollection

class SearchHandlerTest(unittest2.TestCase):
    def setUp(self):
        # Create a WSGI application.
        app = webapp2.WSGIApplication([('/', Collection)])
        # Wrap the app with WebTest's TestApp.
        self.testapp = webtest.TestApp(app)

        self.testbed = testbed.Testbed()
        # Then activate the testbed, which prepares the service stubs for use.
        self.testbed.activate()
        # Next, declare which service stubs you want to use.
        self.testbed.init_datastore_v3_stub()
        os.environ['USER_EMAIL'] = ''
        os.environ['USER_IS_ADMIN'] = ''

    # Test the handler.
    # def test(self):
    #     response = self.testapp.get('/search/')
    #     self.assertEqual(response.status_int, 200)
    #     self.assertEqual(response.normal_body, 'Hello World!')
    #     self.assertEqual(response.content_type, 'text/plain')

# see http://webtest.pythonpaste.org/en/latest/index.html#making-requests
    # def testPost(self):
    #     response = self.testapp.post('/search/SNP/',{'vars': 'values'})

    def testCreateCollection(self):
        os.environ['USER_EMAIL'] = 'test@example.com'
        # os.environ['USER_IS_ADMIN'] = '1'
        response = self.testapp.post('/?name=test')
        self.assertTrue(1 == 1)
