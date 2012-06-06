from controllers.collectionController import Collection
from models.collection import SNPCollection
from models.users import UserData
from test.testUtils import TestCase

# see http://webtest.pythonpaste.org/en/latest/index.html#making-requests
class CollectionTests(TestCase):
    # def setUp(self):
    #     # Create a WSGI application.
    #     app = webapp2.WSGIApplication(controllers.__routes__)
    #     # Wrap the app with WebTest's TestApp.
    #     self.testapp = webtest.TestApp(app)

    #     self.testbed = testbed.Testbed()
    #     # Then activate the testbed, which prepares the service stubs for use.
    #     self.testbed.activate()
    #     # Next, declare which service stubs you want to use.
    #     self.testbed.init_datastore_v3_stub()
    #     # self.testbed.init_user_stub()

    #     # FUCK ME DET TOG LANG TID AT FINDE DET HER!!!
    #     self.testbed.setup_env(
    #         USER_EMAIL = 'test@example.com',
    #         USER_ID = '123',
    #         # USER_IS_ADMIN = '1',
    #         overwrite = True)

    # def tearDown(self):
    #     self.testbed.deactivate()

    def testCreateCollection(self):
        """Test collection creation
        - try to create a new collection called 'test'
        - test that this collection has actually been added to the database afterwards
        """
        test_coll = "test"
        response = self.testapp.post('/collection/', {'name': test_coll})
        user = UserData.current()
        coll = SNPCollection.all().ancestor(user).filter("name =", test_coll).get()


        self.assertIsNotNone(coll)
        # Collection.
        # self.assertTrue(1 == 1)
