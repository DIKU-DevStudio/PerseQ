from controllers.collectionController import Collection
from models.collection import SNPCollection
from models.users import UserData
from test.testUtils import TestCase

# see http://webtest.pythonpaste.org/en/latest/index.html#making-requests
class CollectionTests(TestCase):
 
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
