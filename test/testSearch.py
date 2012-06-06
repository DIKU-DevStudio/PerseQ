import sys
from datetime import date
from google.appengine.ext import db
from google.appengine.ext import testbed

from models.study import Snp, Gene, Study, GWAS, Disease
from test.testUtils import TestCase
from Util import AddStudyDocument
from Util import purge,populate



# test get_or_insert
class SearchHandlerTest(TestCase):
    """def setUpClass():
        super(SearchHandlerTest, self).setUp()
        purge()
        populate()

    def tearDownClass():
        purge()"""

    def setUp(self):
        super(SearchHandlerTest, self).setUp()
        populate()

    def testSearchDocuments(self):
        """Test search document logic
        - try to add a study and fetch it using full-text search
        """
        disease_name = 'Pizza'
        disease = Disease.get_or_insert(disease_name,
            name=disease_name)
        obj = {'pubmed_id':'123',
                'name':'How is the pizza?',
                'disease_ref':disease,
                'disease_trait':'Pizza'}

        study = Study.get_or_insert(obj['pubmed_id'],
                    name=obj['name'],
                    disease_ref = obj['disease_ref'],
                    disease_trait = obj['disease_trait'],
                    pubmed_id = obj['pubmed_id'],
                    date = date.today())
        study.put()

        AddStudyDocument(study)
        result = Study.search_todict('pizza')
        self.assertEqual(result, obj)

    def testSearchPage(self):
        search = '/search/'
        response = self.testapp.get(search)
        self.assertEqual(response.status_int, 200)

    def testSearchStudy(self):
        search = '/search/study/'
        response = self.testapp.get(search, status=404)
        self.assertEqual(response.status_int, 404)

        snp = '4343' # should exist
        response = self.testapp.get(search+snp)
        self.assertEqual(response.status_int, 200)
        self.assertRegexpMatches(response.body, snp)

    def testSearchDisease(self):
        search = '/search/disease/'
        response = self.testapp.get(search, status=404)
        self.assertEqual(response.status_int, 404)

        disease = 'cancer'
        response = self.testapp.get(search+disease)
        self.assertEqual(response.status_int, 200)
        self.assertRegexpMatches(response.body, disease)
