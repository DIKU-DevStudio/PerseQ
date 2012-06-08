import sys
from datetime import date
from google.appengine.ext import db
from google.appengine.ext import testbed

from models.study import Snp, Gene, Study, GWAS, Disease
from test.testUtils import TestCase
from Util import AddStudyDocument, populate



# test get_or_insert
class SearchHandlerTest(TestCase):

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
        result = Study.search_todict('123')
        self.assertEqual(obj['name'], result[0]['name']) #result[0], obj)

        self.tSearchPage()
        self.tSearchStudy()
        self.tSearchDisease()

    def tSearchPage(self):
        search = '/search/'
        response = self.testapp.get('/search/')
        self.assertEqual(response.status_int, 200)

    def tSearchStudy(self):
        search = '/search/study/'
        response = self.testapp.get(search, status=404)
        self.assertEqual(response.status_int, 404)

        study = '17463249' # should exist
        response = self.testapp.get(search+study)
        self.assertEqual(response.status_int, 200)
        self.assertRegexpMatches(response.body, study)

    def tSearchDisease(self):
        search = '/search/disease/'
        response = self.testapp.get(search, status=404)
        self.assertEqual(response.status_int, 404)

        disease = 'cancer'
        response = self.testapp.get(search+disease)
        self.assertEqual(response.status_int, 200)
        self.assertRegexpMatches(response.body, disease)
