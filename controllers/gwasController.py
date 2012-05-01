###
# Controller for main presentation
#
###
from Utilities import AppRequestHandler
import csv
from StringIO import StringIO
from models.study import Study

class gwasReader(AppRequestHandler):
    def get(self):
        self.setTemplate('gwas.html')
        studies = Study.all().fetch(100,0)
        self.out({'studies': studies})

class studyPresenter(AppRequestHandler):
    def get(self, i):
        self.setTemplate('study.html')
        studies = Study.gql('pubmed_id = :1', i).get()
        self.out({'studies': studies})

__routes__ = [('/gwas/', gwasReader),
              ('/study/(.*)', studyPresenter)]
