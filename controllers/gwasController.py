###
# Controller for main presentation
#
###
from util import AppRequestHandler
import csv
from StringIO import StringIO
from google.appengine.api import users
from models.study import Study
from models.annotation import Comment
from datetime import datetime

class gwasReader(AppRequestHandler):
    def get(self):
        self.setTemplate('gwas.html')
        studies = Study.all().fetch(10,0)
        self.out({'studies': studies})

class studyPresenter(AppRequestHandler):
    def get(self, i):
        self.setTemplate('study.html')
        study = Study.gql("WHERE pubmed_id = :1", i).get()
        self.out({'studies': [study]})

    def post(self, i):
        self.setTemplate('study.html')
        study = Study.gql("WHERE pubmed_id = :1", i).get()

        comment = Comment()
        comment.study = study.key()
        comment.body = self.request.get("comment")
        comment.user = users.get_current_user()
        comment.date = datetime.now()
        comment.put()

        self.out({'studies': [study]})


__routes__ = [('/gwas/', gwasReader),
              ('/study/(.*)', studyPresenter)]
