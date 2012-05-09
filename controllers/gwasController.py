###
# Controller for main presentation
#
###
from util import AppRequestHandler
import csv
from models.users import UserData
from models.study import Study
from models.annotation import Comment
from datetime import datetime

import logging
from google.appengine.api import memcache

class gwasReader(AppRequestHandler):
    def get(self):
        self.setTemplate('gwas.html')

        # check memcache for main
        rendered = memcache.get("gwas_main")
        if rendered is None:
            # make large query, to check for speed when cached
            studies = Study.all().fetch(50)
            rendered = self.render({'studies': studies})
            # logging.debug(rendered )
            # add to memchache
            if not memcache.add('gwas_main', rendered):
                logging.error("Memcache set failed.")
        
        self.response.out.write(rendered)

        # self.out({'studies': studies})

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
        comment.user = UserData.current()
        comment.date = datetime.now()
        comment.put()

        self.out({'studies': [study]})


__routes__ = [('/gwas/', gwasReader),
              ('/study/(.*)', studyPresenter)]
