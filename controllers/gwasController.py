###
# Controller for main presentation
#
###
from util import AppRequestHandler
import csv
from models.users import UserData
from models.study import Study, Disease
from models.annotation import Comment
from datetime import datetime

import logging
from google.appengine.api import memcache

class diseaseList(AppRequestHandler):
    _template = "baserender.html"
    """Present a unique list of diseases, each disease linking to a page listing the studies reporting on those diseases"""
    def get(self):
        filter = self.request.get("filter") # returns name of disease to filter on
        if filter is not None:
            # TODO (pm): do a filter query - possibly caching result if queries are not terribly random
            filter = None
            pass
        # snp = self.request.get("filter") # returns name of disease to filter on
        rendered = memcache.get("disease_0:50")
        if rendered is None:
            # make large query, to check for speed when cached
            diseases = Disease.all()
            # generate only the bare-bones list of diseases, ignore everything from base.html etc.
            rendered = self.render({'diseases': diseases, 'filter': filter}, "diseaselistrender.html")

            # add to memchache
            if not memcache.set('disease_0:50', rendered):
                logging.error("Memcache set failed for 'disease_0:50'")

        # use cached data to render page with user-date etc. intact.
        self.out({"rendered":rendered})

class studyList(AppRequestHandler):
    _template = 'baserender.html'
    def get(self):
        # check memcache for main
        filter = self.request.get("filter") # returns name of disease to filter on
        if filter is not None:
            # TODO (pm): do a filter query - possibly caching result if queries are not terribly random
            filter = None
            pass

        rendered = memcache.get("studylist_0:50")
        if rendered is None:
            # make large query, to check for speed when cached
            studies = Study.all().fetch(50)
            rendered = self.render({'studies': studies, 'filter': filter}, "studylistrender.html")
            # logging.debug(rendered )
            # add to memchache
            if not memcache.add('gwas_main', rendered):
                logging.error("Memcache set failed.")
        
        self.out({'rendered':rendered})

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


__routes__ = [('/studies/', studyList),
              ('/study/(.*)', studyPresenter),
              ('/diseases/', diseaseList)]
