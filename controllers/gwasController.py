"""
Controller for main presentation of GWAS data
"""
from util import AppRequestHandler
from models.users import UserData
from models.study import Study, Disease, Gene
from models.annotation import Comment
from datetime import datetime

import logging
from google.appengine.api import memcache
from google.appengine.ext import db

class diseaseList(AppRequestHandler):
    """Present a unique list of diseases, each disease linking to a page listing the studies reporting on those diseases"""
    _template = "baserender.html"
    def get(self):
        # if we need to filter later on..
        filter = self.request.get("filter") # returns name of disease to filter on
        if filter is not None:
            # TODO (pm): do a filter query - possibly caching result if queries are not terribly random
            filter = None
            pass
        # snp = self.request.get("filter") # returns name of disease to filter on
        rendered = memcache.get("diseaselist_0:50")
        if rendered is None:
            # make large query, to check for speed when cached
            diseases = Disease.all()
            # generate only the bare-bones list of diseases, ignore everything from base.html etc.
            rendered = self.render("diseaselistrender.html", diseases=diseases, filter=filter)

            # add to memchache
            if not memcache.set('diseaselist_0:50', rendered):
                logging.error("Memcache set failed for 'diseaselist_0:50'")

        # use cached data to render page with user-date etc. intact.
        self.out(rendered = rendered)

class diseaseView(AppRequestHandler):
    """Present a unique list of diseases, each disease linking to a page listing the studies reporting on those diseases"""
    _template = "baserender.html"
    def get(self, name):
        #name = self.request.get("name") # returns name of disease to filter on
        if name is None:
            # TODO (pm): return error            
            return

        # snp = self.request.get("filter") # returns name of disease to filter on
        rendered = memcache.get(name, namespace="disease")
        if rendered is None:
            # make large query, to check for speed when cached
            disease = db.get(db.Key.from_path("Disease", name))

            # generate only the bare-bones list of diseases, ignore everything from base.html etc.
            rendered = self.render("diseaseview.html", disease=disease)

            # add to memchache
            if not memcache.set(name, rendered, namespace="disease"):
                logging.error("Memcache set failed for 'disease:%s'" % name)

        # use cached data to render page with user-date etc. intact.
        self.out(rendered=rendered)

class studyList(AppRequestHandler):
    """Show a list of studies"""
    _template = 'baserender.html'
    def get(self):
        # check memcache for main
        filter = self.request.get("filter") # returns name of disease to filter on
        if filter is not None:
            # TODO (pm): do a filter query - possibly caching result if queries are not terribly random
            # filter = None
            pass
            # studies = Study.gql("WHERE pubmed_id = :1", i).get()

        rendered = memcache.get("studylist_0:50")
        if rendered is None:
            # make large query, to check for speed when cached
            studies = Study.all().fetch(50)
            rendered = self.render("studylistrender.html",studies=studies, filter=filter)
            # logging.debug(rendered )
            # add to memchache
            if not memcache.add('studylist_0:50', rendered):
                logging.error("Memcache set failed.")
        
        self.out(rendered=rendered)

        # self.out({'studies': studies})

class studyView(AppRequestHandler):
    """View a particular study"""
    def get(self, i):
        self.setTemplate('studyview.html')
        study = Study.gql("WHERE pubmed_id = :1", i).get()
        self.out(study=study)

    # Comment on a study via POST
    def post(self, i):
        self.setTemplate('studyview.html')
        study = Study.gql("WHERE pubmed_id = :1", i).get()

        comment = Comment()
        comment.study = study.key()
        comment.body = self.request.get("comment")
        comment.user = UserData.current()
        comment.date = datetime.now()
        comment.put()

        self.out(study=study)

class genePresenter(AppRequestHandler):
    """View a particular gene"""
    _template = 'gene.html'
    def get(self, gene):
        gene = Gene.gql("WHERE name = :1", gene).get()

        self.out({'gene':gene})

    # Comment on a gene via POST
    def post(self, gene):
        gene = Gene.gql("WHERE name = :1", gene).get()

        comment = Comment()
        comment.gene = gene.key()
        comment.body = self.request.get("comment")
        comment.user = UserData.current().user_id #users.get_current_user()
        comment.date = datetime.now()
        comment.put()

        self.out(gene=gene)

class commentHandler(AppRequestHandler):
    """Util comment actions"""
    def get(self, comment):
        """ Delete comment via GET """
        user = UserData.current()
        comment = Comment.get(comment)
        if user.developer or user.moderator \
                or user.user_id == comment.user.user_id:
            comment.delete()
            self.response.write("Deleted")

class commentEditor(AppRequestHandler):
    def post(self, comment):
        """ Edit comment via POST """
        comment = Comment.get(comment)
        if UserData.current().user_id == comment.user.user_id:
            comment.body = self.request.get('value')
            comment.put()
            self.response.write(comment.body)


__routes__ = [('/studies/', studyList),
              ('/study/(.*)', studyView),
              ('/comment/(.*)/delete', commentHandler),
              ('/comment/(.*)/edit', commentEditor),
              ('/diseases/', diseaseList),
              ('/disease/(.*)', diseaseView),
              ('/gene/(.*)',  genePresenter)]
