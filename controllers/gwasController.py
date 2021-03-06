###
# Controller for main presentation
#
###
from Util import AppRequestHandler
from models.users import UserData
from models.study import Study, Disease, Gene
from models.annotation import Comment
from datetime import datetime

import logging
from google.appengine.api import memcache
from google.appengine.ext import db

# from he3.src.he3.db.tower.pagin import PagedQuery
from he3.db.tower.paging import PagedQuery

class diseaseList(AppRequestHandler):
    """Present a unique list of diseases, each disease linking to a page listing the studies reporting on those diseases"""
    def __init__(self, *args, **kwargs):
        super(diseaseList, self).__init__(*args, **kwargs)
        self._template = "baserender.html"
        self.diseasequery = PagedQuery(Disease.all(), 30)
        self.count = self.diseasequery.page_count()

    def get(self):
        # if we need to myfilter later on..
        myfilter = self.request.get("filter", "") # returns name of disease to myfilter on
        if myfilter != "":
            diseases = Disease.search_todict('"'+myfilter+'"')
            return self.out(rendered = self.render("diseaselistrender.html", diseases = diseases, filter=myfilter))

        # snp = self.request.get("filter") # returns name of disease to myfilter on

        # get number of pages
        
        # get page number
        # if not specified => pagenr = 1
        page = self.request.get("page", "")
        pagenr = None
        if page == "":
            pagenr = 1
        else:
            try:
                pagenr = int(page)
            except:
                pagenr = 1

        # get pagenumber from diseasequery
        diseases = self.diseasequery.fetch_page(pagenr)

        # generate only the bare-bones list of diseases, ignore everything from base.html etc.
        rendered = self.render("diseaselistrender.html", diseases=diseases, page=pagenr, count=self.count)

        # # add to memchache
        # if not memcache.set('diseaselist_0:50', rendered):
        #     logging.error("Memcache set failed for 'diseaselist_0:50'")

        # use cached data to render page with user-date etc. intact.
        self.out(rendered = rendered)

class diseaseView(AppRequestHandler):
    """Present a unique list of diseases, each disease linking to a page listing the studies reporting on those diseases"""
    _template = "baserender.html"
    def get(self, name):
        #name = self.request.get("name") # returns name of disease to myfilter on
        if name is None:
            self.error(404)
            return

        # snp = self.request.get("filter") # returns name of disease to myfilter on
        rendered = memcache.get(name, namespace="disease")
        if rendered is None:
            # make large query, to check for speed when cached
            disease = db.get(db.Key.from_path("Disease", name))
            studies = db.get(disease.studies)
            if studies is None:
                logging.info("ERRRRRRROR: no studies!")
            else:
                logging.info("%s" % studies)

            # generate only the bare-bones list of diseases, ignore everything from base.html etc.
            rendered = self.render("diseaseview.html", disease=disease, studies=studies)

            # add to memchache
            if not memcache.set(name, rendered, namespace="disease"):
                logging.error("Memcache set failed for 'disease:%s'" % name)

        # use cached data to render page with user-date etc. intact.
        self.out(rendered=rendered)



class studyList(AppRequestHandler):
    """Show a list of studies"""
    def __init__(self, *args, **kwargs):
        super(studyList, self).__init__(*args, **kwargs)
        self._template = 'baserender.html'
        self.studyquery = PagedQuery(Study.all(), 10)
        self.count = self.studyquery.page_count()

    def get(self):
        # check memcache for main
        myfilter = self.request.get("filter", "") # returns name of disease to filter on
        if myfilter != "" and myfilter is not None:
            studies = Study.search_todict('"'+myfilter+'"')
            return self.out(rendered = self.render("studylistrender.html", studies = studies, filter=myfilter))

        page = self.request.get("page", "")
        pagenr = None
        if page == "":
            pagenr = 1
        else:
            try:
                pagenr = int(page)
            except:
                pagenr = 1

        # get pagenumber from diseasequery
        studies = self.studyquery.fetch_page(pagenr)

        # rendered = memcache.get("studylist_0:50")
        # if rendered is None:
            # make large query, to check for speed when cached
            # studies = Study.all().fetch(100)
        rendered = self.render("studylistrender.html",studies=studies, myfilter=myfilter, page=pagenr, count=self.count)
            # logging.debug(rendered )
            # add to memchache
            # if not memcache.add('studylist_0:50', rendered):
            #     logging.error("Memcache set failed.")

        self.out(rendered=rendered)

class studyView(AppRequestHandler):
    """View a particular study"""
    def get(self, pubmed_id):
        self.setTemplate('studyview.html')
        study = Study.get_by_key_name(pubmed_id)
        if study is None:
            self.redirect('/studies/')
            return
        self.out(study=study)

    # Comment on a study via POST
    def post(self, pubmed_id):
        self.setTemplate('studyview.html')
        study = Study.get_by_key_name(pubmed_id)
        if study is None:
            self.redirect('/studies/')
            return

        comment = Comment()
        comment.study = study.key()
        comment.body = self.request.get("comment")
        comment.user = UserData.current()
        comment.date = datetime.now()
        comment.put()

        self.out(study=study)

class SNPView(AppRequestHandler):
    _template = "baserender.html"
    def get(self, snpid):
        if snpid == "":
            self.error(404)
            return

        render = memcache.get(snpid, namespace="snp")
        if render is None:
            key = db.Key.from_path("Snp", snpid)
            snp = db.get(key)
            if snp is None:
                self.error(404)
                return
            render = self.render("snpview.html", snp=snp)
            memcache.set(snpid, render, namespace="snp")

        self.out(rendered=render)

class genePresenter(AppRequestHandler):
    """View a particular gene"""
    _template = 'gene.html'
    def get(self, gene):
        gene = Gene.gql("WHERE name = :1", gene).get()
        self.out(gene=gene)

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
              ('/gene/(.*)',  genePresenter),
              ('/snp/(.*)', SNPView)]
