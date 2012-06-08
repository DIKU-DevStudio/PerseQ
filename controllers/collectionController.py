from Util import AppRequestHandler
import logging
from google.appengine.api import memcache
from google.appengine.ext import db
from google.appengine.api import users
from models.collection import SNPCollection
from models.users import UserData

from models.study import Snp

class diseasesGroupedBySNPs(AppRequestHandler):
    _template = "baserender.html"

    def get(self, collection_id):
        # snp_disease_collection for diseases grouped by snp
        # render = memcache.get(collection_id, namespace="snp_disease_collection")
        # if render is None:
        user = UserData.current()
        # coll = db.get(db.Key.from_path("SNPCollection", collection_id, parent=user.key()))
        coll = SNPCollection.all().ancestor(user).filter("name =", collection_id).get()
        if coll is None:
            self.redirect("/collections/")
            return

        snps = db.get(coll.snps)
        render = self.render("snp_disease_collection.html", snps=snps, collection=coll)
        # memcache.set(collection_id, render, namespace="snp_disease_collection")

        self.out(rendered=render)

class CollectionList(AppRequestHandler):
    _template = "baserender.html"
    def get(self):
        # user = users.get_current_user()
        user = UserData.current()
        collections = SNPCollection.all().ancestor(user)
        # check or update cache
        # render from template
        render = self.render("collectionlist.html",collections=collections)

        self.out(rendered=render, )

class Collection(AppRequestHandler):
    """controller handling the creation and presentation of collections"""
    _template = "baserender.html"
    def get(self):
        """view or create specified collection"""
        collection_name = self.request.get("name")
        if collection_name == "":
            self.redirect("/collections/")
            return

        # login is mandatory
        # user = users.get_current_user()
        user = UserData.current()
        # query for particular collection, creating if it does not exist
        coll = SNPCollection.all().ancestor(user).filter("name =", collection_name).get()
        if coll is None:
            # should have a 404-page
            self.redirect("/collections/")
            return

        # check or update cache
        # render from template
        render = self.render("collectionview.html",collection=coll)
        self.out(rendered=render)

    def post(self):
        # logging.info("yeah!")
        user = UserData.current()
        name = self.request.get("name")
        if name == "":
            self.redirect("/collections/")
            return

        coll = SNPCollection.all().ancestor(user).filter("name =", name).get()
        if coll is None:
            coll = SNPCollection(parent=user, name=name)
            coll.put()

        self.redirect("/collection/?name=%s" % name)
        return

class DeleteCollection(AppRequestHandler):
    def get(self, name):
        user = UserData.current()
        coll = SNPCollection.all().ancestor(user).filter("name =", name).get()
        if coll is None:
            logging.info("collection does not exist")
            self.redirect("/collections/")
            return
        logging.info("deleting collection %s" % name)
        coll.delete()
        self.redirect("/collections/")
        return
        
    def post(self, name):
        user = UserData.current()
        coll = SNPCollection.all().ancestor(user).filter("name=", name).get()
        if coll is None:
            logging.info("collection does not exist")
            self.redirect("/collections/")
            return
        logging.info("deleting collection %s" % name)

        coll.delete()
        self.redirect("/collections/")
        return

class AddSNPs(AppRequestHandler):
    _template = "baserender.html"
    def get(self, collection_name):
        # user = users.get_current_user()
        user = UserData.current()
        coll = SNPCollection.all().ancestor(user).filter("name =", collection_name).get()
        if coll is None:
            # logging.info("redirecting")
            # if collection doesn't exist - redirect to list of collections
            self.redirect("/collections/")
            return

        # we support a number of actions
        # addsnps = add a list of snps to a collection
        # action = self.requst.get("action")
        error = []
        warning = []
        # action = "addsnps"
        # if action == "addsnps":
        logging.info("adding snps to coll: %s" % collection_name)
        snps = self.request.get("snps")
        # for testing
        # snps = '10260404,3825776,12284594,rs1805007,12345'
        # comma seperated list of snps
        snpids = [snp.strip() for snp in snps.split(',')]
        # verify input:
        # - if snpid is not a valid integer-value
        #   replace with None and add to list of 'invalid' ids
        #   (might be because its prefixed with RS/rs)
        invalid_keys = []
        for i in xrange(len(snpids)):
            snp = snpids[i]
            try:
                # try to convert to int
                int(snp)
            except:
                invalid_keys.append(snp)
                snpids[i] = None

        # collect errors
        if len(invalid_keys) > 0:
            error.append("invalid keys: %s" % invalid_keys)

        # generate composite keys from the list of SNPs which are not None
        snpids = [snpid for snpid in snpids if not snpid is None]
        keys = [db.Key.from_path("Snp", snpid) for snpid in snpids]
        # logging.debug("keys to add: ")
        # retrieve all Snp-instances from these keys
        # - note: invalid keys / non-existent keys return with 'None'
        #   aka: "the SNP-id might be valid, but it is not in the GWAS catalog"
        snp_keys = db.get(keys)
        # logging.info(snp_keys)
        # logging.info("ids: %s" % [snp.key() for snp in snp_keys if snp])

        # all 'None' keys are from ids unknown to GWAS, compile a list of these

        unknown_keys = [snpids[i] for i in xrange(len(snp_keys)) if snp_keys[i] is None]
        if unknown_keys:
            warning.append("Keys unknown to the GWAS catalog: %s" % unknown_keys)

        # add keys to collection
        # list = [snp.key() for snp in snp_keys if snp]
        # logging.info("WTF: %s" % list)
        for key in [snp.key() for snp in snp_keys if snp is not None]:
            if key not in coll.snps:
                coll.snps.append(key)
        coll.count = len(coll.snps)
        coll.put()
        # try:
        #     coll.addMultipleSNP(snp_keys)
        # except Exception, e:
        #     logging.warning("issues inserting: %s" % e)
        #     error.append("Unknown issue adding snps to collection")

        # render
        logging.info("errors: %s" % error )
        logging.info("warning: %s" % warning)
        render = self.render("collectionview.html", collection=coll)
        self.out(rendered=render, error=error, warning=warning)



__routes__ = [('/collection/',Collection),
              ('/collections/',CollectionList),
              ('/addsnps/(.*)', AddSNPs),
              ('/diseases/groupebysnps/(.*)',diseasesGroupedBySNPs),
              ('/deletecollection/(.*)', DeleteCollection)]
# __routes__ = []