from google.appengine.ext import db
from models.snp import snp
# SNP URLS (one to many) -
#   this is created by the snp_url 

class snp_url(db.Model):
    snp = db.ReferenceProperty(snp, collection_name = 'urls')
    url = db.StringProperty()
    title = db.StringProperty()
    updated = db.DateTimeProperty(auto_now_add=True)

