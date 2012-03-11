from google.appengine.ext import db

class SNP_URL(db.Model):
    snp = db.ReferenceProperty(SNP, collection_name = 'snp_urls'
    url = db.StringProperty()
    title = db.StringProperty()
    updated = db.DateTimeProperty(auto_now_add=True)
