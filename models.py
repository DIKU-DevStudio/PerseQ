from google.appengine.ext import db

# SNP models:
class SNP(db.Model):
	# The essential information
	id = db.IntegerProperty()
    magnitude = db.FloatProperty()

    # The raw SNP pattern - optional?
    seq_data = db.Blob()

    # A little metadata
    updated = db.DateTimeProperty(auto_now_add=True)

    # Domain Tags assigned (many to many)
	domain_tags = db.ListProperty(db.Key)
  
  	# SNP URLS (one to many) -
  	#   this is created by the snp_url 

 class SNP_URL(db.Model):
    snp = db.ReferenceProperty(SNP, collection_name = 'urls'
    url = db.StringProperty()
    title = db.StringProperty()
    updated = db.DateTimeProperty(auto_now_add=True)

# Tags for domains (e.g. 'Dermitology', 'Oncology')
class Domain_Tag(db.Model):
    tag = db.StringProperty()
    updated = db.DateTimeProperty(auto_now_add=True)

    # Property so we can quickly query SNPs with this tag
  @property
    def snps(self):
        return SNP.gql("WHERE domain_tags = :1", self.key())