from google.appengine.ext import db

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
  	#   this is created by the snp_url class