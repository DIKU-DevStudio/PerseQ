from google.appengine.ext import db

# SNP models:
class snp(db.Model):
    # The essential information
    snpid = db.IntegerProperty()
    magnitude = db.FloatProperty()

    # The raw SNP pattern - optional?
    seq_data = db.Blob()

    # A little metadata
    updated = db.DateTimeProperty(auto_now_add=True)

    # Domain Tags assigned (many to many)
    domain_tags = db.ListProperty(db.Key)

    @property
    def tags(self):
        return db.get(self.domain_tags)
