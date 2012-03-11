from google.appengine.ext import db

class Domain_Tag(db.Model):
    tag = db.StringProperty()
    updated = db.DateTimeProperty(auto_now_add=True)

    # Property so we can quickly query SNPs with this tag
  @property
    def members(self):
        return SNP.gql("WHERE domain_tags = :1", self.key())