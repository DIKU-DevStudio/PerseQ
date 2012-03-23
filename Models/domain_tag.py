from google.appengine.ext import db

# Tags for domains (e.g. 'Dermitology', 'Oncology')
class domain_tag(db.Model):
    tag = db.StringProperty()
    updated = db.DateTimeProperty(auto_now_add=True)

    # Property so we can quickly query SNPs with this tag
    @property
    def snps(self):
        return snp.gql("WHERE domain_tags = :1", self.key())
