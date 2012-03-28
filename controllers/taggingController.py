
from Utilities import AppRequestHandler
from models.snp import snp
from models.snp_url import snp_url
from models.domain_tag import domain_tag

class addTag(AppRequestHandler):
    def get(self):
        # Get tags logic
        self.setTemplate('Main/autocomplete.html')
        snpStr = self.request.get("snp")
        snpObj = snp.gql("WHERE snpid = '"+snpStr+"'").get()
        tags = domain_tag.all()
        #tags = snpObj.tags
        self.out({'tags':tags, 'snp':snpObj})
    def post(self):
        # Add tag logic
        snpStr = self.request.get("snp")
        tagStr = self.request.get("tag")
        snpObj = snp.gql('WHERE snpid = '+snpStr).get()
        tagObj = domain_tag.gql("WHERE tag = '"+tagStr+"'").get()
        methodStr = self.request.get("method")

        if methodStr == 'delete':
            snpObj.domain_tags.remove(tagObj.key())
        else:
            snpObj.domain_tags.append(tagObj.key())
        snpObj.put()

        self.out();

__routes__ = [('/tag/', addTag)]
