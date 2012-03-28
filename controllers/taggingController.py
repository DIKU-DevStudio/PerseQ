
from Utilities import AppRequestHandler
from models.snp import snp
from models.snp_url import snp_url
from models.domain_tag import domain_tag

class addTag(AppRequestHandler):
    def get(self):
        # Get tags logic
        snpStr = self.request.get("snp").lstrip('rsRS')
        # Get if exists, make dummy if not
        snpObj = snp.get_by_key_name(snpStr)
        if snpObj == None:
            snpObj = snp()
            snpObj.snpid = int(snpStr)
        tags = domain_tag.all()
        self.out({'tags':tags, 'snp':snpObj})

    def post(self):
        # Add tag logic
        snpStr = self.request.get("snp")
        tagStr = self.request.get("tag")

        snpObj = snp.get_or_insert(snpStr, snpid=int(snpStr))
        if snpObj == None:
            snpObj = snp(snpid = int(snpStr))

        tagObj = domain_tag.get_or_insert(tagStr, tag = tagStr)
        if tagObj == None:
            tagObj = tag(tag=tagStr)
        tagObj.put()

        methodStr = self.request.get("method")
        if methodStr == 'delete':
            snpObj.domain_tags.remove(tagObj.key())
        else:
            snpObj.domain_tags.append(tagObj.key())
        snpObj.put()

        self.out();

__routes__ = [('/tag/', addTag)]
