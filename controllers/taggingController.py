"""Tagging controller is deprecated.
Allows basic tagging of SNPs only"""

from Util import AppRequestHandler
from models.snp import snp
from models.snp_url import snp_url
from models.domain_tag import domain_tag

class addTag(AppRequestHandler):
    """Tag request handler"""
    def get(self):
        """List tags"""
        # Get tags logic
        snpStr = self.request.get("snp").lstrip('rsRS')
        # Get if exists, make dummy if not
        snpObj = snp.get_by_key_name(snpStr)
        if snpObj == None:
            snpObj = snp(snpid = snpStr)
        tags = domain_tag.all()
        self.out({'tags':tags, 'snp':snpObj})

    def post(self):
        """Add tag via POST params; snp and tag"""
        # Add tag logic
        snpStr = self.request.get("snp")
        tagStr = self.request.get("tag")

        snpObj = snp.get_by_key_name(snpStr)
        if snpObj == None:
            snpObj = snp(snpid = snpStr, key_name=snpStr)

        tagObj = domain_tag.get_by_key_name(tagStr)
        if tagObj == None:
            tagObj = domain_tag(tag=tagStr, key_name=tagStr)
        tagObj.put()

        methodStr = self.request.get("method")
        if methodStr == 'delete':
            snpObj.domain_tags.remove(tagObj.key())
        else:
            snpObj.domain_tags.append(tagObj.key())

        snpObj.put()
        self.out();

__routes__ = [('/tag/', addTag)]
