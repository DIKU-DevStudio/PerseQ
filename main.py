#!/usr/bin/env python
#
# Copyright 2007 Google Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
from google.appengine.ext import webapp
from google.appengine.ext.webapp import util
from wikitools import wiki, page, api
import jinja2
from jinja2 import Environment, FileSystemLoader
import os

from xml.dom.minidom import parseString
from Bio import Entrez
from models import snp, snp_url, domain_tag

Entrez.email = 'pamad05+entrez@gmail.com'
# import webapp2
# import pysam

import pprint
pp = pprint.PrettyPrinter(indent=4)

e = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(__file__), 'templates')))

class snpsList(webapp.RequestHandler):
    def get(self):
        t = e.get_template('snplist.html')
        snpids = snp.all()
        snpids.order("snpid")
        
        self.response.out.write(t.render(snps = snpids))

# Given a list of comma-seperated SNP cluster ids (example: "1805007,1805008")
# this method returns a list of article-ids of all articles referencing the given snps
# result for snp="1805007":
# [{u'DbFrom': 'snp', u'IdList': ['1805007'], u'LinkSetDbHistory': [], u'LinkSetDb': [{u'DbTo': 'pubmed', u'Link': [{u'Id': '21700618'}, {u'Id': '20670983'}, {u'Id': '20585627'}, {u'Id': '20042077'}, {u'Id': '19884608'}, {u'Id': '17999355'}, {u'Id': '17952075'}], u'LinkName': 'snp_pubmed_cited'}]}]
class pubmed(webapp.RequestHandler):
    def get(self):
        snp = self.request.get("snp")
        if not snp:
            snp = "1805007"

        handle = Entrez.elink(db="pubmed", dbfrom="snp", id=snp, linkname="snp_pubmed_cited")
        dbs = Entrez.read(handle)
        # print record
        handle.close()

        # we only get results from the snp-database
        if len(dbs) == 0:
            self.response.out.write("No referenced articles")
            return

        ref_ids = []
        for db in dbs:
            if len(db["LinkSetDb"]) == 0:
                continue               
            # linkname="snp_pubmed_cited" means the uery only returns results to pubmed
            for ref in db["LinkSetDb"][0]["Link"]:
                ref_ids.append(ref['Id'])

        if len(ref_ids) > 0:
            self.response.out.write(",".join(ref_ids))
        else:
            self.response.out.write("No referenced articles")   

class dbSNP(webapp.RequestHandler):
    def get(self):
        snp = self.request.get("snp")
        if not snp:
            snp = "1805007"
        
        res = Entrez.efetch("snp", id=snp, rettype="xml", retmode="text")
        dom = parseString(res.read())

        for node in dom.getElementsByTagName("Rs"):
            # print node.hasAttributes()
            rsid = node.getAttribute("rsId")
            print rsid

        self.response.out.write(dom.toprettyxml())
        # self.response.out.write(res.read().replace("\n","<br>"))


class LookUpSNP(webapp.RequestHandler):
    def get(self):
        snp = self.request.get("snp")
        if len(snp) == 0:
            self.response.out.write("No SNP title given (snp='')")
            return
        site = wiki.Wiki("http://bots.snpedia.com/api.php")
        params = {
            'action': 'query',
            'prop': 'revisions',
            'rvprop': 'content',
            'rvlimit': '1',
            'titles':snp
        }
        req = api.APIRequest(site, params)
        result = req.query(querycontinue=False)

        # if pageid == -1 <=> title=snp does not exist
        pageid = int(result['query']['pages'].keys()[0])
        if pageid == -1:
            self.response.out.write("SNP title does not exist (NoPage)")
            return

        # OMFG this is ugly..
        self.response.out.write(result['query']["pages"][str(pageid)]["revisions"][0]["*"].encode('utf-8').replace("{{","<br>").replace("}}", "<br>").replace("\n","<br>"))

class TestJinja(webapp.RequestHandler):
    def get(self):
        t = e.get_template('hello.html')
        self.response.out.write(t.render(msg = 'Hello World!!!'))


def main():
    application = webapp.WSGIApplication([('/', snpsList),
                                          ('/test/', TestJinja),
                                          ('/dbsnp/', dbSNP),
                                          ('/pubmed/', pubmed)
                                          ], debug=True)
    util.run_wsgi_app(application)

if __name__ == '__main__':
    main()
