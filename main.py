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
# this class queries for all the articles referenced from each of the SNP-ids and 
# returns a list of dicts containing three values for each article:
# - 'title' article
# - 'abstracts' is a list of abstracts with a possible label. Ex: {label:"intro", "text":"<abstract_text>"}
# - 'PMID' (PubMed ID) of article
class pubmed(webapp.RequestHandler):
    def get(self):
        snp = self.request.get("snp")
        if snp == "":
            self.response.out.write("No snpId provided.")
            return

        # Query dbSNP for PMIDs of articles referenced by this SNP
        handle = Entrez.elink(db="pubmed", dbfrom="snp", id=snp, linkname="snp_pubmed_cited")
        dbs = Entrez.read(handle)
        # print record
        handle.close()

        # no results to return
        if len(dbs) == 0:
            self.response.out.write("No referenced articles")
            return
        
        # - should just be one, but for the hell of it, let's capture all the cases
        ref_ids = [] # holds the id of each of the referened articles
        
        # for each database with references to this snp
        for db in dbs:
            if len(db["LinkSetDb"]) == 0:
                continue               
            # linkname="snp_pubmed_cited" means just one LinkSetDb - namely PubMed
            for ref in db["LinkSetDb"][0]["Link"]:
                ref_ids.append(ref['Id'])

        # fetch all the articles with ids in ref_ids
        handle = Entrez.efetch("pubmed", id=ref_ids, retmode="xml")
        pubs = Entrez.read(handle)
        handle.close()

        # For each pubmed article, extract title, abstract and id (might not be in the same order as was queried)
        articles = []
        for pub in pubs:
            # TODO: abstracts are somewhat grouped into labels:
            # example - Background, Result and Conclusion
            base_abs = pub["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]
            categories = []
            for abstract in base_abs:
                label = None
                if hasattr(abstract, "attributes"):
                   if abstract.attributes.has_key("Label"):
                       label = abstract.attributes["Label"]

                categories.append({
                    "label" : label,
                    "text" : abstract,
                })

            articles.append({
                "title" : pub["MedlineCitation"]["Article"]["ArticleTitle"],
                "abstracts": categories,
                "pmid": pub["MedlineCitation"]["PMID"]
            })

        # print each article
        t = e.get_template('pubmeds.html')        
        self.response.out.write(t.render(pubmeds = articles))

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

        self.out({'msg':dom.toprettyxml()})


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
