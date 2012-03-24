###
# Controller containing logic for requests to different SNP databases
# 
#
###
from Utilities import AppRequestHandler
from Bio import Entrez
from xml.dom.minidom import parseString

Entrez.email = 'pamad05+entrez@gmail.com'
# Given a list of comma-seperated SNP cluster ids (example: "1805007,1805008")
# this class queries for all the articles referenced from each of the SNP-ids and 
# returns a list of dicts containing three values for each article:
# - 'title' article
# - 'abstracts' is a list of abstracts with a possible label. Ex: {label:"intro", "text":"<abstract_text>"}
# - 'PMID' (PubMed ID) of article
class pubmed(AppRequestHandler):
    def get(self):
        self.setTemplate('data/pubmeds.html')
        snp = self.request.get("snp")
        if snp == "":
            self.out("No SNP id provided.")
            return

        # Query dbSNP for PMIDs of articles referenced by this SNP
        handle = Entrez.elink(db="pubmed", dbfrom="snp", id=snp, linkname="snp_pubmed_cited")
        dbs = Entrez.read(handle)
        # print record
        handle.close()

        # - should just be one, but for the hell of it, let's capture all the cases
        ref_ids = [] # holds the id of each of the referened articles

        # for each database with references to this snp
        for db in dbs:
            if len(db["LinkSetDb"]) == 0:
                continue
            # linkname="snp_pubmed_cited" means just one LinkSetDb - namely PubMed
            for ref in db["LinkSetDb"][0]["Link"]:
                ref_ids.append(ref['Id'])

        # no results to return
        if len(ref_ids) == 0:
            self.out({'msg':"No referenced articles"})
            return

        # fetch all the articles with ids in ref_ids
        handle = Entrez.efetch(db="pubmed", id=ref_ids, retmode="xml")
        pubs = Entrez.read(handle)
        handle.close()

        # For each pubmed article, extract title, abstract and id (might not be in the same order as was queried)
        articles = []
        for pub in pubs:
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
        self.out({'pubmeds': articles})

class dbSNP(AppRequestHandler):
    def get(self):
        snp = self.request.get("snp")
        if snp == "":
            self.out("No SNP_id provided")
            return

        res = Entrez.efetch("snp", id=snp, rettype="xml", retmode="text")
        dom = parseString(res.read())

        for node in dom.getElementsByTagName("Rs"):
            # print node.hasAttributes()
            rsid = node.getAttribute("rsId")
            #print rsid

        self.toXML(dom.toprettyxml(), True)


class LookUpSNP(AppRequestHandler):
    def get(self):
        self.setTemplate('data/snpedia.html')

        snp = self.request.get("snp")
        if len(snp) == 0:
            self.out({'msg':"No SNP title given (snp='')"})
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
            self.out({'msg':"No SNP title given (snp='')"})
            return

        # OMFG this is ugly..
        self.out({'msg':(result['query']["pages"][str(pageid)]["revisions"][0]["*"].encode('utf-8').replace("{{","<br>").replace("}}", "<br>").replace("\n","<br>"))})
