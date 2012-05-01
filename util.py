from google.appengine.ext import webapp
from jinja2 import Environment, FileSystemLoader
import os
import json
# import inspect
import logging
from google.appengine.ext import db

from Bio import Entrez

"""
Date Added to Catalog
PUBMEDID
First Author
Date
Journal
Link
Study
Disease_Trait
Initial Sample Size
Replication Sample Size
Region
Chr_id
Chr_pos
Reported Gene(s)
Mapped_gene Upstream_gene_id
Downstream_gene_id
Snp_gene_ids
Upstream_gene_distance
Downstream_gene_distance
Strongest SNP-Risk Allele
SNPs
Merged
Snp_id_current
Context
Intergenic
Risk Allele Frequency
p-Value
Pvalue_mlog
p-Value (text)
OR or beta
95\% CI (text)
Platform [SNPs passing QC]
CNV"""
from models.study import Study, GWAS, Gene, Snp
import csv
def import_gwas(path="gwascatalog.txt"):
    docs = csv.DictReader(open('gwascatalog.txt','rb'), dialect='excel-tab')
    pubids = {}
    for doc in docs:
        pubid = doc["PUBMEDID"]
        if pubids.has_key(pubid):
            pubids[pubid].append(doc)
        else:
            pubids[pubid] = [doc]

    for pid, line in pubids.iteritems():
        rel = line[0]
        print pid
        study = Study.get_or_insert(pid, 
            name=rel["Study"],
            pubmed_url=rel["Link"],
            init_sample = rel["Initial Sample Size"],
            repl_sample= rel["Replication Sample Size"],
            platform = rel["Platform [SNPs passing QC]"],
            disease_trait = rel["Disease/Trait"],
            pubmed_id = pid)
            
        for rel in line:
            # init gene relation
            gene = None
            if rel["Intergenic"] == "1":
                gene = Gene.get_or_insert(rel["Snp_gene_ids"].strip(),
                    name = rel["Mapped_gene"],
                    geneid = rel["Snp_gene_ids"])
                gene.studies.append(study.key())
                gene.put()

            # init snps..
            snp = None
            if rel["Snp_id_current"] != "":
                snp = Snp.get_or_insert(rel["Snp_id_current"],
                    snpid=rel["Snp_id_current"])
                if rel["Intergenic"] == "1":
                    snp.gene = gene.key()
                snp.studies.append(study.key())
            # maybe an if here, to check or overwriting?
                snp.put()

            gwas = GWAS(study=study.key())
            if gene:
                gwas.gene = gene.key()
            gwas.p_string = rel["p-Value"]
            try:
                int(rel["p-Value"])
                gwas.p_val = int(rel["p-Value"].split("E")[1])    
            except Exception, e:
                gwas.p_val = 0
            gwas.OR_beta = rel["OR or beta"]
            gwas.riscAlleleFrequency = rel["Risk Allele Frequency"]
            gwas.intergenic = rel['Intergenic'] == '2'
            gwas.put()
    print "done"

def tear_down():
    for model in ["Snp", "Gene", "GWAS", "Study"]:
        try:
            while True:
                q = db.GqlQuery("SELECT __key__ FROM %s" % model)
                assert q.count()
                db.delete(q.fetch(200))
                # time.sleep(0.5)
        except Exception, e:
            print e
            pass
    print "done"

# given list of snpids - returns the list of related OMIM IDs
def snp_omim(snpids=None):
    if not snpids:
        return []

    handle = Entrez.elink(db="omim", dbfrom="snp", id=snpids,linkname="snp_omim")
    from_dbs = Entrez.read(handle)
    handle.close()

    omim_ids = [] # holds the id of each of the referened articles        
    # for each database with references to this snp
    for from_db in from_dbs:
        # print db
        if not from_db.has_key("LinkSetDb"):
            continue
        # print db
        for to_db in from_db["LinkSetDb"]:
            for link in to_db["Link"]:
                omim_ids.append(link['Id'])

    return omim_ids

def omim_efetch(db=None, ids=None):
    handle = Entrez.efetch(db=db, id=ids, retmode="xml")
    pubs = Entrez.read(handle)
    handle.close()
    print pubs

class jTemplate():
    _e = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(__file__), 'templates')))

    @staticmethod
    def render(template, variables, printer):
        t = jTemplate._e.get_template(template)
        printer(t.render(variables))

class AppRequestHandler(webapp.RequestHandler):
    _template = None

    def setTemplate(self, template):
        self._template = template

    def out(self, dictionary = {}):
        if(self._template == None):
            # Get template from controller / method names
            actionName = self.__class__.__name__
            self._template = actionName+".html"
        jTemplate.render(self._template, dictionary, self.response.out.write)

    def toJson(self, dictionary, prettify = False):
        data = {"json": json.dumps(dictionary)}
        if prettify:
            jTemplate.render("data/prettify/json.html", data , self.response.out.write);
        else:
            jTemplate.render("data/json.html", data , self.response.out.write);

    def toXML(self, xml, prettify = False):
        data = {'xml':xml}
        if prettify:
            jTemplate.render("data/prettify/xml.html", data, self.response.out.write);
        else:
            jTemplate.render("data/xml.html", data,self.response.out.write );
