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
Mapped_gene
Upstream_gene_id
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
def populate(path="gwascatalog.txt"):
    docs = csv.DictReader(open('gwascatalog.txt','rb'), dialect='excel-tab')
    pubids = {}
    for doc in docs:
        pubid = doc["PUBMEDID"]
        if pubids.has_key(pubid):
            pubids[pubid].append(doc)
        else:
            pubids[pubid] = [doc]

    limit = 300
    i = 0
    for pid, line in pubids.iteritems():
        i += 1
        if i == limit:
            break
        rel = line[0]
        # print pid
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
            up_gene = None
            down_gene = None

            # if intergenic gwas, we go nuts..
            intergenic = rel["Intergenic"] == "2"

            if not intergenic:
                geneid = rel["Snp_gene_ids"].strip()
                if geneid != "":
                # not intergenic => one direct gene
                    names = rel["Mapped_gene"].split(" - ")
                    if len(names) == 1:
                        gene = Gene.get_or_insert(rel["Snp_gene_ids"].strip(),
                            name = rel["Mapped_gene"],
                            geneid = rel["Snp_gene_ids"])

                        gene.studies.append(study.key())
                        gene.put()
            else:
                # up and downstream genes must be set
                down_id = rel["Downstream_gene_id"].strip()
                up_id = rel["Upstream_gene_id"].strip() 

                if up_id != "" and down_id != "":
                    up_down_names = rel["Mapped_gene"].split(" - ")
                    if len(up_down_names) < 2:
                        # gene = NR / NS or whatever..
                        up_down_names = ["N/A", "N/A"]
                    
                    # create upstream gene
                    up_name = up_down_names[0]
                    down_name = up_down_names[1]
                    print down_name, up_name
                    # assert("-" not in up_name)
                    # assert("-" not in down_name)

                    assert(up_id != "")
                    assert(down_id != "")

                    up_gene = Gene.get_or_insert(up_id,
                        name = up_name, 
                        geneid = up_id)
                    # up_gene.study.append(study.key())
                    down_gene = Gene.get_or_insert(down_id,
                        name = down_name, 
                        geneid = down_id)
                    # up_gene.study.append(study.key())

            # init snps..
            snp = None
            if rel["Snp_id_current"] != "":
                # non=blank snp
                snpid = rel["Snp_id_current"].strip()

                snp = Snp.get_or_insert(snpid,
                    snpid=snpid)
                if gene:
                    snp.gene = gene.key()
                snp.studies.append(study.key())
            # maybe an if here, to check for overwriting?
                snp.put()
            
            # init gwas relation
            gwas = GWAS(study=study.key(), intergenic=intergenic)
            if intergenic:
                if down_gene and up_gene:
                    gwas.downstream = down_gene.key()
                    gwas.upstream = up_gene.key()
            else:
                if gene:
                    gwas.gene = gene.key()

            # parse remaining gwas information
            gwas.p_string = rel["p-Value"]
            gwas.snps = rel["Snp_id_current"].strip()
            # parse out the exponent: 6E-8 => -8
            try:
                # test that p-Value is actually a float before parsing out
                float(rel["p-Value"])
                gwas.p_val = int(rel["p-Value"].split("E")[1])    
            except Exception, e:
                print e
                # forces the filter to downgrade this gwas wrt. p-value
                gwas.p_val = 0
            # could be interesting
            gwas.strongest_snp_risk_allele = \
                rel["Strongest SNP-Risk Allele"].strip()
            gwas.OR_beta = rel["OR or beta"].strip()
            gwas.riscAlleleFrequency = \
                rel["Risk Allele Frequency"].strip()

            gwas.put()
    print "done"

def purge():
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
