"""Utility functionality used across multiple controllers"""

# from google.appengine.ext import webapp
import webapp2
from jinja2 import Environment, FileSystemLoader
import os
import simplejson
import logging
from google.appengine.ext import db
from google.appengine.api import memcache, users, search
from datetime import datetime
from models.users import UserData
from models.study import Study, GWAS, Gene, Snp, Disease

from Bio import Entrez
import StringIO


import csv
from google.appengine.api import memcache
def reset():
    """resets the mem-cache:
    (auto-reset on deploy)[http://stackoverflow.com/questions/1983556/how-can-i-have-google-app-engine-clear-memcache-every-time-a-site-is-deployed]
    """
    memcache.flush_all()

def AddDiseaseDocument(d):
    doc = search.Document(doc_id="".join(d.name.split()), # doesnt allow spaces
        fields=[
            search.TextField(name='name', value=d.name),
            ])
    search.Index(name=d._index).add(doc)


def AddStudyDocument(study):
    doc = search.Document(doc_id=study.pubmed_id, # Treat pubmed_id as key
        fields=[
            search.TextField(name='name', value=study.name),
            search.TextField(name='disease_trait', value=study.disease_trait),
            search.TextField(name='id', value=study.pubmed_id)
            ])
    search.Index(name=study._index).add(doc)

def AddSNPDocument(snp):
    doc = search.Document(
        fields=[
            search.TextField(name='snpid', value=snp.snpid),
            ])
    search.Index(name=snp._index).add(doc)

def AddGeneDocument(gene):
    doc = search.Document(doc_id=gene.geneid,
        fields=[
            search.TextField(name='id', value=gene.geneid),
            search.TextField(name='name', value=gene.name),
            search.TextField(name='alias', value=str(gene.alias)),
            ])
    search.Index(name=gene._index).add(doc)

def populate(path="gwascatalog.txt", limit=100):
    """Populate the database with data from gwascatalog.txt - one hell of an import!
    - create models of SNP, Gene, Study and Disease and relations between the objects
    - based on the following field-names in the TSV (excel):
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
    docs = csv.DictReader(open(path,'rb'), dialect='excel-tab')
    pubids = {}
    # read all GWAS into a dictionary, using pubmed_id as key
    # - collecting all GWAS assorted with the same study under same key
    for doc in docs:
        pubid = doc["PUBMEDID"]
        if pubids.has_key(pubid):
            pubids[pubid].append(doc)
        else:
            pubids[pubid] = [doc]

    i = 0
    for pid, line in pubids.iteritems():
        i += 1
        if i == limit:
            break
        rel = line[0]
        
        # create or get disease with disease_name as key
        disease_name = rel["Disease/Trait"].strip().lower()
        disease = Disease.get_or_insert(disease_name,
            name=disease_name)
        AddDiseaseDocument(disease)

        # create or get study with study_id
        study = Study.get_or_insert(pid,
            name=rel["Study"].strip(),
            disease_trait = disease_name,
            disease_ref = disease,
            pubmed_id = pid)
        
        # parse date of study
        date = datetime.strptime(rel["Date"].strip(), "%m/%d/%Y").date()
        # populate model
        study.date = date
        study.pubmed_url=rel["Link"].strip()
        study.init_sample = rel["Initial Sample Size"].strip()
        study.repl_sample= rel["Replication Sample Size"].strip()
        study.platform = rel["Platform [SNPs passing QC]"].strip()
        study.put()
        AddStudyDocument(study)

        for rel in line:
            # A gwas has either a direct gene or a 
            # down-stream and up-stream gene
            gene = None
            up_gene = None
            down_gene = None

            # The retards at GWAS use 1 == intergenic, 2 == not intergenic
            # ... no seriously, _retards_
            intergenic = rel["Intergenic"].strip() == "2"

            if not intergenic:
                geneid = rel["Snp_gene_ids"].strip()
                if geneid != "":
                # not intergenic => one direct gene
                    names = rel["Mapped_gene"].split(" - ")
                    if len(names) == 1:
                        gene = Gene.get_or_insert(geneid,
                            name = rel["Mapped_gene"],
                            geneid = geneid)
                        if not study.key() in gene.studies:
                            gene.studies.append(study.key())
                        if not disease.key() in gene.diseases:
                            gene.diseases.append(disease.key())
                        gene.put()
                        AddGeneDocument(gene)
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
                    # print down_name, up_name
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
            snpid = rel["Snp_id_current"].strip()

            if snpid != "":
                # non=blank snp
                try:
                    # create relation only if a 'single' snp is in the gwas
                    int(snpid)
                    snp = Snp.get_or_insert(snpid,
                        snpid=snpid)
                    if gene:
                        snp.gene = gene
                    if not study.key() in snp.studies:
                        snp.studies.append(study.key())
                    if not disease.key() in snp.diseases:
                        snp.diseases.append(disease.key())
                    snp.put()
                    AddSnpDocument(snp)
                except:
                    # haplotype?
                    snpid = "N/A"
            # if no gene or snp relation is mentioned - ignore and just insert study
            if (gene is None or up_gene is None) and snp is None:
                print "skipping gwas"
                continue
            # init gwas
            gwas = GWAS(study=study,
                intergenic=intergenic,
                disease=disease,
                snp=snp)

            # if SNP is intergenic, save up/down-stream, else gene
            if intergenic:
                if down_gene is not None and up_gene is not None:
                    gwas.downstream = down_gene
                    gwas.upstream = up_gene
            else:
                # could be None
                gwas.gene = gene

            # parse remaining gwas information
            gwas.p_string = rel["p-Value"].strip()
            # could be none
            gwas.snps = snpid

            # parse out the exponent: 6E-8 => -8
            try:
                # test that p-Value is actually a float before parsing out
                float(rel["p-Value"])
                gwas.p_val = int(rel["p-Value"].split("E-")[1])  
            except Exception, e:
                # forces the filter to downgrade this gwas wrt. p-value
                gwas.p_val = 0
            # could be interesting
            gwas.strongest_snp_risk_allele = \
                rel["Strongest SNP-Risk Allele"].strip()
            # gwas.CI = rel["95% CI (text)"].strip()
            gwas.OR_beta = rel["OR or beta"].strip()
            gwas.riscAlleleFrequency = \
                rel["Risk Allele Frequency"].strip()
            gwas.put()

    memcache.flush_all()
    print "done"

def purge():
    """Clear the database to make ready for (re-)population"""
    for model in ["Snp", "Gene", "GWAS", "Study", "Disease"]:
        try:
            while True:
                q = db.GqlQuery("SELECT __key__ FROM %s" % model)
                assert q.count()
                db.delete(q.fetch(200))
        except Exception, e:
            print e
            pass
    print "Datastore cleared"

    for model in [Snp, Gene,GWAS, Study, Disease]:
        index = search.Index(name=model._index)
        while True:
            # Get a list of documents populating only the doc_id field and extract the ids.
            document_ids = [document.doc_id
                            for document in index.list_documents(ids_only=True)]
            if not document_ids:
                break
            # Remove the documents for the given ids from the Index.
            index.remove(document_ids)
    print "Fulltext docs cleared"

def snp_omim(snpids=None):
    """given list of snpids - returns the list of related OMIM IDs"""
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
    """Omin fetch"""
    handle = Entrez.efetch(db=db, id=ids, retmode="xml")
    pubs = Entrez.read(handle)
    handle.close()
    print pubs

class jTemplate():
    """Template helper, sets the template base path and uses a given renderer on
    the template."""

    _e = Environment(loader=FileSystemLoader(os.path.join(
        os.path.dirname(__file__), 'templates')))

    @staticmethod
    def render(template, variables, printer):
        t = jTemplate._e.get_template(template)
        printer(t.render(variables))

env = Environment(loader=FileSystemLoader(os.path.join(
            os.path.dirname(__file__), 'templates')))
class AppRequestHandler(webapp2.RequestHandler):
    """Base class for controllers"""
    _template = None

    def render(self, template = None, **args):
        if template is None:
            raise Exception("No 'template' argument to render from")
        temp = env.get_template(template)
        return temp.render(**args)

    def setTemplate(self, template):
        """Set the template to be used for calls rendering"""
        self._template = template

    def out(self, **dictionary):
        """Render the template, passing a dict as inputs"""
        if self._template is None:
            # Get template from controller / method names
            actionName = self.__class__.__name__
            self._template = actionName+".html"

        dictionary['user'] =  UserData.current()
        dictionary['user_logout'] = users.create_logout_url('/')
        dictionary['user_login'] = users.create_login_url('/')
        dictionary['user_admin'] = users.is_current_user_admin()

        temp = env.get_template(self._template)
        self.response.write(temp.render(dictionary))

    def toJson(self, dictionary, prettify = False):
        """Display JSON data template.
        Prettify flag tells whether to use the google code prettify markup"""
        enc = simplejson.JSONEncoder()
        data = {"json": enc.encode(dictionary)}
        if prettify:
            jTemplate.render("data/prettify/json.html", data , self.response.out.write);
        else:
            jTemplate.render("data/json.html", data , self.response.out.write);

    def toXML(self, xml, prettify = False):
        """Display XML data template.
        Prettify flag tells whether to use the google code prettify markup"""
        data = {'xml':xml}
        if prettify:
            jTemplate.render("data/prettify/xml.html", data, self.response.out.write)
        else:
            jTemplate.render("data/xml.html", data,self.response.out.write );

