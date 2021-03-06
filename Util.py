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
    """Add study to fulltext"""
    doc = search.Document(doc_id=study.pubmed_id, # Treat pubmed_id as key
        fields=[
            search.TextField(name='name', value=study.name),
            search.TextField(name='disease_trait', value=','.join([s.name() for s in study.diseases])),
            search.TextField(name='pubmed_id', value=study.pubmed_id),
            search.TextField(name='repl_sample', value=study.repl_sample),
            search.TextField(name='platform', value=study.platform),
            search.TextField(name='init_sample', value=study.init_sample),
            search.DateField(name='date', value=study.date)
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
        pubid = doc["PUBMEDID"].strip()
        if pubids.has_key(pubid):
            pubids[pubid].append(doc)
        else:
            pubids[pubid] = [doc]

    i = 0
    for pid, lines in pubids.iteritems():
        i += 1
        if i == 200:
            break
        #print i
        # create a new study object for each new iteration
        # - use the first line to initiate the study model
        init = lines[0]

        # create or get study with study_id
        study = Study.get_or_insert(pid,
            name=init["Study"].strip(),
            pubmed_id = pid) # disease_ref = disease,

        # populate study with static data
        # date = datetime.strptime(init["Date"].strip(), "%m/%d/%Y").date()
        study.date = datetime.strptime(init["Date"].strip(), "%m/%d/%Y").date()
        study.pubmed_url=init["Link"].strip()
        study.init_sample = init["Initial Sample Size"].strip()
        study.repl_sample= init["Replication Sample Size"].strip()
        study.platform = init["Platform [SNPs passing QC]"].strip()
        study.put()

        disease_name = None
        disease = None
        for line in lines:
            # if the disease in this GWAS row differs from the other
            # ones in this study, create new disease relation:
            tmp = line["Disease/Trait"].strip().lower()
            if disease_name != tmp:
                disease_name = tmp
                disease = Disease.get_or_insert(disease_name,
                    name=disease_name)
                study.add_disease(disease)
            
            # A gwas has either a direct gene or a 
            # down-stream and up-stream gene
            gene = None
            up_gene = None
            down_gene = None

            # The retards at GWAS use 1 == intergenic, 2 == not intergenic
            # ... no seriously, _retards_
            intergenic = line["Intergenic"].strip() == "2"

            if not intergenic:
                geneid = line["Snp_gene_ids"].strip()
                if geneid != "":
                # not intergenic => one direct gene
                    names = line["Mapped_gene"].split(" - ")
                    if len(names) == 1:
                        gene = Gene.get_or_insert(geneid,
                            name = line["Mapped_gene"],
                            geneid = geneid)
                        if not study.key() in gene.studies:
                            gene.studies.append(study.key())
                        if not disease.key() in gene.diseases:
                            gene.diseases.append(disease.key())
                        gene.put()
            else:
                # up and downstream genes must be set
                down_id = line["Downstream_gene_id"].strip()
                up_id = line["Upstream_gene_id"].strip()

                if up_id != "" and down_id != "":
                    up_down_names = line["Mapped_gene"].split(" - ")
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
            snpid = line["Snp_id_current"].strip()

            if snpid != "":
                # non=blank snp
                try:
                    # create lineation only if a 'single' snp is in the gwas
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
                except:
                    # haplotype?
                    snpid = "N/A"
            # if no gene or snp lineation is mentioned - ignore and just insert study
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
            gwas.p_string = line["p-Value"].strip()
            # could be none
            gwas.snps = snpid

            # parse out the exponent: 6E-8 => -8
            try:
                # test that p-Value is actually a float before parsing out
                float(line["p-Value"])
                gwas.p_val = int(line["p-Value"].split("E-")[1])  
            except Exception, e:
                # forces the filter to downgrade this gwas wrt. p-value
                gwas.p_val = 0
            # could be interesting
            gwas.strongest_snp_risk_allele = \
                line["Strongest SNP-Risk Allele"].strip()
            # gwas.CI = line["95% CI (text)"].strip()
            gwas.OR_beta = line["OR or beta"].strip()
            gwas.riscAlleleFrequency = \
                line["Risk Allele Frequency"].strip()
            gwas.put()

    memcache.flush_all()
    print "done"
    update_search()

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
    purge_search()

def purge_search():
    """Clear all fulltext indexes"""
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

def update_search():
    """Populate fulltext indexes"""
    for gene in Gene.all().run():
        AddGeneDocument(gene)
    for study in Study.all().run():
        AddStudyDocument(study)
    for disease in Disease.all().run():
        AddDiseaseDocument(disease)
    for snp in Snp.all().run():
        AddSNPDocument(snp)
    print "Fulltext docs updated"

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

class PJEncoder(simplejson.JSONEncoder):
    """Custom JSON encoder, since simplejson doesn't handle datetime"""
    def default(self, obj):
        if hasattr(obj, 'isoformat'):
            return obj.isoformat()
        else:
            return simplejson.JSONEncoder.default(self, obj)



class jTemplate():
    """Template helper, sets the template base path and uses a given renderer on
    the template."""

    _e = Environment(loader=FileSystemLoader(os.path.join(
        os.path.dirname(__file__), 'templates')))

    @staticmethod
    def render(template, variables, printer):
        t = jTemplate._e.get_template(template)
        printer(t.render(variables))

class AppRequestHandler(webapp2.RequestHandler):
    """Base class for controllers"""
    _template = None

    def render(self, template = None, **args):
        if template is None:
            raise Exception("No 'template' argument to render from")
        temp = jTemplate._e.get_template(template)
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

        temp = jTemplate._e.get_template(self._template)
        self.response.write(temp.render(dictionary))

    def toJson(self, dictionary, prettify = False):
        """Display JSON data template.
        Prettify flag tells whether to use the google code prettify markup"""
        enc = PJEncoder()
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

