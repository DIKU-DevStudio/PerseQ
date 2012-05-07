from google.appengine.ext import db

class Disease(db.Model):
    name = db.StringProperty(required=True)

# PK = pubmed_id
class Study(db.Model):
    # 3
    date = db.DateProperty()
    # 1 - private key
    pubmed_id = db.StringProperty(required=True)
    # 5
    pubmed_url = db.LinkProperty()
    # 6 - description of study
    name = db.StringProperty(required=True)
    # 7 disease or trait researched
    disease_ref = db.ReferenceProperty(Disease, required=True)
    # to prevent relation traversal
    disease_trait = db.StringProperty(required=True)

    abstract = db.StringProperty()
    # 8 - antal forsogspersoner
    init_sample = db.StringProperty()
    # 9
    repl_sample = db.StringProperty()
    platform = db.StringProperty()

    @property
    def genes(self):
        return Gene.gql("WHERE studies = :1", self.key())

# id NCBI gene_id
# handle relations with ancestors
class Gene(db.Model):
    studies = db.ListProperty(db.Key) # keys of studies
    diseases = db.ListProperty(db.Key) # keys of diseases
    alias = db.StringListProperty()
    name = db.StringProperty()
    geneid = db.StringProperty(required=True)

# id = SNPID
class Snp(db.Model):
    studies = db.ListProperty(db.Key)
    diseases = db.ListProperty(db.Key)
    # kan ligge uden for et gen, derfor ikke required
    gene = db.ReferenceProperty(Gene,
            collection_name='snps')
    snpid = db.StringProperty() # rs1805007

# ID = random
# ancestor=study_id == pubmed_id
class GWAS(db.Model):
    # maybe ancestor instead?
    # - improves consistency in high replication data store
    study = db.ReferenceProperty(Study, required=True)

    disease = db.ReferenceProperty(Disease, required=True)

    # if the snp is _in_ a specific gene - this is the id
    gene = db.ReferenceProperty(Gene, collection_name="gene")
    # if not, these two contain the reference-ids
    upstream = db.ReferenceProperty(Gene, collection_name="upstream")
    downstream = db.ReferenceProperty(Gene, collection_name="downstream")

    # snp ids..
    snps = db.StringProperty()
    # 27 - # 6 * 10^-8
    p_string = db.StringProperty()
    # 6 * 10^-8 => p_val = -8
    p_val = db.IntegerProperty() 
    # 30 - Odds ratio #.##
    # OR_beta = db.FloatProperty()
    OR_beta = db.StringProperty()

    CI = db.StringProperty()
    # 26 - #.##
    riscAlleleFrequency = db.StringProperty()
    # Strongest SNP-Risk Allele
    strongest_snp_risk_allele = db.StringProperty()
    # 25 - 1=no, 2=yes
    intergenic = db.BooleanProperty()
