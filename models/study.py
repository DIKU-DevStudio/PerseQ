from google.appengine.ext import db

# PK = pubmed_id
class Study(db.Model):
    # 3
    date = db.DateProperty()
    # 1 - private key
    pubmed_id = db.IntegerProperty()
    # 5
    pubmed_url = db.LinkProperty()
    # 6 - description of study
    name = db.StringProperty(required=True)
    # 7 disease or trait researched
    disease_trait = db.StringProperty(required=True)

    abstract = db.StringProperty()
    # 8 - antal forsogspersoner
    init_sample = db.StringProperty()
    # 9
    repl_sample = db.StringProperty()
    platform = db.StringProperty()

# id NCBI gene_id
# handle relations with ancestors
class Gene(db.Model):
    studies = db.ListProperty(db.Key) # or simply pubmed_ids..
    alias = db.StringListProperty()
    name = db.StringProperty(required=True)
    geneid = db.IntegerProperty()

# id = SNPID
class Snp(db.Model):
    studies = db.ListProperty(db.Key)
    # kan ligge uden for et gen, derfor ikke required
    gene = db.ReferenceProperty(Gene)
    snpid = db.StringProperty() # rs1805007

# ID = random
# ancestor=study_id == pubmed_id
class GWAS(db.Model):
    # maybe ancestor instead?
    # - improves consistency in high replication data store
    study = db.ReferenceProperty(Study, required=True)

    # if the snp is _in_ a specific gene - this is the id
    gene = db.ReferenceProperty(Gene)
    # if not, these two contain the reference-ids
    upstream = db.StringProperty()
    downstream = db.StringProperty()
    # upstream = db.ReferenceProperty(Gene, )
    # downstream = db.ReferenceProperty(Gene)

    # snp ids.. 
    snps = db.ListProperty(db.Key)
    # 27 - # 6 * 10^-8
    p_string = db.StringProperty()
    # 6 * 10^-8 => p_val = -8
    p_val = db.IntegerProperty() 
    # 30 - Odds ratio #.##
    OR_beta = db.FloatProperty()
    # 26 - #.##
    riscAlleleFrequency = db.FloatProperty()
    # 25 - 1=no, 2=yes
    intergenic = db.BooleanProperty()