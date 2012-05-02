import unittest2
import sys

# sys.path.insert(0, "/Users/pamad/dev/PerseQ/")
# sys.path.insert(0, "/usr/local/google_appengine/")

from Bio import Entrez
Entrez.email = 'pamad05+entrez@gmail.com'
from xml.dom.minidom import parseString

# from google.appengine.api import memcache
from google.appengine.ext import db
from google.appengine.ext import testbed

from models.snp import snp

from models.study import Snp, Gene, Study, GWAS

from util import populate

# test get_or_insert
class DemoTestCase(unittest2.TestCase):

    def setUp(self):
        # First, create an instance of the Testbed class.
        self.testbed = testbed.Testbed()
        # Then activate the testbed, which prepares the service stubs for use.
        self.testbed.activate()
        # Next, declare which service stubs you want to use.
        self.testbed.init_datastore_v3_stub()

        # self.testbed.init_memcache_stub()

    def tearDown(self):
        self.testbed.deactivate()

    def testGetOrInsert(self):
        """- Make sure it does not create duplicates"""
        id1 = snp(key_name="test1", snpid="1805007", magnitude=0.2).put()
        act_id1 = snp.get_or_insert("test1")
        self.assertEqual(act_id1.magnitude, 0.2)
        self.assertEqual(act_id1.snpid, "1805007")
        self.assertEqual(act_id1.key(), id1)

        act_id2 = snp.get_or_insert("test2")
        self.assertEqual(act_id2.magnitude, None)

    def testdbSNP(self):
        """- Retrieve binary data from dbSNP"""
        snp = "1805007"
        res = Entrez.efetch("snp", id=snp, rettype="xml", retmode="text")
        dom = parseString(res.read())

        for node in dom.getElementsByTagName("Rs"):
            # print node.hasAttributes()
            rsid = node.getAttribute("rsId")
            #print rsid

        self.assertEqual(snp, rsid)

    def testPubMed(self):
        """- snp's list of referenced PubMed articles"""
        # articles referenced by snpid = 1805007
        referenced_articles = ['21700618', '20670983', '20585627', '20042077', '19884608', '17999355', '17952075']
        snpid = "1805007"

        handle = Entrez.elink(db="pubmed", dbfrom="snp", id=snpid, linkname="snp_pubmed_cited")
        dbs = Entrez.read(handle)
        handle.close()

        ref_ids = [] # holds the id of each of the referened articles

        # for each database with references to this snp
        for db in dbs:
            if len(db["LinkSetDb"]) == 0:
                continue
            # linkname="snp_pubmed_cited" means just one LinkSetDb - namely PubMed
            for ref in db["LinkSetDb"][0]["Link"]:
                ref_ids.append(ref['Id'])

        self.assertEqual(referenced_articles, ref_ids)

    def testGWAS(self):
        """Test relations in GWAS:
        - study(17603485) <*--*> gene(HNF1B)
        - study(17603485) <*---*> snp(4430796)"""
        populate(limit=100)
        pubmed_id = "17603485"
        snpid = "4430796"

        
        # get reference study
        stud = Study.get_by_key_name(pubmed_id)
        
        # determine the relation to the gene
        names = [gwas.gene.name for gwas in stud.gwas if gwas.gene]
        self.assertTrue("HNF1B" in names)

        # make sure the relation to study is known to gene
        gene = Gene.all().filter('name =', "HNF1B")[0]
        self.assertTrue(stud.key() in gene.studies)

        # make sure that the study is known to snp
        snp = Snp.get_by_key_name(snpid)
        self.assertTrue(stud.key() in snp.studies)


if __name__ == '__main__':
	unittest2.main()