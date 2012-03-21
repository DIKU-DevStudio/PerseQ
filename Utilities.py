from google.appengine.ext import webapp
from Template import jTemplate

from Bio import Entrez

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

class AppRequestHandler(webapp.RequestHandler):
    _template = "sample.html"

    def setTemplate(self, template):
        self._template = template

    def out(self, dictionary = {}, template = None):
        if(template == None):
            template = self._template
        jTemplate.render(template, dictionary, self.response.out.write)

