from google.appengine.ext import webapp
from jinja2 import Environment, FileSystemLoader
import os
import json
import inspect

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
            actionName = get_class_from_frame(inspect.stack()[1][0]).__name__
            ctrlName = inspect.stack()[1][1].split('/')[-1].split('Controller.py')[0]
            self._template = ctrlName+"/"+actionName+".html"
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

def get_class_from_frame(fr):
    args, _, _, value_dict = inspect.getargvalues(fr)
    # we check the first parameter for the frame function is
    # named 'self'
    if len(args) and args[0] == 'self':
        # in that case, 'self' will be referenced in value_dict
        instance = value_dict.get('self', None)
        if instance:
            # return its class
            return getattr(instance, '__class__', None)
    # return None otherwise
    return None
