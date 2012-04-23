###
# Controller for main presentation
#
###
from Utilities import AppRequestHandler
import csv
from StringIO import StringIO

class gwasReader(AppRequestHandler):
    def get(self):
        self.setTemplate('gwas.html')
        reader = csv.DictReader(open('gwascatalog.txt','rb'), dialect='excel-tab')
        text = []
        for row in reader:
            text.append(row)

        self.out({'gwasList': text[:100]})

__routes__ = [('/gwas/', gwasReader)]
