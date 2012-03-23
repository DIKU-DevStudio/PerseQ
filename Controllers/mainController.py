###
# Controller for main presentation
#
###
import Utilities
from Models.snp import snp

class snpsList(Utilities.AppRequestHandler):
    def get(self):
        self.setTemplate('Main/snpsearch.html')
        snpids = snp.all()
        snpids.order("snpid")
        self.out({'snps':snpids})
