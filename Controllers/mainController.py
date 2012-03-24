###
# Controller for main presentation
#
###
from Utilities import AppRequestHandler
from Models.snp import snp


class snpsList(AppRequestHandler):
    def get(self):
        self.setTemplate('Main/snpsearch.html')
        snpids = snp.all()
        snpids.order("snpid")
        self.out({'snps':snpids})

__routes__ = [('/', snpsList)]
