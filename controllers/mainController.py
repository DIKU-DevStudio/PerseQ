###
# Controller for main presentation
#
###
from Utilities import AppRequestHandler
from models.snp import snp


class snpSearch(AppRequestHandler):
    def get(self):
        snpids = snp.all()
        snpids.order("snpid")
        self.out({'snps':snpids})

__routes__ = [('/', snpSearch)]
