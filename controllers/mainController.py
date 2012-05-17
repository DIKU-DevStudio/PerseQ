"""
Controller for main presentation
"""
from util import AppRequestHandler
from models.snp import snp
from models.users import UserData

class snpSearch(AppRequestHandler):
    """Search handler"""
    def get(self):
        snpids = snp.all()
        snpids.order("snpid")
        self.out(snps=snpids, search=self.request.get("q"))

class dashboard(AppRequestHandler):
    """Dashboard data fetching"""
    def get(self):
        user = UserData.current()
        self.out(username=user.nickname)

__routes__ = [('/',dashboard),
              ('/search/', snpSearch)]
