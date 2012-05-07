###
# Controller for main presentation
#
###
from util import AppRequestHandler
from models.snp import snp
from models.users import UserData

class snpSearch(AppRequestHandler):
    def get(self):
        snpids = snp.all()
        snpids.order("snpid")
        self.out({'snps':snpids, 'search':self.request.get("q")})

class dashboard(AppRequestHandler):
    def get(self):
        #TODO: If we want to use a login link on the dashboard, this is google's suggested code:
        # user = users.get_current_user()
        # if user:
        #     greeting = ("Welcome, %s! (<a href=\"%s\">sign out</a>)" %
        #                 (user.nickname(), users.create_logout_url("/")))
        # else:
        #     greeting = ("<a href=\"%s\">Sign in or register</a>." %
        #                 users.create_login_url("/"))
        user = UserData.current()
        self.out()
        
        

__routes__ = [('/',dashboard),
              ('/search/', snpSearch)]
