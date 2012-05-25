###
# Controller for administrator presentation
#
###
from Util import AppRequestHandler

class adminDashboard(AppRequestHandler):
    def get(self):
        self.setTemplate('dashboard.html')
        self.out()


__routes__ = [('/admin/',adminDashboard)]
