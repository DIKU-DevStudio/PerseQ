###
# Controller for administrator presentation
#
###
from util import AppRequestHandler

class adminDashboard(AppRequestHandler):
    def get(self):
        self.setTemplate('dashboard.html')
        self.out()


__routes__ = [('/admin/',adminDashboard)]
