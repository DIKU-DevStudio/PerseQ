from google.appengine.ext import webapp
from Template import jTemplate

class AppRequestHandler(webapp.RequestHandler):
    _template = "sample.html"

    def setTemplate(self, template):
        self._template = template

    def out(self, dictionary = {}, template = None):
        if(template == None):
            template = self._template
        jTemplate.render(template, dictionary, self.response.out.write)

