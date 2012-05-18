from google.appengine.ext import db

class AppModel(db.Model):
    """ General methods for models """
    def __dict__(self):
        """ Allow serialization of model objects """
        return dict([(p, unicode(getattr(self, p))) for p in self.properties()])
