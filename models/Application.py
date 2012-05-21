from google.appengine.ext import db
from google.appengine.api import search

class AppModel(db.Model):
    _index = None

    """ General methods for models """
    def __dict__(self):
        """ Allow serialization of model objects """
        return dict([(p, unicode(getattr(self, p))) for p in self.properties()])

    @classmethod
    def search(self,query):
        """Allow full-text search"""
        if self._index:
            search.Index(name=self._index).search('"'+query+'"')
        else:
            raise AttributeError('Index is not set for search on '+str(self))
