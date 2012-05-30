from google.appengine.ext import db
from google.appengine.api import search

class AppModel(db.Model):
    _index = None

    """ General methods for models """
    def __dict__(self):
        """ Allow serialization of model objects """
        return dict([(p, unicode(getattr(self, p))) for p in self.properties()])

    @classmethod
    def search(self, query):
        """Allow full-text search. Return search.SearchResults object"""
        if self._index:
            return search.Index(name=self._index).search(query)
        else:
            raise AttributeError('Index is not set for search on '+str(self))

    @classmethod
    def search_todict(self, query):
        """Full-text search returning objects as dict"""
        search = self.search(query)
        results = []
        for doc in search.results:
            item = {}
            for field in doc.fields:
                item[field.name] = field.value
            results.append(item)
        return results

