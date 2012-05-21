""" Controller for managing the different kinds of search """
from Util import AppRequestHandler
from models.study import Study, Snp, Gene, Disease

class SearchHandler(AppRequestHandler):
    _model = None
    """ Handle SNP searches """
    def get(self, query):
        """ Request via GET = single """
        result = None
        if query != "":
            result = self._model.search(query)
        self.response.out.write(result)
        #self.toJson(result)

    #def post(self, query):
    #    """Request via POST = Batch"""
    #    pass


class SNPSearch(SearchHandler):
    _model = Snp

class GeneSearch(SearchHandler):
    _model = Gene

class DiseaseSearch(SearchHandler):
    _model = Disease

class StudySearch(SearchHandler):
    _model = Study

__routes__ = [('/search/snp/(.*)', SNPSearch),
              ('/search/gene/(.*)', GeneSearch),
              ('/search/disease/(.*)', DiseaseSearch),
              ('/search/study/(.*)', StudySearch)]
