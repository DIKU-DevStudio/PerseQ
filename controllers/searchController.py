""" Controller for managing the different kinds of search """
from Util import AppRequestHandler
from models.study import Study, Snp, Gene, Disease

class SNPSearch(AppRequestHandler):
    """ Handle SNP searches """
    def get(self, query):
        """ Request via GET = single """
        snp = None
        if query != "":
            snp = Snp.gql('WHERE snpid = :1', query)

        if snp.get():
            ret = {'count':snp.count(),'result':[snp.get()]}
            ret['result'][0]['gene'] = snp.get().gene
        else:
            ret = {'count':0,'result':[]}
        self.toJson(ret)

    def post(self, query):
        """ Request via POST = Batch """
        pass


class GeneSearch(AppRequestHandler):
    """ Handle Gene searches """
    def get(self, query):
        """ Request via GET = single """
        result = None
        if query != "":
            result = Gene.gql('WHERE :1 IN (geneid, name, alias)', query)

        self.toJson(result)

    def post(self, query):
        """ Request via POST = Batch """
        pass

class DiseaseSearch(AppRequestHandler):
    """ Handle Gene searches """
    def get(self, query):
        """ Request via GET = single """
        result = None
        if query != "":
            result = Gene.gql('WHERE :1 IN (geneid, name, alias)', query)
        self.toJson(result)

    def post(self, query):
        """ Request via POST = Batch """
        pass

class StudySearch(AppRequestHandler):
    """ Handle Gene searches """
    def get(self, query):
        """ Request via GET = single """
        result = None
        if query != "":
            result = Gene.gql('WHERE :1 IN (geneid, name, alias)', query)
        self.toJson(result)

    def post(self, query):
        """ Request via POST = Batch """
        pass


#class MainPage(webapp.RequestHandler):
#    """Handles search requests for comments."""
#
#    def get(self):
#        """Handles a get request with a query."""
#        uri = urlparse(self.request.uri)
#        query = ''
#        if uri.query:
#            query = parse_qs(uri.query)
#            query = query['query'][0]
#
#        # sort results by author descending
#        expr_list = [search.SortExpression(
#            expression='author', default_value='',
#            direction=search.SortExpression.DESCENDING)]
#        # construct the sort options 
#        sort_opts = search.SortOptions(
#             expressions=expr_list)
#        query_options = search.QueryOptions(
#            limit=3,
#            sort_options=sort_opts)
#        query_obj = search.Query(query_string=query, options=query_options)
#        results = search.Index(name=_INDEX_NAME).search(query=query_obj)
#
#        if users.get_current_user():
#            url = users.create_logout_url(self.request.uri)
#            url_linktext = 'Logout'
#        else:
#            url = users.create_login_url(self.request.uri)
#            url_linktext = 'Login'
#
#        template_values = {
#            'results': results,
#            'number_returned': len(results.results),
#            'url': url,
#            'url_linktext': url_linktext,
#        }
#
#        path = os.path.join(os.path.dirname(__file__), 'index.html')
#        self.response.out.write(template.render(path, template_values))
#
#


__routes__ = [('/search/snp/(.*)', SNPSearch),
              ('/search/gene/(.*)', GeneSearch),
              ('/search/disease/(.*)', DiseaseSearch),
              ('/search/study/(.*)', StudySearch)]
