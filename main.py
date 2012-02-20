#!/usr/bin/env python
#
# Copyright 2007 Google Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
from google.appengine.ext import webapp
from google.appengine.ext.webapp import util
from wikitools import wiki, page, api
import pprint
pp = pprint.PrettyPrinter(indent=4)



class LookUpSNP(webapp.RequestHandler):
    def get(self):
        snp = self.request.get("snp")
        if len(snp) == 0:
            self.response.out.write("No SNP title given (snp='')")
            return

        site = wiki.Wiki("http://bots.snpedia.com/api.php")
        params = {
            'action': 'query',
            'prop': 'revisions',
            'rvprop': 'content|timestamp',
            'rvlimit': '1',
            'titles':snp
        }
        req = api.APIRequest(site, params)
        result = req.query(querycontinue=False)

        # if pageid == -1 <=> title=snp does not exist
        pageid = int(result['query']['pages'].keys()[0])
        if pageid == -1:
            self.response.out.write("SNP title does not exist (NoPage)")
            return
        
        # should use the proper respone.out.write(result)
        # but for some reason, right now it looks like crap..
        pp.pprint(result['query']['pages'][str(pageid)]['revisions'][0]['*'].encode('utf-8'))
    	# self.response.out.write(result)

def main():
    application = webapp.WSGIApplication([('/', LookUpSNP)], debug=True)
    util.run_wsgi_app(application)


if __name__ == '__main__':
    main()
