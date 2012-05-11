"""Annotation models. Models for user feedback like comments and ratings"""

from google.appengine.ext import db
from models.study import Study, Snp,GWAS # For referencing
from models.users import UserData

class Comment(db.Model):
    """Comments, associated with either studies, snps or gwas entries"""
    user = db.ReferenceProperty(UserData)
    date = db.DateTimeProperty()
    body = db.StringProperty()
    study = db.ReferenceProperty(Study,
                collection_name='comments')
    snp = db.ReferenceProperty(Snp,
                collection_name='comments')
    gwas = db.ReferenceProperty(GWAS,
                collection_name='comments')

