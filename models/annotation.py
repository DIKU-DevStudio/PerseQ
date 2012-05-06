from google.appengine.ext import db
from models.study import Study, Snp # For referencing

class StudyComment(db.Model):
    user = db.UserProperty()
    date = db.DateTimeProperty()
    body = db.StringProperty()
    study = db.ReferenceProperty(Study,
                collection_name='comments')


class SnpComment(db.Model):
    user = db.UserProperty()
    date = db.DateTimeProperty()

    snp = db.ReferenceProperty(Snp,
                collection_name='comments')
