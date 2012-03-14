from models import domain_tag

tags = [
    "Oncology",
    "Dermatology",
    "Neurology",
    "Unknown",
        ]
def init():
    for t in tags:
        thistag = domain_tag(tag=t)
        print 'adding ' + t
        thistag.put()

