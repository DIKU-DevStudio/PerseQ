'''
Custom Datastore Properties for handling dates
'''

#imports
from pytz.gae import pytz
from google.appengine.ext import db

class UtcDateTimeProperty(db.DateTimeProperty):
    '''Marks DateTimeProperty values returned from the datastore as UTC. Ensures
    all values destined for the datastore are converted to UTC if marked with an
    alternate Timezone.
    
    Inspired by 
    http://www.letsyouandhimfight.com/2008/04/12/time-zones-in-google-app-engine/
    http://code.google.com/appengine/articles/extending_models.html
    '''
    
    
    def get_value_for_datastore(self, model_instance):
        '''Returns the value for writing  to the datastore. If value is None,
        return None, else ensure date is converted to UTC. Note Google App 
        Engine already does this. Called by datastore		
        '''
        
        date = super(UtcDateTimeProperty, self)\
            .get_value_for_datastore(model_instance)
        if date:
            if date.tzinfo:
                return date.astimezone(pytz.utc)
            else:
                return date.replace(tzinfo=pytz.utc)
            return None
        
        
    def make_value_from_datastore(self, value):
        '''Returns the value retrieved from the datastore. Ensures all dates
        are properly marked as UTC if not None
        '''
        
        if value is None:
            return None
        else:
            return value.replace(tzinfo=pytz.utc)

