'''
The reference module contains properties that relate to other datastore
entities
'''
from google.appengine.ext import db
from google.appengine.api import datastore_errors
from google.appengine.api import datastore

BadValueError = datastore_errors.BadValueError #pylint:disable=C0103

class ReferenceListProperty(db.ListProperty):
    '''
    The ReferenceListProperty allows a many-to-many relationship between
    2 entity types. As a ListProperty subclass, you interact with the 
    ReferenceListProperty as a list.
    
    ReferenceListProperty, like the ReferenceProperty, creates a corresponding
    collection on the side of the referenced entity
    
    
    CURRENT STATUS
    At time of writing
    - Reverse Collection Property is not being created
    - Referenced entities are always retrieved from the datastore. This should
        be replaced by model-instance caching and possibly even lazy loading so
        individual list members can be instanced as required.
    
    
    '''

    def __init__(self, reference_class=None, verbose_name=None, 
                 collection_name=None, **attrs):
        '''
        Construct ReferenceListProperty
        
        args:
            reference_class is a db.model class to create a list of references
            of. As per ReferenceProperty
            verbose_name is a human readable name
            collection_name is not used (yet!)
        
        raises:
            KindError if the reference class is not a db.model object
            
        '''      
        super(ReferenceListProperty, self).__init__(item_type=db.Key, 
                                                    verbose_name=verbose_name,
                                                    **attrs)
    
        self.collection_name = collection_name
    
        if reference_class is None:
            reference_class = db.Model
        if not (isinstance(reference_class, type) and 
                issubclass(reference_class, db.Model)):
            raise db.KindError('reference_class must be db.Model')
        self.reference_class = self.data_type = reference_class

    def validate(self, value):
        '''
        Validate the ReferenceListProperty values as a list of model instance
        of class reference class that have been saved.

        '''

        if value is not None:
            if not isinstance(value, list):
                raise datastore_errors.BadValueError(
                    'Property %s must be a list' % self.name)
            
            value = self.validate_list_contents(value)
            
            #check required, choices, user validator
            value = super(ReferenceListProperty, self).validate(value)
        
        return value
        
    def validate_list_contents(self, value):
        '''
        validates each member of the list as being not None, not a member of 
        any class except reference_class, and being saved. 
        '''
        for member in value:
            
            if isinstance(member, datastore.Key):
                #skip further checks
                continue
            
            if isinstance(member, db.Model) and not member.has_key():
                raise BadValueError(
                    '%s instance must have a complete key before it can be '
                    'stored as a reference' % self.reference_class.kind())
            
            if not isinstance(member, self.reference_class):
                raise db.KindError('All members of %s must be an instance of %s'
                                   % (self.name, self.reference_class.kind())) 
                
        return value

    def get_value_for_datastore(self, model_instance):
        '''
        Returns a list of keys for storage in the datastore
        '''
        entity_list = self.validate_list_contents(
            self.__get__(model_instance, self.reference_class))

        return [e.key() for e in entity_list if e]
    
    def make_value_from_datastore(self, value):
        '''
        Returns a populated list of entity models
        '''
        return db.get(value);

