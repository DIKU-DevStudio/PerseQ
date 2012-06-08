'''
he3.db.tools.mappers
contains mappers written to work with App Engine MapReduce
http://code.google.com/p/appengine-mapreduce/

Created on Jun 27, 2010
@author: Ben Davies, Helium 3 IT Solutions
'''

import google.appengine.ext.db as db
import mapreduce.operation as op

class Mapper (object):
	'''A standard base class for Mappers defined here'''

	@staticmethod
	def process(data):
		'''Standard method call to invoke mapper on mapped data. Intended to 
		be overriden.
		'''
		pass
	
class ModelHygieneMapper (Mapper):

	
	'''A datastore mapper for performing common maintenance tasks on
	datastore model entities.'''
	
	@staticmethod
	def process(entity):
		'''Checks and repairs model integrity of the passed entity 
		1. Removes dangling references
		2. Sets undefined datastore values to the default
		'''
		
		props = [x for x in entity.__class__.__dict__.values()\
				if isinstance(x, db.Property)]
		changed = False
		
		for prop in props:
			if not prop.get_value_for_datastore(entity):
				prop.__set__(entity, prop.default_value())
				changed = True		
			elif isinstance(prop, db.ReferenceProperty):
				if not db.get(prop.get_value_for_datastore(entity)):
					prop.__set__(entity, None)
					changed = True
						
		if changed: yield op.db.Put(entity)
				