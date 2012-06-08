'''
This module contains classes for improving query performance
By Ben Davies (including source from Nick Johnson and others)
code.google.com/p/he3-appengine-lib for updates and further links 
'''
import google.appengine.ext.db as db

namespace = 'he3'

class PrefetchingQuery(object):
	'''
	This class is a facade to a db.Query and db.GqlQuery object that offers 
	additional functionality to automatically prefetch reference properties
	and parents across a set of fetch results such that individual get() calls 
	for the same datastore object are minimised.  
		
	The idea and technique for prefetching reference properties came from a blog
	post (http://blog.notdot.net/2010/01/ReferenceProperty-prefetching-in-App-Engine)
	by Nick Johnson of Google. This class is little more than a repackaging of the
	concept and code provided by Nick for prefetching property values. 
	
	!WARNING! - Prefetching of Parent Props
	Please note that in order to set the parent value in entities this code uses
	the internal/private model instance _parent. This is an undocumented and 
	unsupported approach to setting the parent. 
	UPDATES TO THE SDK AND/OR APP ENGINE PLATFORM MAY BREAK THIS CODE IN THE 
	FUTURE  
	
	!WARNING! - Shared reference entities
	This code may (and in the most efficient case, definitely will) cause 
	reference properties across instances to reference the same entity. Care 
	should be taken when modifying prefetched entities.
	
	Prop* to Nick! (no pun intended).
	 	
	USAGE:
	
	Instantiate a PrefetchingQuery with an existing db.Query or db.GqlQuery and
	an optional list of properties to prefetch: 
	
	myPrefetchingQuery = PrefetchedQuery(myEntity.all(), (prop1, prop2))
	
	When you call fetch() on the PrefetchingQuery, the reference properties you
	have specified or are used by default are prefetched for each entity being 
	returned by fetch().

	Note that PrefetchingQuery can also prefetch the entity's parent, if that 
	is useful. Simply pass the string 'parent' in the list of reference 
	properties to prefetch. 
		
	You can specify reference properties to prefetch using PrefetchingQuery in
	3 ways (in the following order of precedence):
	
	1. In the constructor, as above.
	
	2. In a class attribute, properties_to_prefetch
		Eg. MyModelEntity.properties_to_prefetch =
		(MyModelEntity.ref_prop1, MyModelEntity.ref_prop2)
	
	3. Not at all, in which case all reference properties (including the 
	parent) will be prefetched.
	
	PrefetchingQuery supports the methods and functionality of the underlying 
	Query object it was passed in the constructor. For example, you can use
	the order(), filter() and ancestor() methods to further refine your query
	(now wrapped in the a PrefetchingQuery), but only if the PrefetchingQuery
	was initialised with a db.Query, not a db.GqlQuery.
	'''
	
	class_property_name = 'properties_to_prefetch'

	def __init__(self, query, properties_to_prefetch = None):
		'''
		Constructor for a PrefetchingQuery.
		@param query: a google.appengine.ext.db.query or db.GqlQuery object
		@param properties_to_prefetch: a list of reference property names
			defining which properties to prefetch on a fetch() call. If not
			supplied, a class attribute is checked. If not present, properties
			to dereference are automatically determined. 
		
		@raise TypeError: raised if query is not an instance of db.Query or 
		db.GqlQuery 
		'''
		
		self._query = query
		
		self._properties_to_prefetch = properties_to_prefetch
		
		if isinstance(query, db.Query): self._query_type = 'Query'
		elif isinstance(query, db.GqlQuery): self._query_type = 'GqlQuery'
		else: raise TypeError('Query type not supported: '\
			 + type(query).__name__)

	def fetch(self, limit, offset=0):
		''' executes query against datastore as per db.Query.fetch(). Afterwards,
		performs reference property prefetching before returning results.
		
		@param limit: Maximum amount of results to retrieve as per 
		db.Query.fetch()
		@param offset: Number of results to skip prior to returning resultset.
		As per db.Query.fetch().
		
		@return: A list of entity results, as per db.Query.fetch()

		NOTE: this method should match the corresponding signature of 
		db.Query.fetch() precisely.
		@see: http://code.google.com/appengine/docs/python/datastore/queryclass.html
		'''

		entities = self._query.fetch(limit,offset)

		if len(entities):
			
			props = (self.properties_to_prefetch 
				or PrefetchingQuery._get_properties_defined_in_class(entities[0])
				or PrefetchingQuery._automatically_determine_refprops(entities[0]))
			
			self._prefetch_refprops(entities, props)
		
		return entities
	
	def filter(self, property_operator, value):
		'''Adds a property condition filter to the query. Only entities with
		properties that meet all of the conditions will be returned by the 
		query. This method should behave identically to the db.Query.filter()
		method. Using this method also clears any caching of the object.
		@attention: This method is only available for Queries used
		to initalise the DereferencedQuery of type db.Query
		@see: http://code.google.com/appengine/docs/python/datastore/queryclass.html
		
		@param property_operator: A string containing the property name, and an 
		optional comparison operator
		@param value: The value to use in the comparison on the right-hand side
		of the expression
		@return: The query with filter added
		@raise TypeError: raised if the query not the correct type
		'''
		self._check_query_type_is('Query')
		self._query = self._query.filter(property_operator, value)
		return self 
		
	
	def order(self, property):
		'''Adds an ordering for the results. Results are ordered starting with
		the first order added. This method should behave identically to the 
		db.Query.order() method. Using this method also clears any caching of 
		the object.
		
		@attention: This method is only available for Queries used
		to initalise the DereferencedQuery of type db.Query
		@see: http://code.google.com/appengine/docs/python/datastore/queryclass.html
		
		@param property: A string, the name of the property to order
		@return: The query with order added
		@raise TypeError: raised if the query not the correct type
		'''
		self._check_query_type_is('Query')
		self._query.order(property)
		return self
	
	def ancestor(self, ancestor):
		'''Adds an ancestor condition filter to the query. Only entities with
		the given entity as an ancestor (anywhere in its path) will be returned 
		by the query. This method should behave identically to the 
		db.Query.ancestor() method. Using this method also clears any caching of 
		the object.
		
		@attention: This method is only available for Queries used
		to initalise the DereferencedQuery of type db.Query
		@see: http://code.google.com/appengine/docs/python/datastore/queryclass.html
		
		@param ancestor: A Model instance or Key instance representing the 
		ancestor.
		@return: Itself after ancestor condition has been added
		@raise TypeError: raised if the query not the correct type
		'''
		self._check_query_type_is('Query')
		self._query.ancestor(ancestor)
		return self
	
	def count(self, limit=1000):
		'''Returns the number of results this query fetches. This method should
		behave identically to the method of the same name of db.Query and 
		db.GqlQuery
		
		@see: http://code.google.com/appengine/docs/python/datastore/queryclass.html
		
		@param limit: The maximum number of results to count.
		@return: Returns the number of result this query fetches
		'''
		return self._query.count(limit)		

	def with_cursor(self, cursor):
		'''Sets a cursors to use for fetches. This method should behave 
		identically to the one existing on db.Query and db.GqlQuery objects.
		
		@see: http://code.google.com/appengine/docs/python/datastore/queryclass.html
		
		@param cursor: Cursor to set.
		@return: None
		'''
		return self._query.with_cursor(cursor)
			
	def cursor(self):
		'''Returns a cursor pointing to the result after the last result fetched.
		This method should behave identically to the one existing on db.Query 
		and db.GqlQuery objects.
		
		@see: http://code.google.com/appengine/docs/python/datastore/queryclass.html

		@return: base-64 encoded string cursor
		'''
		return self._query.cursor()

	def _check_query_type_is(self, required_query_type):
		'''This is a helper method to assert that the query the DereferencedQuery
		was initialised with is of the correct type.
		
		@param required_query_type: Value of self._query_type expected (
		currently only 'Query' or 'GqlQuery')
		@return: nothing 
		@raise TypeError: raised if the query not the correct type
		'''
		
		if self._query_type != required_query_type:
			raise TypeError('Operation not allowed for query type ('\
				+ type(self._query).__name__)		
		
	def _get_properties_to_prefetch(self):
		'''Returns the list of properties to dereference on a fetch() call'''
		return self._properties_to_prefetch
	
	def _set_properties_to_prefetch(self, list_of_prop_names):
		'''Sets the list of property names to derference on a fetch() call. Type
		checking only is done here'''

		if type(list_of_prop_names) <> tuple:
			raise TypeError('Invalid list_of_prop_names parameter. Must be a list')
		self._properties_to_prefetch = list_of_prop_names
	
	def _prefetch_refprops(self, entities, props):
		'''This entire function almost identical to Nick Johnsons blog post 
		referenced above. I need to think *really* hard to follow it. 
		
		Changes from Nicks Code: 
		1. Filtered none values from set of keys passed to db.get(). There was an
			an alternative implementation of this using a filter()ed ref_keys set
			and using the original in the zip operation, talked about in the comments
			of Nicks blog entry.
			(impacts optional properties)
		2. used a technique by Ubaldo Huerta to include parent keys
			see http://groups.google.com/group/google-appengine-python/msg/22c2010a8f102f32
		3. Filtered none values from the set of referential entities returned
			from db.Get(). Skip populating those fields where an entity was not
			returned. If an application does not clean up dangling references
			this would otherwise cause errors at 'x.key()' 
		'''
		fields = [(entity, prop) for entity in entities for prop 
				in props]
		ref_keys = [x.key().parent() if prop == 'parent' else 
				prop.get_value_for_datastore(x) for x, prop in fields]

		ref_entities = dict((x.key(), x) for x in db.get(set(ref_keys)-set((None,)))
						if x is not None)
		for (entity, prop), ref_key in zip(fields, ref_keys): 
			if ref_entities.has_key(ref_key):
				if prop == 'parent': 
					# Big warning ! Using internals of Model (might	break in the future)
					if ref_key:entity._parent = ref_entities[ref_key] 
				else:
					if ref_key:prop.__set__(entity, ref_entities[ref_key])
			else:
				#We couldn't retrieve a referential entity for the current 
				#entity,prop pair. This can happen if a App Engine application
				#deleted entities without cleaning up the entities that reference
				#them. This why simply testing referential entities on retrieval,
				#without purposefully cleaning up dangling references, sucks.
				#</rant>
				pass 
		return entities
	
	@staticmethod
	def _get_properties_defined_in_class(entity_instance):
		'''Returns the properties to prefetch that are defined on the class of
		the provided entity, if present, or returns None'''
		
		class_dict = entity_instance.__class__.__dict__
		class_prop_name = PrefetchingQuery.class_property_name
		
		if class_dict.has_key(class_prop_name): return class_dict[class_prop_name]
		else: return None
	
	@staticmethod
	def _automatically_determine_refprops(entity_instance):
		'''Given an entity instance determine which properties are refprops
		or otherwise available for prefetching by default'''
		
		ref_props = [x for x in entity_instance.__class__.__dict__.values()\
				if isinstance(x, db.ReferenceProperty)]
		ref_props.append('parent')

		return ref_props
		
			
	properties_to_prefetch = property(fget=_get_properties_to_prefetch
										, fset=_set_properties_to_prefetch
										, doc='Properties to dereference on fetch()')			