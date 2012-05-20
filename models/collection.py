from google.appengine.ext import db

# from models.Application import AppModel

# from models.users import UserData

class SNPCollection(db.Model):
	"""SNPCollection is a simple list of SNPs associated with a user
	- each collection-name is unique for that user
	- the list is allowed to be empty"""
	# we do not require the name be unique globally
	name = db.StringProperty(required=True)
	# list of snp-ids
	snps = db.ListProperty(db.Key)
	count = db.IntegerProperty(default=0)

	@property
	def addMultipleSNP(self, snps=None):
		"""Argument: takes a list of keys for SNPs and adds them to the collection
		- making sure not to add duplicates
		- updates the database afterwards"""
		if snps is None:
			return

		# snps must be a list of keys
		if not isinstance(snps, collections.Iterable):
			raise Exception("Must be a list of db.Keys")
			return 

		for snp in snps:
			if snp is not NoneType and not snp in self.snps:
				self.snps.append(snp)
		# self.snps.extend(snps)
		# update in db
		self.count = len(self.snps)
		self.put()

		return


	@property
	def addSNP(self, snp=None):
		"""Adds a _single_ Snp-key to the list of SNPs"""
		if snp is None:
			return None

		if not isinstance(snp, db.Key):
			# logging.warning("addSNP received a non-db.Key instance ")
			raise Exception("argument 'snp' must be an instance of type db.Key")

		if snp in self.snps:
			return None
		
		self.snps.append(snp)
		self.count = len(self.snps)
		self.put()

	@property
	def removeSNP(self, key):
		"""Removes the given SNP from the collection (if it is in the collection)"""
		if not isinstance(key, db.Key):
			raise Exception("key argument must be an instance of 'db.Key'")

		try:
			self.snps.remove(snp_key)
			self.count -= 1
		except:
			return False

		return True


	@staticmethod
	def newCollection(name=None, user=None, snps=None):
		"""creates a new SNPCollection with <name> for a given user

		ARGUMENTS:
		- <user> and <name> are required arguments
		- snps must be a list of SNP.key()'s, or None

		RETURNS:
		- SNPCollection instance upon success
		- None if an error occurred
		"""
		if user is None:
			return None

		if name is None or name == "":
			return None

		if snps is not None and isinstance(snps, collections.Iterable):
			# if snps is a list remember to update count when creating object
			return SNPCollection.get_or_insert(parent=user, name=name,
			snps=snps, count=len(snps))

		return SNPCollection.get_or_insert(parent=user, name=name)