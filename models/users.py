from google.appengine.api import datastore_errors, users
from google.appengine.ext import db

class UserData(db.Model):
    # derived partly from https://khanacademy.kilnhg.com/Code/Website/Group/stable/File/user_models.py?rev=tip


    # do not reference this directly
    user = db.UserProperty()

    # This should be treated as primary key - we populate automatically from Google user ID
    user_id = db.StringProperty()

    user_nickname = db.StringProperty(indexed=False)

    username = db.StringProperty(default="")

    moderator = db.BooleanProperty(default=False)

    developer = db.BooleanProperty(default=False)

    joined = db.DateTimeProperty(auto_now_add=True)

    @property
    def nickname(self):
        """Gets human friendly name"""
        return self.user_nickname

    @staticmethod
    def current():
        """Return current user (as a users.UserData object).
            Creates a user object if there isn't already one
        """
        
        # Get the Google user object
        user = users.get_current_user()
        uid = user.user_id()
        email = user.email()

        # do we already have an object for this user?
        existing = UserData.get_from_request_info(uid, email)
        if existing:
            return existing
        else:
            return UserData.insert_for(uid, email)

    @staticmethod
    def get_from_request_info(user_id,email=None):
        if not user_id:
            return None

        existing = (UserData.get_from_user_id(user_id) or 
                    UserData.get_from_db_key_email(email))
        if existing:
            return existing

    @staticmethod
    def get_from_user_id(user_id):
        if not user_id:
            return None

        query = UserData.all()
        query.filter('user_id =' ,user_id)
        return query.get()

    @staticmethod
    def get_from_db_key_email(email):
        if not email:
            return None

        query = UserData.all()
        query.filter('user =', users.User(email))
        return query.get()

    @staticmethod
    def insert_for(user_id, email,
                   username=None):
        """Create a user with the specified values, if possible, or returns
        an existing user if the user_id has been used by an existing user.

        Arguments:
            user_id:
                The unique user_id of the user. Should be non-empty.
            email:
                The unique email address of the user. This corresponds
                to the "key_email" value of UserData, and will also be used
                for the "email" property, if non-empty. This can be empty
                for child accounts
            username:
                The unique username to set for this user, if any. Optional.
        
        Returns:
            None if user_id or email values are invalid, or the user couldn't
            be created for other reasons. Otherwise, returns the created,
            or existing user if this user_id is being re-used.
        """

        if not user_id:
            # Every account requires a user_id
            return None
        
        if not email:
            # Every account that isn't a child account requires an e-mail.
            return None
        
        # Make default dummy values for the ones that don't matter
        prop_values = {
            'moderator': False,
            'developer': False
        }

        # Forcefully override with important items.
        db_user_email = email or user_id
        user_email = email

        user = users.User(db_user_email)
        key_name = UserData.key_for(user_id)
        for pname, pvalue in {
            'key_name': key_name,
            'user': user,
            'current_user': user,
            'user_id': user_id,
            'user_email': user_email,
            }.iteritems():
            if pname in prop_values:
                logging.warning("UserData creation about to override"
                                " specified [%s] value" % pname)
            prop_values[pname] = pvalue

        if username:
            # Username or passwords are separate entities.
            # That means we have to do this in multiple steps - make a txn.
            def create_txn():
                user_data = UserData.get_by_key_name(key_name)
                if user_data is not None:
                    logging.warning("Tried to re-make a user for key=[%s]" %
                                    key_name)
                    return user_data

                user_data = UserData(**prop_values)
                # Both claim_username and set_password updates user_data
                # and will call put() for us.
                
                return user_data

            user_data = transaction_util.ensure_in_transaction(create_txn,
                                                               xg_on=True)

        else:
            # No username means we don't have to do manual transactions.
            # Note that get_or_insert is a transaction itself, and it can't
            # be nested in the above transaction.
            user_data = UserData.get_or_insert(**prop_values)
            user_data = db.get(user_data.key())   # force-commit for HRD data

        if user_data and not user_data.is_phantom:
            # Record that we now have one more registered user
            if (datetime.datetime.now() - user_data.joined).seconds < 60:
                # Extra safety check against user_data.joined in case some
                # subtle bug results in lots of calls to insert_for for
                # UserData objects with existing key_names.
                user_counter.add(1)

        return user_data
