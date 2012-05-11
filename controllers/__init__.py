import os
__all__ = [pyfile[:-3]
           for pyfile in os.listdir(os.path.dirname(__file__))
           if pyfile[-3:] == '.py' and pyfile != '__init__.py']

__routes__ = []
for module in __all__:
        mod = __import__(module, globals(), locals(), [])
        __routes__ += mod.__routes__
