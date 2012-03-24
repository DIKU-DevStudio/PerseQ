import os
modules = [pyfile[:-3]
           for pyfile in os.listdir(os.path.dirname(__file__))
           if pyfile[-3:] == '.py' and pyfile != '__init__.py']

__routes__ = []
for module in modules:
        mod = __import__(module, globals(), locals(), [])
        __routes__ += mod.__routes__
