from jinja2 import Environment, FileSystemLoader
import os

class jTemplate():
    _e = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(__file__), 'templates')))

    @staticmethod
    def render(template, variables, printer):
        t = jTemplate._e.get_template(template)
        printer(t.render(variables))
