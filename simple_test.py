from rich import print as rprint

from parse_matlab_code import Parser


code = """
c{2}
c{2}
[c{2}]
[c {2}]
"""


parser = Parser()
parse_tree = parser.parse(code)
rprint(parse_tree)