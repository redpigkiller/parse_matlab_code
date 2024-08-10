from rich import print as rprint
from rich.panel import Panel

class ParserError(Exception):
    """Base exception for MATLAB parser errors."""
    pass

class TokenizationError(ParserError):
    """Exception for errors during tokenization."""
    def __init__(self, message, line, column):
        self.message = message
        self.line = line
        self.column = column
        super().__init__(f"Error at line {line}, column {column}: {message}")

class MatlabSyntaxError(ParserError):
    """Exception for syntax errors in MATLAB code."""
    def __init__(self, message, line, column):
        self.message = message
        self.line = line
        self.column = column
        super().__init__(f"Syntax error at line {line}, column {column}: {message}")

class NotImplementedFeatureError(ParserError):
    """Exception for non-implemented"""
    def __init__(self, keyword, line, column):
        self.keyword = keyword
        self.line = line
        self.column = column
        super().__init__(f"Keyword {keyword} at line {line}, column {column} does not implemented yet")

class LarkParserError(ParserError):
    """Exception for Lark parser"""
    def __init__(self, matlab_expression, lark_error_message, line, column):
        self.matlab_expression = matlab_expression
        self.lark_error_message = lark_error_message
        self.line = line
        self.column = column
        super().__init__(f"Lark parser has error parsing code at line {line} col {column}: {matlab_expression}\n{lark_error_message}")

class TreeConversionError(ParserError):
    """Exception for conversion from lark tree to custom tree"""
    def __init__(self, message, line, column):
        self.message = message
        self.line = line
        self.column = column
        super().__init__(f"Conversion failed at line {line}, column {column}: {message}")
        
class IncompleteError(ParserError):
    """Exception for incomplete code"""
    def __init__(self):
        super().__init__(f"The code is incomplete.")

class PaserImplementationError(ParserError):
    """Exception for implementation of the parser"""
    def __init__(self):
        super().__init__(f"This error is occurred due to the bad implementation of the parser.")

def print_tokenizer_error(e: TokenizationError):
    rprint(f"[bright_red bold]Tokenization error at line {e.line}, column {e.column}: {e.message}[/]")
    
def print_syntax_error(e: MatlabSyntaxError):
    rprint(f"[bright_red bold]Syntax error at line {e.line}, column {e.column}:[/]")
    rprint(Panel.fit(e.message))

def print_lark_parser_error(e: LarkParserError):
    rprint(f"[bright_red bold]Lark parser has error parsing code at line {e.line} col {e.column}:[/]")
    rprint(Panel.fit(e.matlab_expression))
    rprint(f"[bright_red bold]With the following error:[/]")
    rprint(Panel.fit(e.lark_error_message))

def print_tree_conversion_error(e: TreeConversionError):
    rprint(f"[bright_red bold]Conversion failed at line at line {e.line}, column {e.column}:[/]")
    rprint(Panel.fit(e.message))
    
def print_not_implemented_error(e: NotImplementedFeatureError):
    rprint(f"[red bold]Keyword {e.keyword} at line {e.line}, column {e.column} does not implemented yet[/]")