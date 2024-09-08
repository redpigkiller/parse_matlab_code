import pkg_resources
from ruamel.yaml import YAML
from typing import *
from enum import Enum, auto
import re

from .parser_error import *

class TokenType(Enum):
    NUMBER = auto()
    IMAG_NUMBER = auto()
    IDENTIFIER = auto()
    KEYWORD = auto()
    OPERATOR = auto()
    PARENTHESIS_OPEN = auto()
    PARENTHESIS_CLOSE = auto()
    BRACKET_OPEN = auto()
    BRACKET_CLOSE = auto()
    BRACES_OPEN = auto()
    BRACES_CLOSE = auto()
    SEMICOLON = auto()
    COLON = auto()
    COMMA = auto()
    CHAR_ARRAY = auto()
    STRING = auto()
    COMMENT = auto()
    NEWLINE = auto()
    AT = auto()
    ELLIPSIS = auto()
    EXCLAMATION = auto()
    QUESTION = auto()
    ASSIGN = auto()
    PERIOD = auto()
    SPACE = auto()

class Token:
    def __init__(self, typ: TokenType, value: str, line: int, column: int) -> None:
        self.typ = typ
        self.val = value
        self.ln = line
        self.col = column

    def __repr__(self) -> str:
        if self.val == "\n":
            return f"Token({self.typ}, '\\n', Ln {self.ln}, Col {self.col})"
        else:
            return f"Token({self.typ}, '{self.val}', Ln {self.ln}, Col {self.col})"

class Tokenizer:
    def __init__(self) -> None:
        yaml = YAML(typ='safe', pure=True)

        file_path = pkg_resources.resource_filename('parse_matlab_code', 'grammar/matlab_keywords.yaml')
        with open(file_path, 'r') as f:
            self.keywords = yaml.load(f)

        file_path = pkg_resources.resource_filename('parse_matlab_code', 'grammar/matlab_operators.yaml')
        with open(file_path, 'r') as f:
            self.operators = yaml.load(f)

        self._match_identifier  = re.compile(r'[a-zA-Z_]\w*')
        self._match_number      = re.compile(r'(?:\d*\.\d+|\d+\.?)[eE][+-]?\d+|0[xX][0-9a-fA-F]+|\d*\.\d+|\d+\.?')

    def tokenize(self, code: str) -> list[Token]:
        # For single quotation ': might be string or hermitian
        self._flag_is_transpose_quotation = False

        # For multi-line comment
        self._multi_line_comment = ""
        self._multi_line_comment_ln_col = (0, 0)
        self._flag_in_multi_line_comment = False

        tokens = []
        lines = code.split('\n')

        # Line and column number
        self._line_num = 1
        self._column_num = 1

        try:
            while self._line_num <= len(lines):
                line = lines[self._line_num-1]

                self._column_num = 1
                while self._column_num <= len(line):
                    token, length = self._get_next_token(line[self._column_num-1:])
                    if token:
                        if token.typ == TokenType.IDENTIFIER or \
                            token.typ == TokenType.PARENTHESIS_CLOSE or \
                            token.typ == TokenType.BRACKET_CLOSE or \
                            token.typ == TokenType.BRACES_CLOSE:
                            self._flag_is_transpose_quotation = True
                        else:
                            self._flag_is_transpose_quotation = False

                        tokens.append(token)
                    
                    self._column_num += length or 1  # Move at least one character forward
                    
                # Add newline if not in multi-line comment block
                if not self._flag_in_multi_line_comment:
                    tokens.append(Token(TokenType.NEWLINE, '\n', self._line_num, len(line) + 1))
                
                # If there is ellipsis at the end, then do not reset transpose quotation flag
                if not tokens[-1] == TokenType.ELLIPSIS:
                    self._flag_is_transpose_quotation = False

                self._line_num += 1
                
        except Exception as e:
            raise TokenizationError(str(e), self._line_num, self._column_num)
        
        if self._flag_in_multi_line_comment:
            raise TokenizationError("Invalid multi-line comments", self._line_num, self._column_num)

        return tokens

    def _get_next_token(self, string: str) -> tuple[Optional[Token], int]:
        
        if self._flag_in_multi_line_comment:
            if string.lstrip().startswith('%}'):
                num_spaces = len(string) - len(string.lstrip())
                self._multi_line_comment += string[:num_spaces + 2]
                # Note: MATLAB seems that it does not accept any string behind %}
                # ret_token = Token(TokenType.COMMENT, self._multi_line_comment, line, column), 2
                ret_token = Token(TokenType.COMMENT, self._multi_line_comment, self._multi_line_comment_ln_col[0], self._multi_line_comment_ln_col[1]), num_spaces + 2
                
                self._multi_line_comment = ""
                self._multi_line_comment_ln_col = (0, 0)
                self._flag_in_multi_line_comment = False
                return ret_token
            
            else:
                self._multi_line_comment += string + "\n"
                return None, len(string)

        # Get the space of the string at front
        if string[0].isspace():
            # Count the number of spaces
            num_spaces = len(string) - len(string.lstrip())
            return Token(TokenType.SPACE, string[:num_spaces], self._line_num, self._column_num), num_spaces

        # Check for comments        
        if self._column_num == 1 and string.startswith('%{'):
            self._multi_line_comment = "%{\n"
            self._multi_line_comment_ln_col = (self._line_num, self._column_num)
            self._flag_in_multi_line_comment = True
            # Note: MATLAB seems that it does not accept any string behind %{, hence, return len(string)
            return None, len(string)
        
        if string.startswith('%'):
            return Token(TokenType.COMMENT, string, self._line_num, self._column_num), len(string)

        # Check for strings
        string_token = self._get_string_token(string, self._line_num, self._column_num)
        if string_token:
            return string_token

        # Check for identifiers and keywords
        match = self._match_identifier.match(string)
        if match:
            value = match.group()
            typ = TokenType.KEYWORD if value in self.keywords else TokenType.IDENTIFIER
            return Token(typ, value, self._line_num, self._column_num), match.end()

        # Check for operators (including '=')
        for op in sorted(self.operators, key=len, reverse=True):
            if string.startswith(op):
                return Token(TokenType.OPERATOR, op, self._line_num, self._column_num), len(op)

        # Check for numbers (including integers, floating point, scientific notation, and hexadecimal)
        match = self._match_number.match(string)
        if match:
            len_of_match = match.end()

            # Check for numbers (imaginary numbers)
            if len_of_match < len(string) and string[len_of_match] in ("i", "j"):
                return Token(TokenType.IMAG_NUMBER, match.group() + string[len_of_match], self._line_num, self._column_num), len_of_match + 1
            else:
                return Token(TokenType.NUMBER, match.group(), self._line_num, self._column_num), len_of_match
        
        # Check for symbols
        if string[0] == '=':
            return Token(TokenType.ASSIGN, '=', self._line_num, self._column_num), 1
        elif string[0] == '(':
            return Token(TokenType.PARENTHESIS_OPEN, '(', self._line_num, self._column_num), 1
        elif string[0] == ')':
            return Token(TokenType.PARENTHESIS_CLOSE, ')', self._line_num, self._column_num), 1
        elif string[0] == '[':
            return Token(TokenType.BRACKET_OPEN, '[', self._line_num, self._column_num), 1
        elif string[0] == ']':
            return Token(TokenType.BRACKET_CLOSE, ']', self._line_num, self._column_num), 1
        elif string[0] == '{':
            return Token(TokenType.BRACES_OPEN, '{', self._line_num, self._column_num), 1
        elif string[0] == '}':
            return Token(TokenType.BRACES_CLOSE, '}', self._line_num, self._column_num), 1
        elif string[0] == ';':
            return Token(TokenType.SEMICOLON, ';', self._line_num, self._column_num), 1
        elif string[0] == ',':
            return Token(TokenType.COMMA, ',', self._line_num, self._column_num), 1
        elif string[0] == ':':
            return Token(TokenType.COLON, ':', self._line_num, self._column_num), 1
        elif string[0] == '@':
            return Token(TokenType.AT, '@', self._line_num, self._column_num), 1
        elif string.startswith('...'):
            return Token(TokenType.ELLIPSIS, '...', self._line_num, self._column_num), 3
        elif string[0] == '!':
            return Token(TokenType.EXCLAMATION, '!', self._line_num, self._column_num), 1
        elif string[0] == '?':
            return Token(TokenType.QUESTION, '?', self._line_num, self._column_num), 1
        elif string[0] == '.':
            return Token(TokenType.PERIOD, '.', self._line_num, self._column_num), 1
        else:
            # If we get here, we have an unknown character, regard it as string input
            first_word = string.split()[0]
            return Token(TokenType.STRING, first_word, self._line_num, self._column_num), len(first_word)

    def _get_string_token(self, string: str, line: int, column: int) -> Optional[tuple[Token, int]]:
        if string[0] not in ("'", '"'):
            return None
        
        if self._flag_is_transpose_quotation:
            return None

        quote = string[0]
        i = 1
        content = []
        while i < len(string):
            if string[i] == quote:
                if i + 1 < len(string) and string[i + 1] == quote:
                    # This is an escaped quote
                    content.append(quote)
                    i += 2
                else:
                    # This is the end of the string
                    break
            else:
                content.append(string[i])
                i += 1

        if i == len(string):
            return None

        full_string = string[:i+1]

        if full_string[0] == "'":
            full_string = f"\"{full_string[1:-1]}\""
            return Token(TokenType.CHAR_ARRAY, full_string, line, column), i + 1
        else:
            return Token(TokenType.STRING, full_string, line, column), i + 1
