import pkg_resources
from pathlib import Path
from typing import *
from rich import print as rprint
from lark import Lark

from .tokenizer import Token, TokenType, Tokenizer
from .tree import *
from .convert_lark_tree import convert_lark_tree

from .parser_error import *


class ASTBuilder:
    _ignored_token_type = (TokenType.SPACE, )

    def __init__(self) -> None:
        # Load grammar
        file_path = pkg_resources.resource_filename('parse_matlab_code', 'grammar/matlab_ebnf.lark')
        with open(file_path, 'r') as f:
            lark_grammar = f.read()

        self._lark_parser = Lark(lark_grammar, parser="lalr")
        self._lark_parser_earley = Lark(lark_grammar, parser="earley")

    def parse(self, tokens: list[Token]) -> tuple[Node, list[Token], list[Token]]:
        self._tokens = tokens
        self._current = 0

        # Remove line continuation and comment at the begining
        comment_tokens, line_continuation_tokens = self._preprocessing()

        return self._build_ast(), comment_tokens, line_continuation_tokens
    
    def _build_ast(self) -> Node:
        """
            Build the parse tree.
            [Return] AST, comment_nodes, line_continuation_nodes
        """
        self._flag_skip_space = True    # Default: skip space

        body = []
        while not self._is_at_end():
            body.append(self._parse_statement())

        # Filter out None elements
        body = [i for i in body if i is not None]

        return Program("Root", body)
    
    def _preprocessing(self) -> tuple[list[Token], list[Token]]:
        """
        Deal with
        1. line continuation: ellipsis ("...")
        2. comment: percent ("%")
        """
        # 1. Preprocess "..."
        ellipsis_tokens: list[Token] = []
        idx = 0
        while idx < len(self._tokens):
            if self._tokens[idx].typ == TokenType.ELLIPSIS:
                ellipsis_token = self._tokens.pop(idx)

                # Capture the text behind ellipsis until encounter a newline
                if idx < len(self._tokens) and self._tokens[idx].typ != TokenType.NEWLINE:
                    tmp_token = self._tokens.pop(idx)
                    str_behind_ellipsis = tmp_token.val

                    while idx < len(self._tokens) and self._tokens[idx].typ != TokenType.NEWLINE:
                        str_behind_ellipsis += self._tokens.pop(idx).val
                else:
                    str_behind_ellipsis = ""
                
                # Check for newline behind ellipsis
                assert idx < len(self._tokens), self._error("Expect a newline!")
                self._tokens.pop(idx)   # Pop newline

                # Append the text
                ellipsis_token.val += str_behind_ellipsis
                ellipsis_tokens.append(ellipsis_token)

            else:
                idx += 1

        # 2. Preprocess "%"
        comment_tokens: list[Token] = []
        idx = 0
        while idx < len(self._tokens):
            if self._tokens[idx].typ == TokenType.COMMENT:
                comment_tokens.append(self._tokens.pop(idx))
            else:
                idx += 1
        
        return comment_tokens, ellipsis_tokens
    
    def _parse_statement(self) -> Optional[Node]:
        if self._match(TokenType.SPACE):
            # Ignore whitespace
            return None
        
        elif self._match(TokenType.NEWLINE, TokenType.SEMICOLON, TokenType.COMMA):
            # Match ";", "," and "\n"
            return EndOfLine(self._previous())
        
        elif self._match(TokenType.KEYWORD):
            # Match keyword: if-else, for-loop, ...
            return self._parse_keyword_statement()
        
        elif self._match(TokenType.IDENTIFIER):
            # Match identifier: might be 
            #   1. assignment such as a = 3
            #   2. function call or array indexing such as a(3)
            #   3. expression such as a + 3
            return self._parse_identifier()
        
        elif self._match(TokenType.BRACKET_OPEN):
            # Match "[": might be 
            #   1. return of a function call such as [a, b] = max(X)
            #   2. a simple array [1 ; 2]
            return self._parse_expression_or_function_return()
        
        else:
            # Match number: must be expression such as 1 + a * 6
            # Match "(": must be expression such as (1 + 2) * 3
            # Match "{": must be cell such as {5, 4}
            return self._parse_expression()

    def _parse_keyword_statement(self) -> Node:
        keyword = self._previous().val

        if keyword == "function":
            return self._parse_function_definition()
        elif keyword == "if":
            return self._parse_if_statement()
        elif keyword == "for":
            return self._parse_for_statement()
        elif keyword == "parfor":
            return self._parse_parfor_statement()
        elif keyword in "spmd":
            return self._parse_spmd_statement()
        elif keyword == "while":
            return self._parse_while_statement()
        elif keyword == "switch":
            return self._parse_switch_statement()
        elif keyword == "try":
            return self._parse_try_statement()
        elif keyword == "break":
            return BreakStatement(Token(TokenType.KEYWORD, "break", self._previous().ln, self._previous().col))
        elif keyword == "continue":
            return ContinueStatement(Token(TokenType.KEYWORD, "continue", self._previous().ln, self._previous().col))
        elif keyword == "return":
            return ReturnStatement(Token(TokenType.KEYWORD, "return", self._previous().ln, self._previous().col))
        elif keyword == "global":
            return self._parser_global_statement()
        elif keyword == "persistent":
            return self._parser_persistent_statement()
        elif keyword == "classdef":
            raise NotImplementedFeatureError(keyword, self._previous().ln, self._previous().col)
        else:
            self._error(f"Unexpected keyword: {keyword}")

    def _parse_identifier(self) -> Node:
        """
        1. Function call: cd test/this_is_test_folder
        2. Function call cd(2), expression such as a + 3, or assignment a = 3 + 8
        """
        # Record current position
        current_backup = self._current

        # Check case 1
        #   Function call: cd test/this_is_test_folder
        # Note: Give that cd = 2 and examples = 3, we have the following results
        #   1. cd./examples         ==> 2./3
        #   2. cd./ examples        ==> 2./3
        #   3. cd ./examples        ==> change directory to examples
        #   4. cd ./ examples       ==> 2./3
        self._flag_skip_space = False
        if self._match(TokenType.SPACE):
            # It is case 3 and 4
            is_function_call = False

            if self._match(TokenType.IDENTIFIER) or \
                self._match(TokenType.CHAR_ARRAY) or \
                self._match(TokenType.STRING) or \
                self._match(TokenType.PERIOD) or \
                self._match(TokenType.BRACKET_OPEN) or \
                self._match(TokenType.BRACKET_CLOSE) or \
                self._match(TokenType.BRACES_OPEN) or \
                self._match(TokenType.BRACES_CLOSE) or \
                self._match(TokenType.AT) or \
                self._match(TokenType.EXCLAMATION) or \
                self._match(TokenType.QUESTION):
                # It is function call
                is_function_call = True

            elif self._match(TokenType.OPERATOR):
                # Find operator, need to determine if it is an expression
                if not self._check(TokenType.SPACE):
                    # No space found, it is case 3
                    is_function_call = True

            if is_function_call:
                self._current = current_backup - 1

                func_name = Identifier(self._consume(TokenType.IDENTIFIER))
                self._consume(TokenType.SPACE)

                # Obtain input arguments
                input_params: list[Node] = []
                token_list: list[Token] = []
                
                self._is_still_in_expr(reset=True)
                while self._is_still_in_expr():
                    if self._match(TokenType.SPACE):
                        input_params.append(Token(TokenType.STRING, ''.join([t.val for t in token_list]), token_list[0].ln, token_list[0].col))
                        token_list = []
                    else:
                        token_list.append(self._advance())
                
                if token_list:
                    input_params.append(Token(TokenType.STRING, ''.join([t.val for t in token_list]), token_list[0].ln, token_list[0].col))
                
                # Restore
                self._flag_skip_space = True
                return FunctionCall(func_name, input_params)

        # Check case 2
        #   Function call cd(2), expression such as a + 3, or assignment a = 3 + 8
        self._current = current_backup - 1
        left_expr: list[Token] = []
        is_find_assign = False

        self._is_still_in_expr(reset=True)
        while self._is_still_in_expr():
            if self._match(TokenType.ASSIGN):
                is_find_assign = True
                break
            left_expr.append(self._advance())

        # Restore
        self._flag_skip_space = True

        if is_find_assign:
            lvalue = self._parse_expression(left_expr)
            rvalue = self._parse_expression()
            return Assignment(lvalue, rvalue)
        else:
            # For expression, there is no "="
            return self._parse_expression(left_expr)
    
    def _parse_expression_or_function_return(self) -> Node:
        # Firstly, we assume it is a function return, it only accepts ","; does not accept ";" and "\n"
        is_expr = False
        output_params: list[Node] = [] # Multiple output arguments

        # Record current position
        self._current_backup = self._current

        # Don't ignore whitespace
        self._flag_skip_space = False

        # For capturing expression
        while not self._match(TokenType.BRACKET_CLOSE):
            # Obtain output parameters
            curr_param: list[Token] = []

            self._is_still_in_expr(reset=True, end_line_condition=(
                TokenType.COMMA,
                TokenType.SEMICOLON,
                TokenType.NEWLINE,
                TokenType.SPACE,
                TokenType.BRACKET_CLOSE
                ))
            while self._is_still_in_expr():
                curr_param.append(self._advance())

            if self._match(TokenType.SEMICOLON, TokenType.NEWLINE):
                # If there is ; or \n, then it is an expression not a function return
                is_expr = True
                break
            
            if curr_param:
                if len(curr_param) == 1 and (curr_param[0].typ == TokenType.OPERATOR and curr_param[0].val == "~"):
                    output_params.append(Ignore(curr_param[0]))
                else:
                    output_params.append(self._parse_expression(curr_param))
                curr_param = []

            self._match(TokenType.COMMA, TokenType.SPACE)

        # Restore
        self._flag_skip_space = True

        if is_expr or not self._check(TokenType.ASSIGN):
            # If is_expr or not detect something like [a, b] = 
            self._current = self._current_backup - 1    # Go back to the first identifier
            return self._parse_expression()
        else:
            # It is function return
            self._consume(TokenType.ASSIGN, error_message="Function's return values can not have multiple rows")
            rvalue = self._parse_expression()
            return Assignment(output_params, rvalue)

    def _parse_function_definition(self) -> Node:
        # Multiple output arguments
        output_params: list[Node] = []

        # Check for the 1st case: [x, y] = func_name
        if self._match(TokenType.BRACKET_OPEN):
            while not self._match(TokenType.BRACKET_CLOSE):
                output_params.append(
                    Identifier(
                        self._consume(TokenType.IDENTIFIER, error_message="Expect a output variable")
                    )
                )
                self._match(TokenType.COMMA)

            self._consume(TokenType.ASSIGN, error_message="Expect a '='")
            func_name = self._consume(TokenType.IDENTIFIER, error_message="Expect a function name")
        
        # Check for
        # the 2nd case: x = func_name, and
        # the 3rd case: func_name
        else:
            # Get function name or single output argument
            temp = self._consume(TokenType.IDENTIFIER, error_message="Expect a function name or a output variable")
            
            if self._match(TokenType.ASSIGN):
                # It is the 2nd case
                output_params.append(Identifier(temp))
                func_name = self._consume(TokenType.IDENTIFIER, error_message="Expect a function name")
            else:
                # It is the 3rd case
                func_name = temp
            
            func_name = Identifier(func_name)   # Turn into Node

        # For input arguments
        input_params: list[Node] = []
        if self._match(TokenType.PARENTHESIS_OPEN):
            # For the case: func(), func(x), func(x, y)
            while not self._match(TokenType.PARENTHESIS_CLOSE):
                if self._check(TokenType.IDENTIFIER):
                    input_params.append(Identifier(self._consume(TokenType.IDENTIFIER, error_message="Expect an input variable")))
                elif self._check(TokenType.OPERATOR, "~"):
                    input_params.append(Ignore(self._match(TokenType.OPERATOR)))
                else:
                    self._error("Function input argument can only be an identifier or '~'")
                self._match(TokenType.COMMA)
        
        # For function body
        body: list[Node] = []
        # Deal with function that without end (for .m file with only one function)
        while not self._is_at_end() and not self._check(TokenType.KEYWORD, "end"):
            body.append(self._parse_statement())
        
        if not self._is_at_end():
            self._consume(TokenType.KEYWORD, "end")

        return FunctionDefinition(func_name, input_params, output_params, body)

    def _parse_if_statement(self) -> Node:
        if_condition = self._parse_expression()
        then_body: list[Node] = []
        elseif_clauses: list[ElseIfClause] = []
        else_body: list[Node] = []
        
        # For if body
        while not (
                self._check(TokenType.KEYWORD, "elseif") or 
                self._check(TokenType.KEYWORD, "else") or
                self._check(TokenType.KEYWORD, "end")
            ):
            then_body.append(self._parse_statement())
        
        # For elseif body
        while self._check(TokenType.KEYWORD, "elseif"):
            self._advance()

            elseif_condition = self._parse_expression()
            elseif_body: list[Node] = []

            while not (
                self._check(TokenType.KEYWORD, "elseif") or 
                self._check(TokenType.KEYWORD, "else") or
                self._check(TokenType.KEYWORD, "end")
                ):
                elseif_body.append(self._parse_statement())
            
            elseif_clauses.append(ElseIfClause(elseif_condition, elseif_body))

        # For else body
        if self._check(TokenType.KEYWORD, "else"):
            self._advance()

            while not self._check(TokenType.KEYWORD, "end"):
                else_body.append(self._parse_statement())
        
        self._consume(TokenType.KEYWORD, "end") # Consume "end"
        
        return IfStatement(if_condition, then_body, elseif_clauses, else_body)

    def _parse_for_statement(self) -> Node:
        # Match extra open parenthesis: for the case:
        # ```
        #   for (( (ii) = 2 : 20))
        # ```
        is_ahead_assignment = True
        loop_index = None
        loop_array = []
        parenthesis_cnt = 0

        self._flag_skip_space = False

        self._is_still_in_expr(reset=True)
        while self._is_still_in_expr():
            if is_ahead_assignment:
                if self._match(TokenType.IDENTIFIER):
                    loop_index = Identifier(self._previous())
                elif self._match(TokenType.ASSIGN):
                    assert loop_index, MatlabSyntaxError("For-loop index not found", self._peek().ln, self._peek().col)
                    is_ahead_assignment = False
                elif self._match(TokenType.PARENTHESIS_OPEN):
                    parenthesis_cnt += 1
                elif self._match(TokenType.PARENTHESIS_CLOSE):
                    parenthesis_cnt -= 1
                else:
                    self._advance()
            else:
                loop_array.append(self._advance())

        self._flag_skip_space = True
        
        loop_array = loop_array[:-parenthesis_cnt] if parenthesis_cnt > 0 else loop_array   # Discard last parenthesis_cnt parenthesis
        loop_array = self._parse_expression(loop_array)

        body = self._parse_loop_body()

        return ForLoop(loop_index, loop_array, body)

    def _parse_parfor_statement(self) -> Node:
        # Match extra open parenthesis: for the case:
        # ```
        #   for (( (ii) = 2 : 20))
        # ```
        is_ahead_assignment = True
        is_found_comma = False
        loop_index = None
        loop_array = []
        parfor_opt = []
        parenthesis_cnt = 0

        self._flag_skip_space = False

        self._is_still_in_expr(reset=True)
        while self._is_still_in_expr():
            if is_ahead_assignment:
                if self._match(TokenType.IDENTIFIER):
                    loop_index = Identifier(self._previous())
                elif self._match(TokenType.ASSIGN):
                    assert loop_index, MatlabSyntaxError("For-loop index not found", self._peek().ln, self._peek().col)
                    is_ahead_assignment = False
                elif self._match(TokenType.PARENTHESIS_OPEN):
                    parenthesis_cnt += 1
                elif self._match(TokenType.PARENTHESIS_CLOSE):
                    parenthesis_cnt -= 1
                else:
                    self._advance()
            elif is_found_comma:
                parfor_opt.append(self._advance())
            else:
                if self._match(TokenType.COMMA):
                    is_found_comma = True
                    assert parenthesis_cnt > 0, MatlabSyntaxError("Expect at least one parenthesis", self._peek().ln, self._peek().col)
                else:
                    loop_array.append(self._advance())

        self._flag_skip_space = True

        loop_array = self._parse_expression(loop_array)
        parfor_opt = parfor_opt[:-parenthesis_cnt] if parenthesis_cnt > 0 else parfor_opt   # Discard last parenthesis_cnt parenthesis
        parfor_opt = self._parse_expression(parfor_opt)

        body = self._parse_loop_body()

        return ParforLoop(loop_index, loop_array, parfor_opt, body)

    def _parse_spmd_statement(self) -> Node:
        body = self._parse_loop_body()

        return SPMDStatement(body)

    def _parse_while_statement(self) -> Node:
        condition = self._parse_expression()
        body = self._parse_loop_body()

        return WhileLoop(condition, body)

    def _parse_loop_body(self) -> list[Node]:
        # Function body
        body: list[Node] = []
        while not self._check(TokenType.KEYWORD, "end"):
            body.append(self._parse_statement())
        self._consume(TokenType.KEYWORD, "end")

        return body
    
    def _parse_switch_statement(self) -> Node:
        def skip_redundancy() -> None:
            # Match all redundents
            while self._match(TokenType.NEWLINE, TokenType.COMMA, TokenType.SEMICOLON):
                pass
        
        condition = self._parse_expression()
        case_body: list[Node] = []
        case_clauses: list[ElseIfClause] = []
        otherwise_body: list[Node] = []

        skip_redundancy()   # Clear all redundancy

        # For case body
        while self._check(TokenType.KEYWORD, "case"):
            self._advance()

            case_cond = self._parse_expression()
            case_body: list[Node] = []

            while not (
                self._check(TokenType.KEYWORD, "case") or 
                self._check(TokenType.KEYWORD, "otherwise") or
                self._check(TokenType.KEYWORD, "end")
                ):
                case_body.append(self._parse_statement())
            case_clauses.append(CaseClause(case_cond, case_body))

            skip_redundancy()   # Clear all redundancy

        skip_redundancy()   # Clear all redundancy

        # For otherwise body
        if self._check(TokenType.KEYWORD, "otherwise"):
            self._advance()

            while not self._check(TokenType.KEYWORD, "end"):
                otherwise_body.append(self._parse_statement())
        
        # For end
        self._consume(TokenType.KEYWORD, "end")

        return SwitchStatement(condition, case_clauses, otherwise_body)

    def _parse_try_statement(self) -> Node:
        try_body: list[Node] = []
        catch_body: list[Node] = []
        
        # For try body
        while not self._check(TokenType.KEYWORD, "catch"):
            try_body.append(self._parse_statement())
        
        # For catch exception
        self._consume(TokenType.KEYWORD, "catch")
        if self._match(TokenType.IDENTIFIER):
            exception = self._previous()
        else:
            exception = None

        # For catch body
        while not self._check(TokenType.KEYWORD, "end"):
            catch_body.append(self._parse_statement())
        self._consume(TokenType.KEYWORD, "end")
        
        return TryCatchStatement(try_body, exception, catch_body)

    def _parser_global_statement(self) -> Node:
        global_var: list[Node] = []
        while self._check(TokenType.IDENTIFIER):
            global_var.append(Identifier(self._consume(TokenType.IDENTIFIER)))
        return GlobalStatement(global_var)

    def _parser_persistent_statement(self) -> Node:
        persistent_var: list[Node] = []
        while self._check(TokenType.IDENTIFIER):
            persistent_var.append(Identifier(self._consume(TokenType.IDENTIFIER)))
        return PersistentStatement(persistent_var)

    def _parse_expression(self, token_list: Optional[list[Token]] = None) -> Optional[Node]:
        # Collecting expression tokens until newline
        if token_list is None:
            # Don't ignore whitespace
            self._flag_skip_space = False

            token_list = []
            self._is_still_in_expr(reset=True)
            while self._is_still_in_expr():
                token_list.append(self._advance())
            
            # Restore
            self._flag_skip_space = True

        elif len(token_list) == 0:
            return None

        # Note: by testing, the MATLAB grammar might have the following behavior when dealing with cell array:
        #   Assume c is a cell array 
        #   1. All the statement without space, such as" c{1}" is valid
        #   2. "c {1}" is invalid, and "c {1} + 3" is also invalid
        #   3. However, "2 + c {2}" is valid.
        #   4. Moreover, "(c {1})", "[c {1}]", "{c {1}}" are valid.
        #   5. For the contruction of array or cell array, "[c {1}]" is different from "[c{1}]"
        #
        #   Hence, for parsing, if we add "," between identifier and ("[" or "{") if this pattern
        #   is in "[]" or "{}", for example, makes "[c {1}]" to "[c,{1}]" 

        tmp_token_list = token_list
        token_list = []
        open_count = 0

        idx = 0
        while idx < len(tmp_token_list) - 2:
            if tmp_token_list[idx].typ in (TokenType.BRACKET_OPEN, TokenType.BRACES_OPEN):
                open_count += 1
                token_list.append(tmp_token_list[idx])

            elif tmp_token_list[idx].typ in (TokenType.BRACKET_CLOSE, TokenType.BRACES_CLOSE):
                open_count -= 1
                token_list.append(tmp_token_list[idx])
            
            elif open_count > 0 and \
                tmp_token_list[idx].typ == TokenType.IDENTIFIER and \
                tmp_token_list[idx+1].typ == TokenType.SPACE and \
                tmp_token_list[idx+2].typ == TokenType.BRACES_OPEN:
                # Match pattern: (identifier, space, "{")
                token_list.append(tmp_token_list[idx])

                # Change space token to comma token
                fake_comma_token = tmp_token_list[idx+1]
                fake_comma_token.typ = TokenType.COMMA
                fake_comma_token.val = ","
                idx += 1    # Step over space token
                token_list.append(fake_comma_token)
            else:
                token_list.append(tmp_token_list[idx])
            idx += 1
        
        # Append last 2 tokens
        while idx < len(tmp_token_list):
            token_list.append(tmp_token_list[idx])
            idx += 1
        
        # Convert list of tokens back to string
        token_str = self._token_list_to_string(token_list)
        
        # Use lark parser to parse matlab expression
        try:
            # Parse the expression using lalr parser (faster)
            parse_tree = self._lark_parser.parse(token_str)
        except Exception as e:
            # If failed, parse the expression using earley parser (slower)
            try:
                parse_tree = self._lark_parser_earley.parse(token_str)
            except Exception as e:
                raise LarkParserError(token_str, str(e), token_list[0].ln, token_list[0].col)

        try:
            # Convert lark tree to custom tree
            parse_tree = convert_lark_tree(parse_tree)
        except Exception as e:
            raise TreeConversionError(str(e), self._peek().ln, self._peek().col)

        return parse_tree
    
    # #################### For _parse functions ####################
    # _check: check current token for a TokenType and a value. Return True if matched, otherwise False
    # _match: match current token for a list of TokenType. If matched, return True and move on to the next token; else, return False
    # _advance: move on to the next token and return current token (before moving on to the next token)
    # _consume: consume a specfied TokenType and value (optional), if failed, raise an error
    # _peek: return current token
    # _previous: return previous token
    # _consume_until: consume until a certain condition is fullfill

    def _check(self, typ: TokenType, value: str = None) -> bool:
        """
        The check method looks at the current token without consuming it. It's used to test
        if the current token matches a specific type and optionally a specific value.

        Args:
            typ: The TokenType to check for.
            value(optional): The specific value to check for.

        Returns:
            True if the current token matches, False otherwise.

        Raises:
            IncompleteError: When the code is at the end.
        """
        self._check_whitespace()

        if self._is_at_end():
            raise IncompleteError()
        
        if value:
            return self._peek().typ == typ and self._peek().val == value
        else:
            return self._peek().typ == typ

    def _match(self, *typs: tuple[TokenType]) -> bool:
        """
        The match method is similar to check, but it consumes the token if there's a match.
        It can check for multiple token types.

        Args:
            *typs: One or more TokenTypes to check for.

        Returns:
            True if the current token matches any of the given types (and consumes it), False otherwise.

        Raises:
            None
        """
        self._check_whitespace()

        if self._is_at_end():
            return False

        if self._peek().typ in typs:
            self._advance()
            return True
        return False
    
    def _advance(self) -> Token:
        """
        The advance method consumes the current token and moves to the next one.

        Args:
            None

        Returns:
            The token that was just consumed.

        Raises:
            IncompleteError: When the code is at the end.
        """
        self._check_whitespace()

        if self._is_at_end():
            raise IncompleteError()
        
        self._current += 1
        return self._previous()

    def _consume(self, typ: TokenType, value: str = None, error_message: str = None) -> Token:
        """
        The consume method is used when we expect a specific token type (and optionally value) to be present.
        It combines checking and advancing.

        Args:
            typ: The TokenType we expect.
            value (optional): The specific value we expect.
            error_message (optional): The error message if the current token doesn't match the expected type/value.

        Returns:
            The consumed token if it matches.

        Raises:
            An error if the current token doesn't match the expected type/value.
        """
        self._check_whitespace()

        if self._check(typ, value):
            return self._advance()
        
        # If the current token doesn't match the expected type/value
        if error_message:
            self._error(error_message)
        else:
            self._error(f"Expected {typ.name}{f' ({value})' if value else ''}, found {self._peek()}")

    def _peek(self) -> Token:
        """
        The peek method return the current token.

        Args:
            None

        Returns:
            The current token.

        Raises:
            PaserImplementationError: When the _parse function wants to peek a out of range token
        """
        self._check_whitespace()

        if self._is_at_end():
            raise PaserImplementationError()

        return self._tokens[self._current]

    def _previous(self) -> Token:
        """
        The previous method return the previous token.

        Args:
            None

        Returns:
            The previous token.

        Raises:
            PaserImplementationError: When the _parse function wants to return a out of range token
        """
        if not (1 <= self._current <= len(self._tokens)):
            raise PaserImplementationError()
        return self._tokens[self._current - 1]
    
    def _error(self, message: str) -> None:
        raise MatlabSyntaxError(message, self._peek().ln, self._peek().col)

    # #################### Helper functions ####################
    def _check_whitespace(self) -> None:
        if self._flag_skip_space:
            while self._current < len(self._tokens):
                curr_token = self._tokens[self._current]
                if curr_token.typ in self._ignored_token_type:
                    self._current += 1
                else:
                    break

    def _is_at_end(self) -> bool:
        return self._current >= len(self._tokens)
    
    def _token_list_to_string(self, token_list: list[Token]) -> str:
        token_string: list[str] = []
        for token in token_list:
            if token.typ == TokenType.CHAR_ARRAY:
                # For char array, convert 'a' to 65 and 'ab' to [65, 66]
                unicode_string = token.val[1:-1]
                char_arr = [str(ord(c)) for c in unicode_string]
                if len(char_arr) == 1:
                    expr_str = char_arr[0]
                else:
                    expr_str = "[" + ", ".join(char_arr) + "]"
                token_string.append(expr_str)

            else:
                # Others
                token_string.append(token.val)

        return ''.join(token_string)

    def _is_still_in_expr(self,
                          reset: bool=False,
                          end_line_condition: tuple[TokenType] = (TokenType.NEWLINE, TokenType.SEMICOLON, TokenType.COMMA,)
        ) -> bool:
        """
        Return non-zero if exit expression
        """
        if reset:
            self._parenthesis_cnt = 0
            self._bracket_cnt = 0
            self._braces_cnt = 0
            self._end_condition = end_line_condition
            return True
        else:
            if self._is_at_end():
                return False
            
            typ = self._peek().typ  # Get current token's type

            if typ in self._end_condition:
                # end_condition is fulfill
                if self._parenthesis_cnt == 0 and self._bracket_cnt == 0 and self._braces_cnt == 0:
                    # a valid expression
                    return False
            
            # Update counters
            if typ == TokenType.PARENTHESIS_OPEN:
                self._parenthesis_cnt += 1
            elif typ == TokenType.PARENTHESIS_CLOSE:
                self._parenthesis_cnt -= 1
            elif typ == TokenType.BRACKET_OPEN:
                self._bracket_cnt += 1
            elif typ == TokenType.BRACKET_CLOSE:
                self._bracket_cnt -= 1
            elif typ == TokenType.BRACES_OPEN:
                self._braces_cnt += 1
            elif typ == TokenType.BRACES_CLOSE:
                self._braces_cnt -= 1
            return True


class Parser:
    def __init__(self) -> None:
        self._tokenizer = Tokenizer()
        self._ast_builder = ASTBuilder()

        self.parse_tree = None
        self.comment_tokens = None
        self.line_continuation_tokens = None

    def parse(self, code: str) -> Node:
        # Tokenize
        tokens = self._tokenizer.tokenize(code)
        
        # Build AST
        parse_tree, comment_tokens, line_continuation_tokens = self._ast_builder.parse(tokens)

        # Store data
        self.parse_tree = parse_tree
        self.comment_tokens = comment_tokens
        self.line_continuation_tokens = line_continuation_tokens

        return parse_tree

