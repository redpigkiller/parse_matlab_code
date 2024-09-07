from lark import Tree, Token as LarkToken

from .tokenizer import Token, TokenType
from .tree import *

from .parser_error import TreeConversionError

def convert_lark_tree(tree, ln_ref: int, col_ref: int) -> Node:
    """
    Converts a Lark parse tree to custom Node objects for MATLAB expressions.
    
    Args:
        tree: A Lark Tree object or Token representing the parsed MATLAB expression.
    
    Returns:
        A custom Node object representing the MATLAB expression.
    
    Raises:
        ValueError: If an unhandled node type is encountered.
    """

    # Handle Lark Tokens (terminal symbols)
    if isinstance(tree, LarkToken):
        if tree.type == 'NUM':
            return Number(Token(TokenType.NUMBER, tree.value, tree.line+ln_ref-1, col_ref if tree.line > 1 else tree.column+col_ref-1))
        elif tree.type == 'IMAG_NUM':
            return ImagNumber(Token(TokenType.IMAG_NUMBER, tree.value, tree.line+ln_ref-1, col_ref if tree.line > 1 else tree.column+col_ref-1))
        elif tree.type == 'IDENTIFIER':
            return Identifier(Token(TokenType.IDENTIFIER, tree.value, tree.line+ln_ref-1, col_ref if tree.line > 1 else tree.column+col_ref-1))
        elif tree.type == 'TILDE':
            return Ignore(Token(TokenType.OPERATOR, tree.value, tree.line+ln_ref-1, col_ref if tree.line > 1 else tree.column+col_ref-1))
        elif tree.type == 'STR':
            return String(Token(TokenType.STRING, tree.value, tree.line+ln_ref-1, col_ref if tree.line > 1 else tree.column+col_ref-1))
        elif tree.type == 'COLON':
            return COLON(Token(TokenType.COLON, tree.value, tree.line+ln_ref-1, col_ref if tree.line > 1 else tree.column+col_ref-1))
        # elif tree.type == 'WS':
        #     return WhiteSpace(Token(TokenType.SPACE, tree.value, tree.line+ln_ref-1, col_ref if tree.line > 1 else tree.column+col_ref-1))
        else:
            # If we encounter an unhandled node type, raise an error with details
            raise TreeConversionError(f"Unhandled terminal in lark grammar file: {tree.type}", tree.line+ln_ref-1, col_ref if tree.line > 1 else tree.column+col_ref-1)
    
    # If it's not a Lark Tree, return it as is
    if not isinstance(tree, Tree):
        return tree

    # Handle different node types
    if tree.data == 'start':
        # The start node just wraps the main expression, so we recurse into its child
        return convert_lark_tree(tree.children[0], ln_ref, col_ref)

    elif tree.data == 'number':
        # Numeric literals
        return convert_lark_tree(tree.children[0], ln_ref, col_ref)

    elif tree.data == 'string_array':
        # Character arrays and string literals
        return convert_lark_tree(tree.children[0], ln_ref, col_ref)
    
    elif tree.data == 'simple_identifier':
        # Simple identifiers (variable names)
        return convert_lark_tree(tree.children[0], ln_ref, col_ref)

    elif tree.data in ('addition', 'subtraction', 'multiplication', 'matrix_multiplication', 
                       'right_division', 'left_division', 'matrix_right_division', 'matrix_left_division',
                       'power', 'matrix_power',
                       'element_wise_or', 'element_wise_and', 'short_circuit_or', 'short_circuit_and',
                       'less_than', 'less_than_or_equal_to', 'greater_than', 'greater_than_or_equal_to',
                       'equal_to', 'not_equal_to'):
        # Binary operations
        return BinaryOperation(
            convert_lark_tree(tree.children[0], ln_ref, col_ref),   # Left operand
            tree.data,                                              # Operator
            convert_lark_tree(tree.children[1], ln_ref, col_ref)    # Right operand
        )
    
    elif tree.data == 'colon_operator':
        # Colon operator (used in ranges and slicing)
        start_tree = tree.children[0]
        if start_tree.data == 'colon_operator':
            start = convert_lark_tree(tree.children[0].children[0], ln_ref, col_ref)
            step  = convert_lark_tree(tree.children[0].children[1], ln_ref, col_ref)
        else:
            start = convert_lark_tree(tree.children[0], ln_ref, col_ref)
            step  = []

        return ColonArray(start, convert_lark_tree(tree.children[1], ln_ref, col_ref), step)

    elif tree.data in ('unary_plus', 'unary_minus', 'logical_negation',
                       'transpose', 'hermitian'):
        # Unary operations
        return UnaryOperation(
            tree.data,                                              # Operator
            convert_lark_tree(tree.children[0], ln_ref, col_ref)    # Operand
        )

    elif tree.data == 'parentheses':
        # Parentheses just wrap an expression, so we recurse into the child
        return convert_lark_tree(tree.children[0], ln_ref, col_ref)

    elif tree.data == 'function_handle':
        # Function handles (@function_name)
        return FunctionHandle(tree.children[0].value)
    
    elif tree.data == 'anonymous_function':
        # Anonymous functions
        input_params = tree.children[0]
        if input_params:
            input_params = [convert_lark_tree(input_param, ln_ref, col_ref) for input_param in input_params.children]
        else:
            input_params = []

        return AnonymousFunction(
            input_params,                                               # Parameters
            [convert_lark_tree(tree.children[1], ln_ref, col_ref)]      # Function body
        )

    elif tree.data == 'function_call':
        # Function calls
        return FunctionCall(
            convert_lark_tree(tree.children[0], ln_ref, col_ref),    # Function identifier
            convert_lark_tree(tree.children[1], ln_ref, col_ref)     # Arguments
        )

    elif tree.data == 'cell_element':
        # Cell array indexing
        return CellArrayAccess(
            convert_lark_tree(tree.children[0], ln_ref, col_ref),    # Cell array identifier
            convert_lark_tree(tree.children[1], ln_ref, col_ref)     # Indices
        )
    
    elif tree.data == 'indexing':
        # Array/Cell indexing
        return IndexExpression(
            [convert_lark_tree(index, ln_ref, col_ref) for index in tree.children]   # Indices
        )

    elif tree.data == 'struct_element':
        # Structure field access
        base = convert_lark_tree(tree.children[0], ln_ref, col_ref)
        field_name = [convert_lark_tree(field, ln_ref, col_ref) for field in tree.children[1:]]
        return StructAccess(base, field_name)
    
    elif tree.data == 'field_name':
        # Structure's field name
        return convert_lark_tree(tree.children[0], ln_ref, col_ref)
    
    elif tree.data in 'matlab_array':
        # Array definitions
        if len(tree.children) == 0 or tree.children[0] is None:
            return MatrixExpression([[]])
        else:
            return MatrixExpression(
                [
                    [convert_lark_tree(elem, ln_ref, col_ref) for elem in row.children]
                    for row in tree.children[0].children
                ]
            )

    elif tree.data in 'matlab_cell':
        # Cell array definitions
        if len(tree.children) == 0 or tree.children[0] is None:
            return CellArrayExpression([[]])
        else:
            return CellArrayExpression([
                [convert_lark_tree(elem, ln_ref, col_ref) for elem in cell_row.children]
                for cell_row in tree.children[0].children
            ])
    
    else:
        # If we encounter an unhandled node type, raise an error with details
        raise TreeConversionError(f"Unhandled node type: {tree.data}\nTree: {tree.pretty()}", ln_ref, col_ref)
