from dataclasses import dataclass, field
from typing import *

from .tokenizer import Token

@dataclass
class Node:
    pass

@dataclass
class Program(Node):
    name: str
    body: list[Node]

##################################################
# Basis
##################################################
@dataclass
class Number(Node):
    value: Token

@dataclass
class ImagNumber(Node):
    value: Token

@dataclass
class String(Node):
    value: Token

@dataclass
class Identifier(Node):
    value: Token

@dataclass
class Ignore(Node):
    value: Token

@dataclass
class WhiteSpace(Node):
    value: Token

@dataclass
class COLON(Node):
    value: Token

@dataclass
class Comment(Node):
    value: Token

@dataclass
class Ellipsis(Node):
    """
    The text behind ellipsis, for example,
    ```
        a = b + ... add two
            2;
    ```
    then, text = "... add two".
    """
    value: Token

@dataclass
class EndOfLine(Node):
    value: Token

##################################################
# Flow
##################################################
@dataclass
class FunctionDefinition(Node):
    name: Node
    input_params: list[Node]
    output_params: list[Node]
    body: list[Node]

@dataclass
class AnonymousFunction(Node):
    parameters: list[Node]
    body: list[Node]

@dataclass
class ElseIfClause(Node):
    condition: Node
    body: list[Node]

@dataclass
class IfStatement(Node):
    condition: Node
    then_body: list[Node]
    elseif_clauses: list[ElseIfClause] = field(default_factory=list)
    else_body: list[Node] = None

@dataclass
class ForLoop(Node):
    identifier: Node
    expression: Node
    body: list[Node]

@dataclass
class ParforLoop(Node):
    identifier: Node
    expression: Node
    option: Node
    body: list[Node]

@dataclass
class SPMDStatement(Node):
    body: list[Node]

@dataclass
class WhileLoop(Node):
    condition: Node
    body: list[Node]

@dataclass
class CaseClause(Node):
    condition: Node
    body: list[Node]

@dataclass
class SwitchStatement(Node):
    expression: Node
    cases: list[CaseClause]
    otherwise: Optional[list[Node]] = None

@dataclass
class TryCatchStatement(Node):
    try_body: list[Node]
    exception: Node
    catch_body: list[Node]

@dataclass
class BreakStatement(Node):
    value: Token

@dataclass
class ContinueStatement(Node):
    value: Token

@dataclass
class ReturnStatement(Node):
    value: Token

##################################################
# Matlab variables
##################################################
@dataclass
class GlobalStatement(Node):
    identifiers: list[Node]

@dataclass
class PersistentStatement(Node):
    identifiers: list[Node]

@dataclass
class Assignment(Node):
    lvalue: list[Node]
    rvalue: Node

##################################################
# For matlab special variables
##################################################
@dataclass
class FunctionHandle(Node):
    value: Token

@dataclass
class FunctionCall(Node):
    # Note that, array and cell indexing is also regarded as function call
    identifier: Node
    arguments: Node

@dataclass
class ArrayAccess(Node):
    identifier: Node
    arguments: Node

@dataclass
class CellArrayAccess(Node):
    identifier: Node
    arguments: Node

@dataclass
class IndexExpression(Node):
    indices: list[Node]

@dataclass
class StructAccess(Node):
    identifier: Node
    arguments: list[Node]

@dataclass
class MatrixExpression(Node):
    elements: list[list[Node]]

@dataclass
class CellArrayExpression(Node):
    elements: list[list[Node]]

@dataclass
class ColonArray(Node):
    start: Node
    stop: Node
    step: Optional[Node] = None

##################################################
# For operations
##################################################
@dataclass
class BinaryOperation(Node):
    left: Node
    operator: Token
    right: Node

@dataclass
class UnaryOperation(Node):
    operator: Token
    operand: Node
