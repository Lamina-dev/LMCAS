#!/usr/bin/env python3
"""
用法：
    python generate_cpp.py expr1 [expr2 ...]
    或
    python generate_cpp.py -f exprs.txt

将数学表达式转换为 LMCAS 兼容的 C++ 代码，自动生成主函数、变量声明和基本操作演示（打印、求导、化简）。
支持特性：
- 基础运算：+, -, *, /, ^
- 常用函数：sin, cos, ln, exp, sqrt, abs 等
- 数据结构：matrix, vector
- 隐式乘法：2x, x(x+1)
"""

import sys
import re
from typing import Optional, Iterator, Any, List, Set, Tuple

# ==========================================
# 词法分析 (Lexer)
# ==========================================
TOKEN_REGEX = re.compile(r"\s*(?:(\d+\.\d*|\d*\.\d+|\d+)|([A-Za-z_]\w*)|(.))")

class Token:
    def __init__(self, typ, value):
        self.typ = typ
        self.value = value
    def __repr__(self):
        return f"Token({self.typ!r}, {self.value!r})"

def tokenize(s) -> Iterator[Token]:
    for num, ident, other in TOKEN_REGEX.findall(s):
        if num:
            yield Token('NUMBER', num)
        elif ident:
            yield Token('IDENT', ident)
        elif other:
            yield Token('SYMBOL', other)
    yield Token('EOF', '')

class ParserError(Exception):
    pass

# ==========================================
# 抽象语法树 (AST)
# ==========================================

class Node:
    def to_cpp(self) -> str:
        raise NotImplementedError
    def collect_vars(self, s: Set[str]):
        pass

class Number(Node):
    def __init__(self, value, denom=None):
        self.value = value
        self.denom = denom # 如果存在，这表示一个分数 value/denom
        
    def to_cpp(self):
        if self.denom:
             return f"SymbolicExpr::number(Rational({self.value}, {self.denom}))"
        if '.' in str(self.value):
             return f"SymbolicExpr::number({self.value})"
        return f"SymbolicExpr::number(BigInt({self.value}))"

class Variable(Node):
    def __init__(self, name):
        self.name = name
    def to_cpp(self):
        return self.name
    def collect_vars(self, s: Set[str]):
        s.add(self.name)

class Binary(Node):
    def __init__(self, left, right):
        self.left = left
        self.right = right
    def collect_vars(self, s: Set[str]):
        self.left.collect_vars(s)
        self.right.collect_vars(s)

class Add(Binary):
    def to_cpp(self):
        return f"SymbolicExpr::add({self.left.to_cpp()}, {self.right.to_cpp()})"

class Subtract(Binary):
    def to_cpp(self):
        return (
            f"SymbolicExpr::add({self.left.to_cpp()}, "
            f"SymbolicExpr::multiply(SymbolicExpr::number(-1), {self.right.to_cpp()}))"
        )

class Multiply(Binary):
    def to_cpp(self):
        return f"SymbolicExpr::multiply({self.left.to_cpp()}, {self.right.to_cpp()})"

class Divide(Binary):
    def to_cpp(self):
        if isinstance(self.left, Number) and isinstance(self.right, Number) and not '.' in str(self.left.value) and not '.' in str(self.right.value):
            return f"SymbolicExpr::number(Rational({self.left.value}, {self.right.value}))"
        return (
            f"SymbolicExpr::multiply({self.left.to_cpp()}, "
            f"SymbolicExpr::power({self.right.to_cpp()}, SymbolicExpr::number(-1)))"
        )

class Power(Binary):
    def to_cpp(self):
        return f"SymbolicExpr::power({self.left.to_cpp()}, {self.right.to_cpp()})"

class ListNode(Node):
    def __init__(self, elements):
        self.elements = elements
    def collect_vars(self, s: Set[str]):
        for e in self.elements:
            e.collect_vars(s)
    def to_cpp(self):
        return "{" + ", ".join(e.to_cpp() for e in self.elements) + "}"

class Function(Node):
    def __init__(self, name, args):
        self.name = name
        self.args = args # list of Nodes
    
    def collect_vars(self, s: Set[str]):
        for arg in self.args:
            arg.collect_vars(s)

    def to_cpp(self):
        mapping = {
            'sin': 'sin', 'cos': 'cos', 'tan': 'tan',
            'cot': 'cot', 'sec': 'sec', 'csc': 'csc',
            'asin': 'arcsin', 'acos': 'arccos', 'atan': 'arctan',
            'arcsin': 'arcsin', 'arccos': 'arccos', 'arctan': 'arctan',
            'sinh': 'sinh', 'cosh': 'cosh', 'tanh': 'tanh',
            'ln': 'ln', 'exp': 'exp', 'sqrt': 'sqrt',
            'log': 'log', 'abs': 'abs',
            'matrix': 'matrix', 'vector': 'vector'
        }
        
        if self.name not in mapping:
             # Default to general function call or error?
             # For now, treat as unsupported
             raise ParserError(f"Unsupported function: {self.name}")

        cpp_func = mapping[self.name]

        if cpp_func == 'log':
            if len(self.args) == 1:
                return f"SymbolicExpr::ln({self.args[0].to_cpp()})"
            elif len(self.args) == 2:
                return f"SymbolicExpr::log({self.args[0].to_cpp()}, {self.args[1].to_cpp()})"
            else:
                raise ParserError("log expects 1 or 2 arguments")

        if cpp_func == 'vector':
             rows = ", ".join(f"{{{e.to_cpp()}}}" for e in self.args)
             return f"SymbolicExpr::matrix({{{rows}}})"

        if cpp_func == 'matrix':
            rows_cpp = []
            for arg in self.args:
                if not isinstance(arg, ListNode):
                     raise ParserError("Elements of matrix must be lists (rows)")
                rows_cpp.append("{" + ", ".join(e.to_cpp() for e in arg.elements) + "}")
            return f"SymbolicExpr::matrix({{{', '.join(rows_cpp)}}})"

        if len(self.args) != 1:
             raise ParserError(f"Function {self.name} expects 1 argument")
        
        return f"SymbolicExpr::{cpp_func}({self.args[0].to_cpp()})"

# ==========================================
# 语法分析 (Parser)
# ==========================================

class Parser:
    def __init__(self, tokens: Iterator[Token]):
        self.tokens = iter(tokens)
        self.lookahead: Optional[Token] = None
        self._advance()

    def _advance(self):
        try:
            self.lookahead = next(self.tokens)
        except StopIteration:
            self.lookahead = None

    def _peek(self, typ, value=None):
        if self.lookahead is None:
            return False
        if self.lookahead.typ != typ:
            return False
        if value is not None and self.lookahead.value != value:
            return False
        return True

    def _accept(self, typ, value=None):
        if self._peek(typ, value):
            tok = self.lookahead
            self._advance()
            return tok
        return None

    def _expect(self, typ, value=None):
        tok = self._accept(typ, value)
        if not tok:
            raise ParserError(f"Expected {typ} '{value if value else ''}', got {self.lookahead}")
        return tok

    def parse(self):
        node = self.expr()
        if self.lookahead and self.lookahead.typ != 'EOF':
             raise ParserError(f"Unexpected token at end: {self.lookahead}")
        return node
    
    def expr(self):
        node = self.term()
        while True:
            if self._accept('SYMBOL', '+'):
                node = Add(node, self.term())
            elif self._accept('SYMBOL', '-'):
                node = Subtract(node, self.term())
            else:
                break
        return node

    def term(self):
        node = self.factor()
        while True:
            if self._accept('SYMBOL', '*'):
                node = Multiply(node, self.factor())
            elif self._accept('SYMBOL', '/'):
                node = Divide(node, self.factor())
            else:
                if self._peek('NUMBER') or self._peek('IDENT') or self._peek('SYMBOL', '('):
                    node = Multiply(node, self.factor())
                else:
                    break
        return node
    
    def factor(self):
        node = self.primary()
        if self._accept('SYMBOL', '^'):
            node = Power(node, self.factor())
        return node

    def primary(self):
        if self._accept('SYMBOL', '('):
            node = self.expr()
            self._expect('SYMBOL', ')')
            return node
        
        if self._accept('SYMBOL', '['):
            elements = []
            if not self._peek('SYMBOL', ']'):
                while True:
                    elements.append(self.expr())
                    if not self._accept('SYMBOL', ','):
                         break
            self._expect('SYMBOL', ']')
            return ListNode(elements)

        tok = self._accept('NUMBER')
        if tok:
            return Number(tok.value)
        
        tok = self._accept('IDENT')
        if tok:
            name = tok.value
            if self._accept('SYMBOL', '('):
                args = []
                if not self._peek('SYMBOL', ')'):
                    while True:
                        args.append(self.expr())
                        if not self._accept('SYMBOL', ','):
                            break
                self._expect('SYMBOL', ')')
                return Function(name, args)
            else:
                return Variable(name)
        
        raise ParserError(f"Unexpected token in primary: {self.lookahead}")

# ==========================================
# 代码生成
# ==========================================

def generate_code(exprs: List[str]):
    asts = []
    all_vars = set()
    errors = []

    for idx, text in enumerate(exprs):
        try:
            parser = Parser(tokenize(text))
            ast = parser.parse()
            ast.collect_vars(all_vars)
            asts.append((idx + 1, text, ast))
        except ParserError as e:
            errors.append(f"// Error parsing expr {idx+1}: {e} (Source: {text})")

    lines = [
        '#include <iostream>',
        '#include <memory>',
        '#include <vector>',
        '#include "../include/symbolic.hpp"',
        '',
        'int main() {',
        '    // Setup aliases',
        '    using Expr = std::shared_ptr<SymbolicExpr>;',
        ''
    ]

    if all_vars:
        lines.append('    // Variables')
        for v in sorted(list(all_vars)):
            lines.append(f'    auto {v} = SymbolicExpr::variable("{v}");')
        lines.append('')
    
    if asts:
        lines.append('    // Expressions')
        for idx, text, ast in asts:
             cpp_code = ast.to_cpp()
             var_name = f"expr{idx}"
             
             lines.append(f'    // {idx}. {text}')
             lines.append(f'    auto {var_name} = {cpp_code};')
             lines.append(f'    std::cout << "Expr {idx}: " << {var_name}->to_string() << std::endl;')
             
             if 'x' in all_vars:
                 lines.append(f'    auto diff{idx} = {var_name}->differentiate("x");')
                 lines.append(f'    std::cout << "  Diff(x): " << (diff{idx} ? diff{idx}->to_string() : "null") << std::endl;')
             
             lines.append(f'    auto simp{idx} = {var_name}->simplify();')
             lines.append(f'    std::cout << "  Simp:    " << (simp{idx} ? simp{idx}->to_string() : "null") << std::endl;')
             lines.append(f'    std::cout << "----------------------------------------" << std::endl;')
             lines.append('')

    if errors:
        lines.append('')
        lines.append('    // Parsing Errors encountered:')
        for err in errors:
            lines.append(f'    {err}')

    lines.append('    return 0;')
    lines.append('}')

    return "\n".join(lines)

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return

    exprs = []
    if sys.argv[1] == '-f':
        if len(sys.argv) < 3:
            print("Error: Missing filename after -f", file=sys.stderr)
            sys.exit(1)
        try:
            with open(sys.argv[2], 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        exprs.append(line)
        except Exception as e:
            print(f"Error reading file: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        exprs = sys.argv[1:]

    cpp_source = generate_code(exprs)
    # Output to stdout. User can redirect to file if needed.
    # e.g. python generate_cpp.py "x+1" > output.cpp
    try:
        sys.stdout.reconfigure(encoding='utf-8') # type: ignore
    except:
        pass
    print(cpp_source)

if __name__ == '__main__':
    main()
