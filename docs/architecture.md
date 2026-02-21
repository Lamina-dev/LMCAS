# 核心架构

LMCAS 采用表达式树（AST）+访问者模式，支持多类型符号节点。

- SymbolicExpr：符号表达式包装类
- SymbolicNode：表达式树节点基类
  - NumberNode、VariableNode、AddNode、MultiplyNode、PowerNode、FunctionNode、MatrixNode 等
- Visitor：操作表达式的访问者（如微分、积分、化简、展开等）

底层数值库采用 LAMMP，支持高精度整数、有理数、复数。

详细模块请见 [modules/README.md](modules/README.md)。
