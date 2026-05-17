# 符号表达式模块

主要头文件：include/symbolic.hpp

## 概览
SymbolicExpr 是 LMCAS 中的核心类，用于包装符号表达式树（AST）。它通过 std::shared_ptr 持有 SymbolicNode 根节点，提供了包括内存管理、表达式创建、数学运算等核心功能。

### 核心特性
- **不可变性 (Immutability)**：所有的修改操作（如展开、求导）均返回新的 SymbolicExpr 实例，保证线程安全与引用透明。
- **自动内存管理**：使用 std::shared_ptr 进行引用计数，当表达式不再被引用时自动释放内存。
- **类型安全**：通过 SymbolicExpr::Type 枚举和动态类型检查确保操作的合法性。

## 成员函数 (Member Functions)

以下方法需在 SymbolicExpr 实例上调用：

### xpand()
- **描述**：对表达式进行代数展开，应用乘法分配律并合并同类项。
- **算法**：递归遍历 AST，将 NumberNode 与 VariableNode 视为基础，将 MultiplyNode 分配到 AddNode 上。
- **示例**：(a+b)^2 -> ^2 + 2ab + b^2

### substitute(var, value)
- **描述**：将表达式中的指定符号变量替换为另一个表达式。
- **参数**：
    - ar: 目标变量名 (std::string)
    - alue: 用于替换的表达式 (std::shared_ptr<SymbolicExpr>)
- **算法**：Visitor 模式遍历，匹配 VariableNode 名称。

### differentiate(var)
- **描述**：对指定变量求导。
- **算法**：依据求导法则（链式法则、乘积法则等）递归构建导数表达式。
- **支持**：多项式、三角函数、指数对数等基本初等函数的求导。

### integrate(var)
- **描述**：对指定变量进行符号积分。
- **算法**：调用 integration 模块，使用启发式策略（由于 Risch 算法过于复杂，目前采用模式匹配与分部积分结合的策略）。

### series(var, point, order)
- **描述**：计算泰勒 (Taylor) 或麦克劳林 (Maclaurin) 级数展开。
- **参数**：
    - ar: 展开变量。
    - point: 展开点 (SymbolicExpr)。
    - order: 展开阶数 (int)。

### simplify()
- **描述**：对表达式进行化简。
- **算法**：综合运用常数折叠、恒等式变换（如 sin^2 + cos^2 = 1）、最大公约数提取等策略。

## 静态算法 (Static Methods)

以下算法为 SymbolicExpr 类的静态成员函数，主要用于多项式与矩阵运算：

### 多项式
- **SymbolicExpr::poly_gcd(a, b)**: 计算两个多项式的最大公约数（GCD）。
- **SymbolicExpr::poly_resultant(a, b, var)**: 计算两个多项式的结式（Resultant）。

### 线性代数 (Matrix)
- **SymbolicExpr::determinant(mat)**: 计算矩阵行列式。支持递归展开与 LU 分解。
- **SymbolicExpr::inverse(mat)**: 计算逆矩阵。
- **SymbolicExpr::eigenvalues(mat)**: 计算特征值（通过特征多项式）。
- **SymbolicExpr::solve_system(eqs, vars)**: 求解线性方程组。

## 模式匹配逻辑

### 实现原理
LMCAS 的模式匹配基于 AST 的结构比对，支持：
- **通配符 (Wildcards)**：匹配任意子树。
- **交换律 (Commutativity)**： + b 可匹配  + a。
- **结合律 (Associativity)**：(a + b) + c 可匹配  + b + c。

### 算法细节
1. **类型检查**：优先比对节点类型（Add, Multiply, Function 等）。
2. **递归比对**：对子节点进行递归匹配。
3. **哈希剪枝**：利用表达式哈希值快速判断是否可能匹配（详见 HashData 结构体）。
4. **全排列尝试**：对于满足交换律的节点（如加法、乘法），在匹配时尝试子节点的全排列组合（注：这是性能瓶颈之一）。

## 节点类型 (SymbolicNode)
- **NumberNode**: 存储数值 (BigInt, Rational, double)。
- **VariableNode**: 存储变量名。
- **AddNode / MultiplyNode**: n元运算节点。
- **PowerNode**: 幂运算。
- **FunctionNode**: 通用函数节点 (sin, cos, exp, ln, lambertw, rootof 等)。
- **MatrixNode**: 矩阵容器。

## 已知问题与局限
- **大表达式性能**：深层嵌套的递归化简可能导致栈溢出或性能显著下降。
- **特殊函数积分**：目前的积分器对非初等积分的支持有限。
- **非线性方程组**：solve_system 主要针对线性系统，非线性支持较弱。
