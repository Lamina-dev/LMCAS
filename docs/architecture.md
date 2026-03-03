# LMCAS 架构设计文档

本文档详细描述了 Lamina Computer Algebra System (LMCAS) 的技术架构、核心设计决策及模块交互方式。

## 1. 系统概述

LMCAS 是一个基于 C++17 标准开发的轻量级符号计算系统。其核心目标是提供一个类型安全、内存自动管理且易于扩展的代数计算引擎。系统采用分层架构，底层依赖高精度数值库 LAMMP，上层通过抽象语法树（AST）实现符号运算。

### 架构分层

```mermaid
graph TD
    User[用户代码 / 测试] --> API[SymbolicExpr API]
    API --> AST[抽象语法树 (AST)]
    API --> Visitors[访问者 (Visitors)]
    
    AST --> Nodes[SymbolicNode 及其子类]
    Visitors --> Algorithm[核心算法 (Simplify, Solve, Integrate)]
    
    Algorithm --> Poly[多项式系统 (Polynomial Bridge)]
    Algorithm --> Matrix[线性代数模块]
    
    Poly --> LAMMP[LAMMP (BigInt, Rational)]
    Nodes --> LAMMP
```

## 2. 核心设计原则

### 2.1 不可变性 (Immutability)

LMCAS 中的符号表达式是**不可变**的。
- 所有的代数操作（如 `add`, `multiply`, `sin`）都不会修改现有对象，而是返回一个新的 `SymbolicExpr` 实例。
- **优势**：
    - 线程安全：无锁即可在多线程环境共享表达式。
    - 结构共享：子表达式（如 `x`）可以在多个父表达式中被复用，显著降低内存占用。
    - 调试友好：表达式状态不会在某个深层函数调用中被意外篡改。

### 2.2 自动内存管理 (RAII)

系统广泛使用 `std::shared_ptr<SymbolicNode>` 管理节点生命周期。
- `SymbolicExpr` 本质上是一个轻量级的智能指针包装器。
- 当没有任何表达式引用某个子树时，该子树的内存会自动释放。
- 通过 `std::enable_shared_from_this` 支持安全的内部引用传递。

### 2.3 访问者模式 (Visitor Pattern)

为了解决表达式类型扩展与操作扩展的冲突，LMCAS 采用访问者模式将**数据结构**与**算法**解耦。
- **SymbolicNode** 定义了数据结构（数值、变量、加法、乘法等）。
- **SymbolicVisitor** 定义了操作接口（求值、打印、微分、化简等）。
- 新增算法只需实现一个新的 Visitor，无需修改现有的 Node 类。

## 3. 符号表达层 (Symbolic Layer)

### 3.1 类层次结构

核心类位于 `include/symbolic_ast.hpp` 和 `include/symbolic.hpp`。

*   **SymbolicExpr**: 用户面临的主要接口类。不包含虚函数，持有一个指向 `SymbolicNode` 的 `shared_ptr`。
*   **SymbolicNode**: 所有 AST 节点的抽象基类。
    *   **原子节点 (Leaf Nodes)**:
        *   `NumberNode`: 存储数值（`double`, `BigInt`, `Rational`）。
        *   `VariableNode`: 存储变量名（字符串）。
    *   **复合节点 (Composite Nodes)**:
        *   `AddNode`: N 元加法。
        *   `MultiplyNode`: N 元乘法。
        *   `PowerNode`: 幂运算 (`base ^ exponent`)。
        *   `FunctionNode`: 标准数学函数 (`sin`, `cos`, `exp`, `log` 等)。
        *   `MatrixNode`: 符号矩阵。
        *   `RelationalNode`: 关系表达式 (`=`, `<`, `>`)，用于方程和不等式。

### 3.2 规范化 (Canonical Form)

为了简化比较和化简逻辑，LMCAS 在构造表达式时会自动应用规范化规则（主要通过 `NormalizationVisitor`）：
- **扁平化**: `(a + b) + c` 会被自动转化为 `a + b + c`。
- **排序**: 加法和乘法操作数会根据哈希值与度数（Degree）进行确定性排序。
    - 例如：`y + x + 2` 总是被存储为 `2 + x + y`。
- **合并**: `x + x` 自动转化为 `2*x`；`x * x` 自动转化为 `x^2`。
- **常数折叠**: `2 + 3` 自动计算为 `5`。

这确保了数学上等价的表达式（如 `a+b` 和 `b+a`）在内存中拥有相同的结构，可以直接通过结构比较判断相等性。

## 4. 核心算法实现机制

### 4.1 多项式桥接 (Polynomial Bridge)

为了利用高效的多项式算法（如 GCD、因式分解、求根），LMCAS 实现了一个独特的抽象层：`Polynomial<T>`。

通常多项式系数是数字，但在 LMCAS 中，我们定义了 **`SymbolicPolyCoeff`**。
这意味着我们可以将任意表达式视为关于变量 $x$ 的多项式，其系数本身可以是包含其他变量（$y, z$）的复杂符号表达式。

**工作流**：
1.  用户求解方程 `ax^2 + bx + c = 0`。
2.  `symbolic_to_poly(expr, "x")` 将表达式树转换为 `Polynomial<SymbolicPolyCoeff>`。
3.  系统识别出系数 $a, b, c$。
4.  调用通用多项式求解器（利用求根公式或数值方法）。
5.  结果被重新包装为 `SymbolicExpr` 返回。

### 4.2 模式匹配 (Pattern Matcher)

`Matcher` 模块允许基于规则的重写。它支持通配符匹配。
- 规则：`sin(x)^2 + cos(x)^2 -> 1`
- 匹配器会遍历 AST，寻找结构匹配的子树，并将通配符绑定到实际子表达式上，然后应用变换。

### 4.3 方程求解 (Solver)

求解器采用混合策略：
1.  **隔离变量**: 对于简单的 $f(x) = c$，尝试逆运算（如 $e^x = 5 \implies x = \ln 5$）。
2.  **多项式转换**: 将方程视为多项式，求解根。
3.  **线性系统**: 对于线性方程组，提取系数矩阵并进行高斯消元。


## 5. 构建与扩展

项目使用 CMake 构建系统。
- 核心库编译为 `lmcas.dll` / `liblmcas.so`。
- 添加新函数的步骤：
    1. 在 `SymbolicNode` 中定义新类型（或复用 `FunctionNode`）。
    2. 在 `Extension` 模块注册导数规则。
    3. 在 `NormalizationVisitor` 添加化简规则。

