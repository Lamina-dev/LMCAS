# 微积分模块

主要头文件：
- `include/integration.hpp`
- `include/visitors/differentiation_visitor.hpp`
- `include/symbolic.hpp`

## 概览 (Overview)

微积分模块提供了符号微分（Differentiation）、符号积分（Integration）以及级数展开（Series Expansion）的核心算法。
用户通常不需要直接实例化 `Integrator` 或 `DifferentiationVisitor` 类，而是通过 `SymbolicExpr` 的成员函数来调用这些功能。

---

## 微分 (Differentiation)

### 入口函数
`SymbolicExpr::differentiate(const std::string& var_name)`

### 实现机制
微分操作通过访问者模式 (`DifferentiationVisitor`) 实现。该访问者遍历表达式树，根据节点的类型应用相应的导数法则：
- **基本法则**
    - 常数 (`NumberNode`) -> 0
    - 变量 (`VariableNode`) -> 1 (若是求导变量) 或 0
    - 加法 (`AddNode`) -> 逐项求导 (Sum Rule)
    - 乘法 (`MultiplyNode`) -> 乘积法则 (Product Rule: `(uv)' = u'v + uv'`)
    - 幂 (`PowerNode`) -> 幂法则与链式法则结合 (Power Rule / Chain Rule)
    - 函数 (`FunctionNode`) -> 链式法则 (Chain Rule: `f(g(x))' = f'(g(x)) * g'(x)`)

### 特性与局限
- **支持**：显式函数求导、多层嵌套函数求导。
- **局限**：目前的实现主要针对单变量显函数，对隐函数求导 (`Implicit Differentiation`) 支持仅作为实验性功能存在 (`include/symbolic_implicit_diff.hpp`)。

---

## 积分 (Integration)

### 入口函数
`SymbolicExpr::integrate(const std::string& var_name)`

### 核心算法 (Integrator)
由于通用的 Risch 算法实现极其复杂，LMCAS 目前采用基于启发式规则（Heuristic Rules）和模式匹配的积分策略。

#### 策略流程
1.  **查表法 (Table Lookup)**：匹配基本积分公式（如 `x^n`, `e^x`, `sin(x)` 等）。
2.  **线性性质**：`∫(a*f + b*g) = a*∫f + b*∫g`。
3.  **换元法 (Substitution)**：
    - 自动检测形如 `f(g(x)) * g'(x)` 的结构，令 `u = g(x)` 简化积分。
4.  **分部积分 (Integration by Parts, IBP)**：
    - 针对乘积结构 `u * v`，尝试应用 `∫u dv = uv - ∫v du`。
    - 使用启发式算法选择 `u` 和 `dv` 以降低被积函数的复杂度。
5.  **部分分式分解 (Partial Fraction Decomposition)**：
    - 对有理函数进行分解，将其拆分为简单的分式和多项式之和。

### 积分历史 (IntegralHistory)
为了防止在分部积分或换元过程中出现循环调用（即积分结果变回原函数或变得更复杂），积分器维护了一个 `IntegralHistory` 栈，检测到重复模式时会强制终止当前路径。

---

## 级数展开 (Series Expansion)

### 入口函数
`SymbolicExpr::series(const std::string& var, const std::shared_ptr<SymbolicExpr>& point, int order)`

### 实现
级数展开基于泰勒公式实现。系统自动计算目标函数在展开点 `point` 处的高阶导数（直到 `order` 阶），并构造多项式：
`f(x) ≈ Σ (f^(n)(a) / n!) * (x-a)^n`

### 注意事项
- 级数展开依赖于高阶导数的计算，对于非常复杂的函数，计算高阶导数可能会非常耗时。
