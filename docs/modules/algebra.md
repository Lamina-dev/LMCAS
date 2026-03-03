# 代数与方程模块

主要头文件：
- `include/polynomial.hpp`
- `include/poly_utils.hpp`
- `include/solver.hpp`

## 概览
本模块提供了底层的多项式运算 (`Polynomial<T>`) 和高层的方程求解器 (`Solver`)。

---

## 多项式运算 (Polynomial Operations)

### `Polynomial<T>` 类
这是一个通用的多项式模板类，支持多种系数类型（如 `BigInt`、`Rational`、甚至 `SymbolicExpr`）。

#### 核心算法
- **加减乘除**：按项合并 (`term-wise`)，自动处理系数运算。
- **除法带余 (`div_mod`)**：计算多项式除法的商和余数。
- **伪除法 (`pseudo_div_mod`)**：用于整系数多项式，避免引入分数。
- **最大公约数 (GCD)**：基于欧几里得算法及其变种（如子结式 PRS 算法）。
- **因式分解 (Factorization)**：目前支持简单的无平方因子分解 (`Square-free factorization`) 和部分特定形式的因式分解。

#### 特性
- **稀疏存储**：仅存储非零项，适合高次稀疏多项式。
- **符号支持**：当 `T = SymbolicExpr` 时，可处理含参数的多项式。

---

## 方程求解 (Equation Solving)

### 入口函数
`Solver::solve(eqn, var)` 或 `SymbolicExpr::solve(eqn, var)`

### 算法策略
求解器根据方程的类型选择不同的策略：
1.  **线性方程 (`Linear`)**：直接移项求解 `ax + b = 0 => x = -b/a`。
2.  **二次方程 (`Quadratic`)**：应用求根公式。
3.  **多项式方程 (`Polynomial`)**：
    - 尝试因式分解。
    - 对于高次方程，尝试寻找有理根或使用数值逼近（如果需要）。
4.  **三角方程 / 超越方程**：
    - 尝试通过换元 (`t = sin(x)`) 转化为多项式方程。
    - 使用 Lambert W 函数求解形如 `x * e^x = y` 的方程（实验性）。

### 线性方程组
对于线性方程组，使用高斯消元法 (`Gaussian Elimination`) 或克拉默法则 (`Cramer's Rule`) 进行求解。

## 局限性
- **五次及以上方程**：无通用求根公式，通常无法获得解析解。
- **混合超越方程**：如 `x + sin(x) = 1`，通常只能求数值解。
