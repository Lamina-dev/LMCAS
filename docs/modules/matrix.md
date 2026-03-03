# 矩阵与向量模块

主要头文件：
- `include/symbolic.hpp`
- `include/symbolic_matrix.hpp`

## 概览 (Overview)
LMCAS 支持符号矩阵（Symbolic Matrix）运算。矩阵可以通过嵌套的 `std::vector` 构造，内部存储为 `MatrixNode`。

---

## 核心线性代数运算 (Linear Algebra)

这些功能主要通过 `SymbolicExpr` 类的静态成员函数提供，是处理一般矩阵问题的首选方式。

### `SymbolicExpr::determinant(mat)`
- **描述**：计算矩阵的行列式。
- **算法**：对于小矩阵使用拉普拉斯展开，对于大矩阵使用高斯消元或 LU 分解。

### `SymbolicExpr::inverse(mat)`
- **描述**：计算逆矩阵。
- **算法**：基于伴随矩阵或高斯-约旦消元法。
- **注意**：如果行列式为零（奇异矩阵），将无法求逆。

### `SymbolicExpr::eigenvalues(mat)`
- **描述**：计算矩阵的特征值。
- **算法**：首先计算特征多项式 `det(A - λI)`，然后求解该多项式的根。
- **局限**：仅适用于特征多项式可解的情况（通常是低阶矩阵或特殊结构矩阵）。

### `SymbolicExpr::poly_resultant(a, b, var)`
- **描述**：虽然属于多项式操作，但在处理矩阵特征值等问题时常用于消元。

---

## 几何变换 (Geometric Transformations)

这些辅助函数位于 `lamina` 命名空间下（定义在 `include/symbolic_matrix.hpp`），用于快速生成几何变换矩阵。

### `lamina::matrix_rotation(theta, dim)`
- **描述**：生成二维或三维旋转矩阵。
- **参数**：
    - `theta`: 旋转角度（弧度）。
    - `dim`: 维度（默认为 2）。

### `lamina::matrix_reflection(angle, dim)`
- **描述**：生成反射矩阵。
- **参数**：
    - `angle`: 反射轴/面的角度。

### `lamina::matrix_scaling(sx, sy, dim)`
- **描述**：生成缩放矩阵。
- **参数**：
    - `sx`, `sy`: X, Y轴的缩放因子。

---

## 向量运算 (Vector Operations)

向量在 LMCAS 中通常表示为 `n x 1` 的列矩阵。
- **点积 (Dot Product)**：通过转置矩阵相乘实现 `A^T * B`。
- **叉积 (Cross Product)**：仅适用于 3D 向量，位于 `include/symbolic_vector_geometry.hpp`。
