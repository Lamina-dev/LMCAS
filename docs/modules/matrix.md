# 说明

### matrix_multiply(A, B)
- 标准矩阵乘法，支持符号元素。
- 算法：按行列遍历，逐项相乘累加。
- 特性：支持符号/数值混合。
- 局限：大规模稀疏矩阵效率低。

### matrix_determinant(A)
- 行列式，递归展开或LU分解。
- 算法：小型矩阵递归展开，大型矩阵LU分解。
- 特性：支持符号矩阵。
- 局限：大规模稀疏矩阵效率低。

### matrix_inverse(A)
- 逆矩阵，伴随矩阵法或高斯消元。
- 算法：计算伴随矩阵，除以行列式；或高斯消元法。
- 特性：支持符号矩阵。
- 局限：行列式为零时无逆。

### matrix_rotation(theta, dim)
- 生成旋转矩阵，支持2/3维。
- 算法：构造标准旋转矩阵。
- 特性：支持符号角度。
- 局限：高维旋转未支持。

### matrix_reflection(angle, dim)
- 生成反射矩阵。
- 算法：构造标准反射矩阵。
- 特性：支持符号角度。
- 局限：高维反射未支持。

### matrix_scaling(sx, sy, dim)
- 生成缩放矩阵。
- 算法：构造标准缩放矩阵。
- 特性：支持符号缩放因子。
- 局限：高维缩放未支持。

### matrix_eigenvalues(A)
- 特征值，特征多项式法。
- 算法：构造特征多项式，求解根。
- 特性：支持符号矩阵。
- 局限：大型矩阵效率低。

### matrix_eigenvectors(A)
- 特征向量，解线性方程组。
- 算法：构造特征值对应方程组，求解。
- 特性：支持符号矩阵。
- 局限：大型矩阵效率低。
# 矩阵与向量模块

主要头文件：`include/symbolic_matrix.hpp`

## 主要函数

- `matrix_multiply(A, B)`：标准矩阵乘法，支持符号元素。
- `matrix_determinant(A)`：递归展开或LU分解。
- `matrix_inverse(A)`：伴随矩阵法或高斯消元。
- `matrix_rotation(theta, dim)`：生成旋转矩阵，支持2/3维。
- `matrix_reflection(angle, dim)`：生成反射矩阵。
- `matrix_scaling(sx, sy, dim)`：生成缩放矩阵。
- `matrix_eigenvalues(A)`：特征多项式法，支持小型矩阵。
- `matrix_eigenvectors(A)`：特征向量，解线性方程组。

## 特性

- 支持符号矩阵、数值矩阵混合。
- 支持高阶线性变换。

## 已知问题

- 大规模稀疏矩阵、特征向量算法效率低。
- 三维向量几何（如直线/平面交点）功能待补充。
