# LMCAS Lamina 计算机代数系统

## 结构

```text
LMCAS/
├── include/            # 已安装的 C++17 公共 API
├── src/                # 按 core/assumptions/algebra/calculus/analysis/facade 分层的实现
├── tests/              # LMCAS 回归、属性、包消费与头文件隔离测试
├── benchmarks/         # 可选基准
├── cmake/              # 安装包与结构约束
├── LMMC/               # C11 数值库；内部携带 LMMP
├── CMakePresets.json   # 受支持的开发、消毒器与安装配置
└── CMakeLists.txt
```

`src/` 中每个 `.cpp` 必须且只能属于一个组件目标。组件依赖按
`core → assumptions → algebra → calculus/linear algebra → analysis → facade`
方向声明；最终仅安装共享库 `lmcas`。

## 构建

### 要求

- CMake 3.26+
- Ninja
- 支持 C++17 和 C11 的 GCC、Clang 或 MinGW 工具链
- 初始化后的 `LMMC/LMMP` 子模块

### 开发与测试

```bash
cmake --preset strict-debug
cmake --build --preset strict-debug
ctest --preset strict-debug
```

该配置启用 LMCAS 与 LMMC 测试、严格警告和 `-Werror`。Linux 另提供
`linux-asan-ubsan` 与 `linux-tsan` 预设。

### 安装包

```bash
cmake --preset package
cmake --build --preset package
```

默认安装到 `build/package-install`。`tests/package_consumer` 验证安装后的
CMake targets、公共头文件隔离以及 C++/LMMC/LMMP 消费路径。

Windows 运行动态链接的测试或消费程序时，需让 `build/<preset>/bin`
（或安装目录的 `bin`）位于 `PATH`。

## 核心功能说明

所有 C++ 类型与函数均使用 `LMCAS` 命名空间，包括符号表达式、数值类型、
代数模板和 visitor。`include/expr.hpp` 提供 `LMCAS::Expr`、`LMCAS::sym`、
`LMCAS::parse_expr` 等表达式接口；公共导出宏为 `LMCAS_API`，
定义于 `include/lmcas_export.hpp`。消费方应使用这些库名称，而非规范名称。

### 符号系统 (include/symbolic.hpp)
LMCAS::SymbolicExpr 是 LMCAS 的核心类，表示一个不可变的符号表达式树。
- **创建表达式**: 使用静态工厂方法，如 LMCAS::SymbolicExpr::add, LMCAS::SymbolicExpr::variable, LMCAS::SymbolicExpr::number, LMCAS::SymbolicExpr::sin 等。
- **内存管理**: 基于 std::shared_ptr 的自动内存管理。
- **不可变性**: 所有对表达式的操作（如相加、求导）都会返回一个新的 LMCAS::SymbolicExpr 对象，原对象保持不变。

### 计算机代数算法
核心算法作为 LMCAS::SymbolicExpr 的成员函数或静态方法提供：
*   **求导**: expr->differentiate("x") - 对变量 x 求导。
*   **积分**: expr->integrate("x") - 对变量 x 进行符号积分。
*   **化简**: expr->simplify() - 调用化简引擎对表达式进行代数化简。
*   **展开**: expr->expand() - 展开多项式或乘积。
*   **代入**: expr->substitute("y", val) - 将变量 y 替换为表达式 val。

### 矩阵运算
符号矩阵操作通过 `include/symbolic_matrix.hpp` 的 checked 自由函数提供：
*   `LMCAS::matrix_determinant_checked`: 计算行列式。
*   `LMCAS::matrix_inverse_checked`: 计算逆矩阵。
*   `LMCAS::matrix_eigenvalues_checked`: 计算特征值。

### 数值系统
*   **LMCAS::BigInt**: 任意精度整数（基于 LMMP）。
*   **LMCAS::Rational**: 任意精度有理数。
*   **LMCAS::Irrational**: 简单的无理数包装（如 sqrt(2), pi, e）。

## 二进制兼容性

LMMC 提供固定 `double` 实数类型的 C ABI。LMCAS 的 C++ API 暴露
`std::string`、`std::vector` 与智能指针，因此不承诺跨编译器、标准库、
编译选项或补丁版本的二进制兼容性。CMake 包版本匹配采用精确版本；
升级 LMCAS 后应重新编译 C++ 消费方。

## 文档
- [符号表达式 API](include/symbolic.hpp)
- [表达式接口](include/expr.hpp)
- [符号矩阵 API](include/symbolic_matrix.hpp)

## 许可证
GNU Lesser General Public License v3.0 (LGPL-3.0)

## 贡献
- Lamina MP (LMMP) - Jecricho Knox - Lamina-dev
- Lamina CAS - Ziyang Bai - Lamina-dev
- All contributors are contributed to the Lamina project.
