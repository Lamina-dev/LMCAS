# LMCAS Lamina 计算机代数系统

## 结构

```text
LMCAS/
├── src/                # 源代码实现 (.cpp)
│   ├── symbolic.cpp    # 符号表达式核心实现
│   ├── integration.cpp # 积分模块实现
│   ├── solver.cpp      # 方程求解器实现
│   └── ...
├── include/            # 头文件声明 (.hpp)
│   ├── symbolic.hpp    # 符号表达式类定义 (SymbolicExpr)
│   ├── integration.hpp # 积分模块接口
│   ├── solver.hpp      # 求解器接口
│   ├── value.hpp       # 通用值类型定义
│   ├── rational.hpp    # 有理数类
│   └── ...
├── tests/              # 单元测试
├── benchmarks/         # 性能测试
├── LAMMP/              # 子模块：Lamina MP (BigInt)
└── CMakeLists.txt      # CMake 构建配置
```

## 构建

### 要求
*   **CMake**: 3.14+
*   **编译器**: MSVC 2019+, GCC 9+, Clang 10+ (需支持 C++17)

### Visual Studio (Windows)

```powershell
mkdir build
cd build

cmake ..

cmake --build . --config Debug

# 运行测试
.\bin\Debug\test_proof.exe
```

注意：如果在运行可执行文件时遇到 DLL 缺失错误，请确保 (lib)LammpCore.dll 和 (lib)lmcas.dll 位于可执行文件同级目录，或将 uild/bin/Debug 添加到系统 PATH。

### Linux / macOS

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
./bin/test_proof
```

## 核心功能说明

### 符号系统 (include/symbolic.hpp)
SymbolicExpr 是 LMCAS 的核心类，表示一个不可变的符号表达式树。
- **创建表达式**: 使用静态工厂方法，如 SymbolicExpr::add, SymbolicExpr::variable, SymbolicExpr::number, SymbolicExpr::sin 等。
- **内存管理**: 基于 std::shared_ptr 的自动内存管理。
- **不可变性**: 所有对表达式的操作（如相加、求导）都会返回一个新的 SymbolicExpr 对象，原对象保持不变。

### 计算机代数算法
核心算法作为 SymbolicExpr 的成员函数或静态方法提供：
*   **求导**: expr->differentiate("x") - 对变量 x 求导。
*   **积分**: expr->integrate("x") - 对变量 x 进行符号积分。
*   **化简**: expr->simplify() - 调用化简引擎对表达式进行代数化简。
*   **展开**: expr->expand() - 展开多项式或乘积。
*   **代入**: expr->substitute("y", val) - 将变量 y 替换为表达式 val。

### 矩阵运算
支持符号矩阵操作（通常通过 SymbolicExpr 静态方法调用）：
*   SymbolicExpr::determinant(mat): 计算行列式。
*   SymbolicExpr::inverse(mat): 计算逆矩阵。
*   SymbolicExpr::eigenvalues(mat): 计算特征值。

### 数值系统
*   **BigInt**: 任意精度整数（基于 LAMMP）。
*   **Rational**: 任意精度有理数。
*   **Irrational**: 简单的无理数包装（如 sqrt(2), pi, e）。

## 文档
- [说明文档](docs/README.md)
- [核心架构](docs/modules/symbolic.md)

## 许可证
GNU Lesser General Public License v3.0 (LGPL-3.0)

## 贡献
- Lamina MP LAMMP - Jecricho Knox - Lamina-dev
- Lamina CAS - Ziyang Bai - Lamina-dev
- All contributors are contributed to the Lamina project.
