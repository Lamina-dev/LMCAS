# LMCAS Lamina 计算机代数系统

## 结构

```text
LMCAS/
├── src/                # 源代码实现 (.cpp)
│   ├── cas.cpp         # CAS 核心算法实现
│   ├── symbolic.cpp    # 符号表达式实现
│   └── ...
├── include/            # 头文件声明 (.hpp)
│   ├── symbolic.hpp    # 符号表达式类定义
│   ├── cas.hpp         # CAS 算法接口
│   ├── value.hpp       # 通用值类型定义
│   ├── rational.hpp    # 有理数类
│   └── ...
├── tests/              # 单元测试
├── benchmarks/         # 性能测试
├── LAMMP/              # 第三方依赖：大数运算库
└── CMakeLists.txt      # CMake 构建配置
```

## 构建


### 要求
*   **CMake**: 3.14+
*   **编译器**: MSVC 2019+, GCC 9+, Clang 10+

### Visual Studio


```powershell
mkdir build
cd build

cmake ..

cmake --build . --config Debug

.\bin\Debug\test_proof.exe
```

如果直接运行 `.exe` 失败，请确保 `LammpCore.dll` 和 `lmcas.dll` 位于可执行文件同级目录。

### Linux / macOS 

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
./bin/test_proof
```

## 说明

### 符号系统 (`include/symbolic.hpp`)
`SymbolicExpr` 是所有符号表达式的基类。使用 `std::shared_ptr` 以管理生命周期，通过静态工厂方法创建实例（如 `SymbolicExpr::add`, `SymbolicExpr::variable`）。
设计上采用了不可变数据结构，修改表达式会返回新的对象。

### 计算机代数算法 (`include/cas.hpp`)
`CAS` 类包含核心处理算法：
*   `diff(expr, var)`: 对表达式 `expr` 关于变量 `var` 求导。
*   `simplify(expr)`: 对表达式进行代数化简。

### 数值系统
*   **BigInt**: 基于 LAMMP 库的大整数实现。
*   **Rational**: 有理数分数表示。
*   **Irrational**: 简单的无理数包装。

### 错误
*   请忽略所有的 Intelisense 错误提示，那完全就是一个垃圾。

## 许可证
GNU Lesser General Public License v3.0 (LGPL-3.0)

## 贡献
Lamina MP LAMMP - Jecricho Knox - Lamina-dev
Lamina CAS - Ziyang Bai - Lamina-dev
All contributors are contributed to the Lamina project.
Here, I(ZiyangBai) extend my sincere gratitude and heartfelt thanks to them!