# 快速入门 (Quick Start)

本文档将指导你如何快速集成 LMCAS 到你的 C++ 项目中。

## 1. 环境准备

确保你的系统已安装：
- **CMake** (3.14 或更高版本)
- **C++ 编译器** (支持 C++17，如 GCC 9+, Clang 10+, MSVC 2019+)

## 2. 获取代码与构建

```bash
git clone https://github.com/your-username/LMCAS.git
cd LMCAS
mkdir build && cd build
cmake ..
cmake --build . --config Release
```

构建完成后，你会得到 `lmcas` 动态库（或静态库）以及 `LammpCore` 依赖库。

## 3. 编写第一个程序

创建一个名为 `demo.cpp` 的文件：

```cpp
#include "lmcas/include/symbolic.hpp" // 请根据你的 include 路径调整
#include <iostream>

using namespace lmcas; // 假设所有类都在 lmcas 或默认命名空间下（实际请检查 symbolic.hpp）

int main() {
    // 1. 创建变量 x
    auto x = SymbolicExpr::variable("x");

    // 2. 构建表达式: f(x) = x^2 + sin(x)
// 注意：LMCAS 目前主要通过静态工厂方法构建
    auto expr = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::sin(x)
    );

    std::cout << "Original: " << expr->to_string() << std::endl;

    // 3. 求导: f'(x)
    auto deriv = expr->differentiate("x");
    std::cout << "Derivative: " << deriv->to_string() << std::endl;

    // 4. 化简导数 (可选)
    auto simplified = deriv->simplify();
    std::cout << "Simplified: " << simplified->to_string() << std::endl;

    return 0;
}
```

## 4. 链接与运行

在 `CMakeLists.txt` 中链接 `lmcas`：

```cmake
add_executable(my_demo demo.cpp)
target_link_libraries(my_demo PRIVATE lmcas)
```

运行程序，预期输出：
```text
Original: (x^2 + sin(x))
Derivative: ((2 * x) + cos(x))
Simplified: 2*x + cos(x)
```
