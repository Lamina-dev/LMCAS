# 如何使用 LMCAS

## 1. 安装与编译

1. 克隆仓库：
	```sh
	git clone https://github.com/yourname/LMCAS.git
	cd LMCAS
	```
2. 使用 CMake 构建（推荐 Release 模式）：
	```sh
	mkdir build && cd build
	cmake .. -DCMAKE_BUILD_TYPE=Release
	cmake --build .
	```
3. 运行测试（可选）：
	```sh
	ctest
	```

## 2. 基本用法

LMCAS 主要以 C++ 库形式使用，也可通过内置脚本接口或命令行工具调用。

### 2.1 C++ 代码集成

```cpp
#include "symbolic.hpp"
using namespace lmcas;

int main() {
	auto x = SymbolicExpr::variable("x");
	auto expr = SymbolicExpr::add(
		SymbolicExpr::add(
			SymbolicExpr::multiply(x, x),
			SymbolicExpr::multiply(SymbolicExpr::number(2), x)
		),
		SymbolicExpr::number(1)
	);
	auto d = expr->differentiate("x"); // 求导
	auto s = expr->integrate("x");    // 积分
	std::cout << "expr: " << expr->to_string() << std::endl;
	std::cout << "diff: " << d->to_string() << std::endl;
	std::cout << "integral: " << s->to_string() << std::endl;
	return 0;
}
```

> 注意：LMCAS 不支持直接用 C++ 运算符拼接表达式，必须用 SymbolicExpr::add/multiply/number/variable 等工厂函数。

### 2.2 命令行工具

部分构建目标（如 test_suite、benchmark_all）可直接运行，或自定义 main 入口。

```sh
./build/test_suite
```

### 2.3 脚本接口

支持 Value/LmModule/LmCppFunction 组合，便于嵌入 REPL 或脚本环境。

## 3. 常见问题

- 编译报错：请检查 C++17 支持，并确认已通过 `--recurse-submodules` 拉取 LAMMP 依赖。
- 表达式输出乱码：请确保终端支持 UTF-8。
- 更多问题见 [FAQ](faq.md)。

## 4. 进阶用法

- 参考 docs/modules/ 下各模块文档，了解符号运算、微积分、矩阵、复数等高级功能。
- 可扩展自定义函数、模块，支持高精度数值与符号混合。
# LMCAS 项目文档

本文档、本项目几乎所有注释和大部分单元测试使用ChatGPT-4.1自动化编写，“涵盖LMCAS项目的核心架构、功能说明、模块详解、开发指南等内容，旨在为用户和开发者提供全面的参考资料。”（这句话是AI自动补全的），经过人工核验真实性但并未人工修订，如有错误敬请谅解。

## 目录结构

- [快速入门](quickstart.md)
- [核心架构](architecture.md)
- [功能说明](features.md)
- [模块详解](modules/README.md)
- [开发与扩展](dev/README.md)
- [FAQ](faq.md)


## 简介
LMCAS（Lamina Computer Algebra System）是一个基于C++17的轻量级符号计算引擎，支持表达式化简、微积分、代数、矩阵、复数等功能。


## 文档树状结构

```
docs/
├── quickstart.md         # 快速入门
├── architecture.md       # 核心架构
├── features.md           # 功能说明
├── modules/
│   ├── README.md         # 模块总览
│   ├── symbolic.md       # 符号表达式模块
│   ├── calculus.md       # 微积分模块
│   ├── algebra.md        # 代数与方程模块
│   ├── matrix.md         # 矩阵与向量模块
│   ├── complex.md        # 复数模块
│   └── utils.md          # 工具与辅助模块
├── dev/
│   ├── README.md         # 开发与扩展总览
│   ├── build.md          # 构建与环境
│   ├── test.md           # 测试与调试
│   └── contribute.md     # 贡献指南
└── faq.md                # 常见问题
```


## 使用说明
请根据需求阅读对应模块文档，开发者可参考 dev/ 下的文档进行二次开发。
