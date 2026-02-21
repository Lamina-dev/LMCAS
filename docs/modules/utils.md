# 说明

### Value 类型推断与混合运算
- 基于 std::variant 自动推断类型，支持多种数值/符号混合。
- 算法：构造时自动识别类型，运算时按类型分派。
- 特性：支持符号与数值混合运算。
- 局限：嵌套容器类型推断复杂。

### 容器嵌套
- 支持 std::set、std::vector、std::vector<vector> 等嵌套。
- 算法：递归遍历容器，自动推断元素类型。
- 特性：支持矩阵、集合等复杂结构。
- 局限：深度嵌套时性能瓶颈。

### LmModule
- 模块化结构，类似 Python dict。
- 算法：内部用 unordered_map 存储子项，支持动态添加/删除。
- 特性：支持嵌套模块，便于脚本化扩展。
- 局限：深度嵌套时查找效率下降。

### LmCppFunction
- C++函数封装，支持 Value 类型参数与返回值。
- 算法：内部用 std::function<Value(std::vector<Value>)>。
- 特性：支持高阶函数、匿名函数。
- 局限：类型推断需手动管理。

### 辅助特性
- 表达式打印、调试输出、数值类型转换。
- 算法：to_string 自动识别类型并格式化输出。
- 特性：支持多种类型自动转换。
- 局限：边界类型未完全覆盖。

# 工具与辅助模块

主要头文件：`include/value.hpp`

## Value 类

- 通用值类型，采用 `std::variant` 实现，支持多种数据结构：
	- 基本类型：`int`, `double`, `bool`, `std::string`
	- 高精度类型：`BigInt`, `Rational`, `Irrational`
	- 符号类型：`std::shared_ptr<SymbolicExpr>`
	- 容器类型：`std::set<Value>`, `std::vector<Value>`, `std::vector<std::vector<Value>>`
	- 模块与函数：`std::shared_ptr<LmModule>`, `std::shared_ptr<LmCppFunction>`
	- 结构体与Lambda：`std::shared_ptr<lmStruct>`, `std::shared_ptr<LambdaDeclExpr>`

- 构造函数重载，自动推断类型。
- 算法：类型安全的自动转换，支持 `as_number()`、`as_symbolic()`、`is_numeric()` 等类型判断与转换。
- 特性：
	- 支持符号与数值混合运算。
	- 支持嵌套容器（如矩阵、集合）。
	- 支持脚本式调用（如模块、函数类型）。
- 已知问题：
	- 类型推断复杂，部分边界类型（如嵌套容器）未完全覆盖。
	- 性能瓶颈主要在频繁的 variant 访问与类型转换。

## LmModule

- 模块化类型，支持类似 Python dict 的结构。
- 算法：内部采用 `std::unordered_map<std::string, Value>` 存储子项。
- 特性：
	- 支持嵌套模块，便于脚本式扩展。
	- 支持动态添加/删除子项。
- 已知问题：
	- 模块嵌套过深时，查找效率下降。

## LmCppFunction

- C++函数封装，支持 Value 类型参数与返回值。
- 算法：内部采用 `std::function<Value(std::vector<Value>)>`。
- 特性：
	- 支持高阶函数、匿名函数。
	- 可与 SymbolicExpr 混合调用。
- 已知问题：
	- 函数类型不支持类型推断，需手动管理参数类型。

## 辅助特性

- 表达式打印：支持 to_string()，自动识别类型并格式化输出。
- 调试输出：支持 err_stream，便于日志与调试。
- 数值类型转换：支持 as_number()、as_symbolic()，自动处理符号与数值类型。

## 已知问题汇总

- 类型系统复杂，部分边界类型未完全覆盖。
- 性能瓶颈主要在 variant 访问与递归容器。
- 脚本式模块与函数扩展能力有限，需进一步完善。
