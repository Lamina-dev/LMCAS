# LSR 支持边界

本文记录当前 LMCAS 与 LMMC 面向 LSR 的公开支持面。LSR 中属于 Lamina 解析器、类型检查器、模块加载器或运行时对象系统的条款，不在 LMCAS/LMMC 内伪实现。

## LMCAS

LMCAS 通过 `include/lsr_expr.hpp` 暴露 `lamina::lsr` 命名空间，当前覆盖：

- `Expr` 构造：符号、整数、有理数、显式近似实数、`pi`、`e`、`phi`、`i`、`I`、显式复数表达式；LSR 门面对旧 AST 中的变量名 `i`/`I` 在复数求值和等价规范化时按虚数单位处理。
- 基础 `Expr` 运算：`add`、`sub`、`mul`、`div`、`neg`、`eq` 构造表达式并返回 `Result`。
- 保留名：`i`、`I`、`pi`、`π`、`e`、`phi` 不能通过 `lamina::lsr::sym` 创建普通符号。
- `std.math` Expr 入口：`sqrt`、`pow`、`sin`、`cos`、`tan`、`asin`、`acos`、`atan`、`exp`、`log`、`log10`、`floor`、`ceil`、`round`、`clamp` 构造符号表达式并返回 `Result`。
- 安全变换入口：`simplify`、`expand`、`differentiate` 包装旧 CAS 变换并返回 `Result`，空输入、空变量、预算耗尽和旧 API 异常不会穿透到 Lamina 绑定层。
- 显式求值：`evalf` 与 `eval_complex` 返回 `Result`，未绑定符号、非有限值、预算耗尽和不支持表达式不会被静默近似。
- 复数边界：`real`、`imag`、`conj`、`abs` 支持显式复数表达式与可证明的实数输入；无法拆分的一般复函数返回不确定诊断。
- 集合结果：`ExprSet` 提供去重、成员、子集、并、交、差和对称差；`roots`、`solve`、`solve_expr_set` 将有限解集降为 `set<Expr>` 风格结果。
- 结构匹配：`expr_match` 暴露表达式模式匹配，支持显式通配符、函数/幂节点、重复通配符一致性、交换律加法/乘法匹配、确定性绑定顺序和 `Result` 错误传播。
- 数值塔域集合：`{Z}`、`{Q}`、`{R}`、`{C}`、`{Expr}` 的包含与子集关系在显式数值、显式复数表达式、LSR 虚数单位别名、由这些成员经加法/乘法/正整数幂构成的可证明复数表达式，以及任意 `Expr` 顶层成员上可判定。
- 数学等价：`equivalent_core` 支持 Core、Trig-Basic、ExpLog-Basic profile，并受显式预算限制。

## LMMC

LMMC 通过 `include/lmmc/lsr_stdlib.h` 暴露 C ABI 适配层，当前覆盖：

- `std.math`：`pi`、`e`、`phi`、初等实函数、复数 `i/I`、`complex`、`real`、`imag`、`conj`、复数 `abs`。
- table 键辅助：有限 `num`、有限 `complex`、`text` 与 `bool` 的相等和哈希，数值哈希规范化 `+0.0/-0.0`。
- `std.constants`：LSR-002 中列出的物理与化学常量，含名称、数值、单位和按索引枚举。
- `std.stats`：均值、中位数、方差、标准差、分位数、协方差和相关系数。
- `std.random`：显式 RNG 与默认 RNG 的种子、均匀采样、整数采样、正态采样和向量抽样。
- `std.units`：标准库形状的数值单位转换、通用单位转换、SI 基准剥离、当前刻度剥离、普通数值及单位字符串的无量纲判断和 LSR-008 错误名。
- `std.linalg`：shape、eye、diag、转置、伴随、行列式、逆、秩、迹、左右求解、eig/svd 表视图、向量与矩阵广播运算、关系广播。
- 诊断映射：`lmmc_lsr_error_name` 将 LMMC 状态码映射到 LSR 标准库诊断名。

## 非 LMCAS/LMMC 责任

以下 LSR 条款必须由 Lamina 编译器、运行时或包管理器实现：

- 词法与语法：复数字面量间距、`import/use`、模块路径、`match`、lambda、管道、推导式。
- 类型检查：集合元素统一、table 键类型约束、单位维度推导、函数签名匹配、`sym i` 的语言级报错。
- 运行时对象系统：vector、matrix、table、set 的引用计数、写时复制、`clone_if_shared`、扩展所有权注解。
- 标准 I/O：`std.io` 的文本输入输出和格式化。
- 扩展打包：`lampm.json`、`lib.lm` 加载、入口符号解析和接口绑定错误。

## 验证入口

当前安装与消费验证使用：

```sh
cmake --build build-dep-clean-mingw --target lmcas -- -j4
cmake --install build-dep-clean-mingw --prefix build-dep-clean-install
cmake -S tests/package_consumer -B build-dep-clean-consumer -DCMAKE_PREFIX_PATH="$PWD/build-dep-clean-install"
cmake --build build-dep-clean-consumer --config Debug -j4
```

LMMC LSR 适配层验证使用：

```sh
cmake --build build-lmmc-standalone-linked --target lmmc_test_lsr_stdlib -- -j4
ctest --test-dir build-lmmc-standalone-linked -R "lsr_stdlib" --output-on-failure
```
