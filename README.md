# LMCAS Lamina 计算机代数系统

## 结构

```text
LMCAS/
├── include/            # 已安装的 C++17 公共 API
├── src/                # 九个显式组件的分层实现
│   └── internal/       # AST、visitor 与算法私有头；不安装
├── tests/              # LMCAS 回归、属性、包消费与头文件隔离测试
├── benchmarks/         # 可选基准
├── cmake/              # 安装包与结构约束
├── LMMC/               # C11 数值库；内部携带 LMMP
├── CMakePresets.json   # 受支持的开发、消毒器与安装配置
└── CMakeLists.txt
```

`src/` 中每个 `.cpp` 必须且只能属于一个组件目标。组件从低层到高层为
`exact_numeric → core → assumptions → matrix_kernel → algebra → calculus/linear_algebra → analysis → facade`；
`calculus` 与 `linear_algebra` 同层。跨组件依赖只能指向更低层，
由编译后的符号引用和传递 include 图检查；最终仅安装共享库 `lmcas`。

## 构建

### 要求

- CMake 3.26+
- Ninja
- 支持 C++17 和 C11 的 GCC、Clang、AppleClang 或 MinGW 工具链；Windows 使用
  MinGW 或 GNU-driver Clang，MSVC ABI 前端（包括 clang-cl）不受支持
- 初始化后的 `LMMC/LMMP` 子模块
- 支持的平台：Windows、Linux 与 macOS

### 开发与测试

```bash
cmake --preset strict-debug
cmake --build --preset strict-debug
ctest --preset strict-debug
```

该配置启用 LMCAS 与 LMMC 测试、严格警告和 `-Werror`。Linux 另提供
`linux-asan-ubsan` 与 `linux-tsan` 预设。

全部预设串行复用 `build/strict-debug`，默认构建与测试并发为 2。
相同配置下直接运行增量构建，使用 `--target lmcas` 或具体测试目标缩小编译范围。
切换开发、消毒器和安装配置时重新配置同一目录；预设显式更新测试、消毒器和
汇编选项。嵌入 Lamina 时优先使用根工程的 Ninja 构建树，复用已有 CAS/MC 产物。
CAS 顶层使用 `LMCAS_LMMP_ASM` 控制内嵌数学栈；每次重新配置会向下同步当前选择，
更新 LMMP 的有效缓存值。

macOS 使用原生 AppleClang + Ninja。arm64 构建使用 `LMCAS_LMMP_ASM=AUTO`，
x86_64 构建使用 `LMCAS_LMMP_ASM=GENERIC`。编译级架构证据仅接受单一
`CMAKE_OSX_ARCHITECTURES=arm64` 或 `x86_64`；启用
`LMCAS_ENABLE_QUALITY_GATES` 时还需 Homebrew LLVM 的 `llvm-nm`（也接受
`llvm-nm-<版本>` 名称）：

```bash
cmake --preset strict-debug \
  -DCMAKE_OSX_ARCHITECTURES=arm64 \
  -DCMAKE_NM="$(brew --prefix llvm)/bin/llvm-nm" \
  -DLMCAS_ENABLE_QUALITY_GATES=ON
```

需要生成 universal binary 时可使用多架构 `CMAKE_OSX_ARCHITECTURES`，但须关闭
编译级证据门；该门要求每条第一方编译命令只有一个且一致的 `-arch`。

普通 CMake 配置默认关闭开发测试。自定义构建通过 `LMCAS_BUILD_TESTS=ON`
启用 LMCAS 测试，通过 `LMCAS_BUILD_LMMC_TESTS=ON` 启用内嵌 LMMC 测试。

| 测试范围 | 开发框架 | 自动获取版本 |
| --- | --- | --- |
| LMMC C11 测试 | cmocka | `cmocka-2.0.2` |
| LMCAS C++17 测试 | GoogleTest | `v1.18.0` |
| LMCAS 性质测试 | 官方 RapidCheck | `2c3c4365aca21ef4e612768fdc95b1ce8b39a651` |

CMake 优先查找已安装的框架包，缺失时获取上述固定版本。离线开发可设置
`FETCHCONTENT_SOURCE_DIR_CMOCKA`、`FETCHCONTENT_SOURCE_DIR_GOOGLETEST`、
`FETCHCONTENT_SOURCE_DIR_RAPIDCHECK`，指向对应的本地框架源码目录。
内部测试共享测试专用库 `lmcas_test_internal`；公共 API 测试继续链接产品库
`lmcas`。测试专用库单独编译源码来导出内部符号，产品库的符号可见性与安装内容
保持不变。Windows 的分配失败和线程析构测试复用生产组件目标文件，在测试
可执行文件内运行 LMCAS 代码，以观测本进程的分配和析构时序。自动获取的
GoogleTest 和 RapidCheck 构建为共享库；已安装框架沿用其包提供的链接方式。
MinGW 使用薄归档引用已有组件目标文件，减少调试信息的重复占用。
Windows 测试目标将依赖的 DLL 复制到可执行文件目录。
框架仅由测试可执行目标私有链接，产品安装、导出及客户消费依赖保持为
LMCAS、LMMC、LMMP 与工具链标准运行库。

CTest 自动发现 GoogleTest 用例，使用 `-L gtest` 筛选；LMMC 测试使用
`-L lmmc` 筛选，由 cmocka 报告组内用例。三个开发测试预设自动生成
`build/strict-debug/junit.xml`，在报告内保留完整用例名。普通 CTest 调用可通过
`--output-junit <path>` 指定统一 JUnit 报告。编译级架构反例使用标准库 `unittest`。

开发测试预设和 CI 设置 `RC_PARAMS=seed=42`，固定性质测试的随机序列。
性质失败会输出收缩后的反例及 `reproduce` 参数；将输出的
`RC_PARAMS="reproduce=..."` 传给对应测试可执行文件即可重放。

### 安装包

```bash
cmake --preset package
cmake --build --preset package
```

默认安装到 `build/package-install`。`tests/package_consumer` 验证安装后的
CMake targets、公共头文件隔离以及 C++/LMMC/LMMP 消费路径。

消费方使用 `find_package(LMCAS CONFIG REQUIRED)` 并链接 `LMCAS::lmcas`；
C++17、LMMC 的 C11 要求及 LMMP 链接依赖由 targets 传递。无需添加源码、
私有头或构建目录到 include 路径。

Windows 运行动态链接的测试或消费程序时，需让 `build/strict-debug/bin`
（或安装目录的 `bin`）位于 `PATH`。
macOS 安装库使用 `@rpath` install name。安装消费检查会移动整个安装前缀，
并在不注入 `DYLD_LIBRARY_PATH` 或 `DYLD_FALLBACK_LIBRARY_PATH` 的情况下
运行消费程序，以验证包不依赖原构建或安装位置。

## 核心功能说明

所有 C++ 类型与函数均使用 `LMCAS` 命名空间，包括符号表达式、数值类型、
代数模板和内部 visitor；visitor 不属于安装 API。`include/expr.hpp` 提供 `LMCAS::Expr`、`LMCAS::sym`、
`LMCAS::parse_expr` 等表达式接口；公共导出宏为 `LMCAS_API`，
定义于 `include/lmcas_export.hpp`。消费方应使用这些库名称，而非规范名称。

### 符号系统 (include/symbolic.hpp)
LMCAS::SymbolicExpr 是 LMCAS 的核心类，表示一个不可变的符号表达式树。
- **创建表达式**: 使用静态工厂方法，如 LMCAS::SymbolicExpr::add, LMCAS::SymbolicExpr::variable, LMCAS::SymbolicExpr::number, LMCAS::SymbolicExpr::sin 等。
- **内存管理**: 基于 std::shared_ptr 的自动内存管理。
- **不可变性**: 所有对表达式的操作（如相加、求导）都会返回一个新的 LMCAS::SymbolicExpr 对象，原对象保持不变。

### 计算机代数算法
表达式变换通过 `include/expr.hpp` 的自由函数返回 `ExprResult`。
多步计算传入同一个 `ComputationContext`，共享取消状态与累计工作预算：
*   **求导**: `LMCAS::differentiate(expr, "x", context)`。
*   **积分**: `integrator.integrate_checked(*expr, "x", context)`；定积分使用 `integrate_def_checked`。
*   **化简**: `LMCAS::simplify(expr, context)`。
*   **展开**: `LMCAS::expand(expr, context)`。
*   **代入**: `LMCAS::substitute(expr, "y", val, context)`；Binding 列表按顺序替换。
*   **四次方程**: 封闭式根的候选式含近似数时，默认的 `SolveOptions::return_rootof` 使用精确 `RootOf` 表示保留全部根。

每次访问 `value()` 前检查结果；失败时传播完整的 `CasError`。
变换的默认重载在最外层创建一次 context。`max_ast_nodes` 限制单个受检输入、
候选 AST 或输出森林的按引用展开节点数，共享子树的每条引用分别计数。
`max_steps` 限制跨步骤的累计工作量。形状缓存使用弱引用，生命周期限定于一次变换预算。

`equivalent(*lhs, *rhs, context, options)` 提供公开等价证明入口。
结构相等先通过取消和规模检查；不同结构且含 binary64 来源的输入返回
`Inconclusive`，残差检查返回 `UnprovedResidual`。精确零证书保留
BigInt/Rational 来源信息，`is_zero/is_one` 按字面数值严格比较。

积分器以 `derivative(candidate) - original_integrand` 的精确零残差证书接受闭式，
常数、线性拆分、分段表达式、内建/自定义策略和用户表项共享此条件。
候选未获证书时继续回退，同模式表项和仿射换元也按优先级继续匹配。
循环积分先消去私有占位符，再认证闭合结果；尚未解决的义务保留原始积分。
`SymbolicExpr::make_integral` 专门构造未求值积分，执行计算使用 checked 积分入口。
`abs(x)` 等实变量积分族通过 context 中的实数域声明获得所需事实；
先以 `assume_domain_checked("x", Domain::Real)` 声明，再用
`context.set_assumptions(...)` 关联假设。域条件待定时保留原始积分义务。

### 求积估计与变换条件

`quadrature_simpson_numeric`、`quadrature_gaussian_numeric` 与
`adaptive_simpson_numeric` 返回 `Result<QuadratureEstimate>`。
其中 `value` 是求积值，`estimated_error` 是由网格加密或相邻阶规则差异得到的
绝对误差估计；数值求值的 `ApproxReal::absolute_error` 继续表示认证界。
原 `numerical_integrate_numeric` 调用迁移到 `quadrature_simpson_numeric`。

Lamina 的三种 `std.cas.integrate_*` 返回 `(value, estimated_error)`，
分别通过 `.1`、`.2` 读取。Laplace、逆 Laplace、Fourier、逆 Fourier、
Z 变换与卷积返回 `(expression, conditions, roc)`；后两项是表达式数组，
保留结果成立所需的参数条件和收敛域约束。`std.cas.raw` 及相应子模块使用
同样的 tuple，并由 `Result` 承载失败。只消费变换表达式时显式取 `.1`，
同时保留其条件供后续计算使用。

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
编译选项或库更新之间的二进制兼容性。LMCAS 与 LMMC 当前均不设版本号，
不提供 CMake 包版本兼容性选择；更新 LMCAS 后仍应重新编译 C++ 消费方。

### LMCAS / LMMC API 迁移

升级后请重新编译调用方，并按下表替换旧接口：

| 旧 API | 新 API |
| --- | --- |
| `BigInt::ToString` / `Abs` / `IsNegative` | `to_string` / `abs` / `is_negative` |
| `BigInt::nPr` / `nCr` | `npr` / `ncr`，以及接受 `BigInt` 参数的 checked 入口 |
| `Rational::to_BigInt` | `to_bigint` |
| `SymbolicPolyCoeff::ToString` | `to_string` |
| `FGLMPoly::LM` / `LC` | `lead_monomial` / `lead_coeff` |
| `LMCAS::I()` | `LMCAS::imaginary_unit()` |
| 宽整数直接传给 `number` / `integer` | 显式构造 `BigInt` 后传入 |
| 无理数的旧机器整数平方根入口 | `Irrational::sqrt_checked(const BigInt&, ..., ComputationContext&)` |
| `Value::as_number` / `as_rational` / `as_irrational` / `as_symbolic` | 对应的 `*_checked` 方法并传播 `Result` 错误 |
| `Value::vector_add` / `dot_product` / `matrix_multiply` | 对应的 `*_checked` 方法并传播维度或输入错误 |
| `SymbolicExpr::get_type` 与嵌套 `Type` | 具体 predicate 或 `expr_match` |
| `SymbolicExpr::get_operands` | 公开变换与 matching API |
| `SymbolicExpr::get_number_value` | `is_number` + `get_number`，或 `evaluate_numeric` |
| `SymbolicExpr::get_identifier` | `symbol_name(expr)` 返回的借用 `optional<string_view>` |
| `SymbolicExpr::simplify/expand/differentiate/substitute` | `expr.hpp` 中对应的自由函数，检查 `ExprResult` |
| `AssumptionContext::simplify` | 通过 `context.set_assumptions(...)` 提供假设，再调用 `LMCAS::simplify(expr, context)` |
| `SymbolicExpr::integrate/integral`、`Integrator::integrate/integrate_def` | `Integrator::integrate_checked/integrate_def_checked`，显式传入 context |
| `equivalent_core` | `equivalent(lhs, rhs, context, options)` |
| `BuiltInIntegrationStrategy`、`PartialFractionStrategy` | 统一的 `IntegrationStrategy` 认证；有理分解由 `RationalDecompositionStrategy` 提供 |

数学整数与无理根的 radicand 使用 LMMP-backed `BigInt`；无理根的系数和
常数使用 `Rational`。窄化使用 `try_to_int64` / `try_to_uint64`，失败不饱和。
计数、容量、资源预算和随机状态继续使用有界机器类型。

AST 与 visitor 头统一在 `src/internal/`，不再提供公共 AST 聚合头。
LMMC 的 scalar/tensor 改名及窄子头见 [LMMC API 迁移](LMMC/README.md#lmmc-api-migration)；
语言层 `math.I` 和 tensor 操作名保持不变。

LMCAS 与 LMMC 不安装版本头文件或 CMake 包版本元数据。
ELF SONAME 分别为 `liblmcas.so` 和 `liblmmc.so`；macOS install name 分别为
`@rpath/liblmcas.dylib`、`@rpath/liblmmc.dylib`，LMMP 保留其版本身份
`@rpath/liblmmp.1.dylib`。LMMC 使用匿名导出节点，
保留公共符号白名单及其余符号的隐藏规则。LMMP 保留其独立版本。
配置时生成构建目录中的 `Doxyfile`，不设置项目版本号。
配置后可运行 `doxygen build/strict-debug/Doxyfile` 生成文档。

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
