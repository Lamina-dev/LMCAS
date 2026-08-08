# LMCAS 项目文档

LMCAS（Lamina Computer Algebra System）是一个 C++17 符号计算库，包含表达式、代数、微积分、矩阵、复数与 LSR 适配接口。LMMC 子模块提供数值、统计、随机与线性代数基础能力。

## 目录

- [快速入门](quickstart.md)
- [核心架构](architecture.md)
- [功能说明](features.md)
- [LSR 支持边界](lsr-support.md)
- [模块详解](modules/README.md)
- [开发与扩展](dev/README.md)
- [FAQ](faq.md)

## 使用说明

普通用户可从快速入门开始；需要对接 Lamina 语言或标准库时，先阅读 LSR 支持边界，确认相关行为属于 LMCAS、LMMC、Lamina 编译器还是 Lamina 运行时。
