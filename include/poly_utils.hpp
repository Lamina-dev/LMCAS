/**
 * @file poly_utils.hpp
 * @brief 符号表达式与多项式之间的转换工具，以及高斯消元。
 */
#pragma once

#include "symbolic.hpp"
#include "polynomial.hpp"
#include "symbolic_ast.hpp"

namespace lamina {

/**
 * @brief 将符号表达式转换为一元多项式
 * @tparam T 系数类型
 * @param expr 符号表达式
 * @param var 主变量名
 * @return 对应的一元多项式
 */
template <typename T>
Polynomial<T> symbolic_to_poly(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var);

/**
 * @brief 将一元多项式转换为符号表达式
 * @tparam T 系数类型
 * @param poly 一元多项式
 * @return 对应的符号表达式
 */
template <typename T>
std::shared_ptr<SymbolicExpr> poly_to_symbolic(const Polynomial<T>& poly);

/** @brief 符号表达式系数包装器，用于以符号表达式作为多项式系数 */
struct SymbolicPolyCoeff {
    std::shared_ptr<SymbolicExpr> val;  ///< 内部符号表达式

    /** @brief 默认构造，值为 0 */
    SymbolicPolyCoeff() : val(SymbolicExpr::number(0)) {}

    /**
     * @brief 从整数构造
     * @param v 整数值
     */
    explicit SymbolicPolyCoeff(int v) : val(SymbolicExpr::number(v)) {}

    /**
     * @brief 从符号表达式构造
     * @param v 符号表达式
     */
    SymbolicPolyCoeff(std::shared_ptr<SymbolicExpr> v) : val(std::move(v)) {}

    /** @brief 判等（通过化简差为零判断） */
    bool operator==(const SymbolicPolyCoeff& other) const {

        if (!val || !other.val) return false;

        if (val == other.val) return true;

        auto diff = SymbolicExpr::add(val, SymbolicExpr::multiply(other.val, SymbolicExpr::number(-1)));
        return diff->simplify()->is_zero();
    }

    bool operator!=(const SymbolicPolyCoeff& other) const {
        return !(*this == other);
    }

    /** @brief 加法 */
    SymbolicPolyCoeff operator+(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(SymbolicExpr::add(val, other.val));
    }

    /** @brief 减法 */
    SymbolicPolyCoeff operator-(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(
            SymbolicExpr::add(val, SymbolicExpr::multiply(other.val, SymbolicExpr::number(-1)))
        );
    }

    /** @brief 乘法 */
    SymbolicPolyCoeff operator*(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(SymbolicExpr::multiply(val, other.val));
    }

    /** @brief 除法 */
    SymbolicPolyCoeff operator/(const SymbolicPolyCoeff& other) const {
        return SymbolicPolyCoeff(SymbolicExpr::divide(val, other.val));
    }

    /** @brief 取负 */
    SymbolicPolyCoeff operator-() const {
        return SymbolicPolyCoeff(SymbolicExpr::multiply(val, SymbolicExpr::number(-1)));
    }

    /**
     * @brief 转换为字符串
     * @return 符号表达式的字符串表示
     */
    std::string ToString() const {
        return val ? val->to_string() : "0";
    }

    /** @brief 取绝对值（符号系数不做实际操作） */
    friend SymbolicPolyCoeff abs(const SymbolicPolyCoeff& s) {

        return s;
    }
};

/**
 * @brief 从符号表达式中提取指定类型的系数值
 * @tparam T 目标系数类型
 * @param c 符号表达式
 * @return 提取的系数值
 */
template<typename T>
T extract_coeff_value(const std::shared_ptr<SymbolicExpr>& c);

template<>
inline SymbolicPolyCoeff extract_coeff_value<SymbolicPolyCoeff>(const std::shared_ptr<SymbolicExpr>& c) {
    return SymbolicPolyCoeff(c);
}

template<>
inline BigInt extract_coeff_value<BigInt>(const std::shared_ptr<SymbolicExpr>& c) {
    auto simp = c->simplify();
    if (auto n = std::dynamic_pointer_cast<NumberNode>(simp->root)) {
        if (std::holds_alternative<BigInt>(n->value)) return std::get<BigInt>(n->value);
        if (std::holds_alternative<Rational>(n->value)) {

            Rational r = std::get<Rational>(n->value);
            if (r.is_integer()) return r.to_BigInt();
        }
        if (std::holds_alternative<lmmc_real_t>(n->value)) return BigInt((long long)std::get<lmmc_real_t>(n->value));
    }

    if (simp->is_zero()) return BigInt(0);

    if (simp->is_one()) return BigInt(1);

    return BigInt(0);
}

template<>
inline Rational extract_coeff_value<Rational>(const std::shared_ptr<SymbolicExpr>& c) {
    auto simp = c->simplify();
    if (auto n = std::dynamic_pointer_cast<NumberNode>(simp->root)) {
        if (std::holds_alternative<Rational>(n->value)) return std::get<Rational>(n->value);
        if (std::holds_alternative<BigInt>(n->value)) return Rational(std::get<BigInt>(n->value));
        if (std::holds_alternative<lmmc_real_t>(n->value)) return Rational::from_double(std::get<lmmc_real_t>(n->value));
    }
    if (simp->is_zero()) return Rational(0);
    if (simp->is_one()) return Rational(1);
    return Rational(0);
}

/**
 * @brief 判断符号节点是否依赖指定变量
 * @param node 符号节点
 * @param var 变量名
 * @return 包含该变量返回 true
 */
inline bool depends_on_var(const std::shared_ptr<SymbolicNode>& node, const std::string& var) {
    if (!node) return false;

    struct DepVisitor : public SymbolicVisitor {
        std::string v;
        bool found = false;
        void visit(NumberNode&) override {}
        void visit(VariableNode& n) override { if (n.name == v) found = true; }
        void visit(AddNode& n) override { for(auto& op : n.operands) { if(found) return; op->accept(*this); } }
        void visit(MultiplyNode& n) override { for(auto& op : n.operands) { if(found) return; op->accept(*this); } }
        void visit(PowerNode& n) override { n.base->accept(*this); if(!found) n.exponent->accept(*this); }
        void visit(FunctionNode& n) override { for(auto& arg : n.arguments) { if(found) return; arg->accept(*this); } }
        void visit(MatrixNode& n) override {
            if (std::holds_alternative<MatrixNode::DenseStorage>(n.storage)) {
                for (auto& item : std::get<MatrixNode::DenseStorage>(n.storage)) {
                    if (item) { item->accept(*this); if (found) return; }
                }
            } else {
                for (auto& [idx, item] : std::get<MatrixNode::SparseStorage>(n.storage)) {
                    item->accept(*this); if (found) return;
                }
            }
        }
    } visitor;
    visitor.v = var;
    node->accept(visitor);
    return visitor.found;
}

/**
 * @brief 判断符号表达式是否包含指定变量
 * @param expr 符号表达式
 * @param var 变量名
 * @return 包含该变量返回 true
 */
inline bool contains(const SymbolicExpr& expr, const std::string& var) {
    return depends_on_var(expr.root, var);
}

/**
 * @brief 递归地将符号节点转换为一元多项式
 * @tparam T 系数类型
 * @param node 符号节点
 * @param var 主变量名
 * @return 对应的一元多项式
 */
template <typename T>
Polynomial<T> symbolic_to_poly_recursive(const std::shared_ptr<SymbolicNode>& node, const std::string& var) {
    if (!node) return Polynomial<T>(var);

    if (!depends_on_var(node, var)) {
        return Polynomial<T>({extract_coeff_value<T>(std::make_shared<SymbolicExpr>(node))}, var);
    }

    if (auto v = std::dynamic_pointer_cast<VariableNode>(node)) {
        if (v->name == var) {
            return Polynomial<T>({T(0), T(1)}, var);
        }
    }

    if (auto add = std::dynamic_pointer_cast<AddNode>(node)) {
        Polynomial<T> res(var);
        for (auto& op : add->operands) {
            res = res + symbolic_to_poly_recursive<T>(op, var);
        }
        return res;
    }

    if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(node)) {
        Polynomial<T> res({T(1)}, var);
        for (auto& op : mul->operands) {
            res = res * symbolic_to_poly_recursive<T>(op, var);
        }
        return res;
    }

    if (auto pow = std::dynamic_pointer_cast<PowerNode>(node)) {
        if (auto exp_num = std::dynamic_pointer_cast<NumberNode>(pow->exponent)) {
            int e_val = 0;
            if (std::holds_alternative<BigInt>(exp_num->value)) e_val = std::get<BigInt>(exp_num->value).to_int();
            else if (std::holds_alternative<double>(exp_num->value)) e_val = (int)std::get<double>(exp_num->value);
            else if (std::holds_alternative<Rational>(exp_num->value)) e_val = (int)std::get<Rational>(exp_num->value).to_double();

            if (e_val == 0) return Polynomial<T>({T(1)}, var);
            if (e_val > 0) {
                auto base_poly = symbolic_to_poly_recursive<T>(pow->base, var);
                Polynomial<T> res({T(1)}, var);
                for (int i = 0; i < e_val; ++i) res = res * base_poly;
                return res;
            }
        }
    }

    return Polynomial<T>(var);
}

template <typename T>
Polynomial<T> symbolic_to_poly(const std::shared_ptr<SymbolicExpr>& expr, const std::string& var) {
    if (!expr || !expr->root) return Polynomial<T>(var);
    return symbolic_to_poly_recursive<T>(expr->root, var);
}

template <typename T>
std::shared_ptr<SymbolicExpr> poly_to_symbolic(const Polynomial<T>& poly) {
    if (poly.is_zero()) return SymbolicExpr::number(0);

    std::vector<std::shared_ptr<SymbolicExpr>> terms;

    for (size_t i = 0; i < poly.coeffs.size(); ++i) {
        if (poly.coeffs[i] == T(0)) continue;

        auto coeff_node = SymbolicExpr::number(0);

        if constexpr (std::is_same_v<T, BigInt>) {
            coeff_node = SymbolicExpr::number(poly.coeffs[i]);
        } else if constexpr (std::is_same_v<T, Rational>) {
            coeff_node = SymbolicExpr::number(poly.coeffs[i]);
        } else {

            coeff_node = SymbolicExpr::number(poly.coeffs[i]);
        }

        if (i == 0) {
            terms.push_back(coeff_node);
        } else {

            auto var_node = std::make_shared<SymbolicExpr>(std::make_shared<VariableNode>(poly.variable_name));
            std::shared_ptr<SymbolicExpr> var_part;
            if (i == 1) {
                var_part = var_node;
            } else {
                var_part = SymbolicExpr::power(var_node, SymbolicExpr::number((int)i));
            }

            if (poly.coeffs[i] == T(1)) {
                terms.push_back(var_part);
            } else if (poly.coeffs[i] == T(-1)) {
                 terms.push_back(SymbolicExpr::multiply(SymbolicExpr::number(-1), var_part));
            } else {
                terms.push_back(SymbolicExpr::multiply(coeff_node, var_part));
            }
        }
    }

    std::reverse(terms.begin(), terms.end());

    if (terms.empty()) return SymbolicExpr::number(0);
    if (terms.size() == 1) return terms[0];

    auto res = terms[0];
    for (size_t i = 1; i < terms.size(); ++i) {
        res = SymbolicExpr::add(res, terms[i]);
    }
    return res;
}

/**
 * @brief 对增广矩阵进行高斯消元
 * @param A 增广矩阵（符号表达式元素）
 * @param m 行数
 * @param n 列数
 * @param pivot_col_for_row 输出各行的主元列号
 * @param sign 输出行交换的符号（+1 或 -1）
 */
void gaussian_eliminate(std::vector<std::vector<std::shared_ptr<SymbolicExpr>>>& A, size_t m, size_t n, std::vector<size_t>& pivot_col_for_row, int& sign);

}
