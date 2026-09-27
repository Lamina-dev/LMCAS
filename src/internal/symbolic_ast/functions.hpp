/** @file internal/symbolic_ast/functions.hpp */
#pragma once
#include "internal/symbolic_ast/arithmetic.hpp"

namespace LMCAS {

/**
 * @brief 函数节点，表示数学函数调用（三角函数、对数、特殊函数等）。
 */
class FunctionNode : public SymbolicNode {
public:
    enum class FuncType {
        Sin, Cos, Tan, Cot, Sec, Csc,       /**< 三角函数。 */
        ArcSin, ArcCos, ArcTan,              /**< 反三角函数。 */
        Sinh, Cosh, Tanh,                    /**< 双曲函数。 */
        Ln, Log, Abs, Sqrt,                  /**< 对数、绝对值、平方根。 */
        Exp,                                 /**< 指数函数。 */
        LambertW,                            /**< Lambert W 函数。 */
        Atan2,                               /**< 双参数反正切。 */
        Infinity,                            /**< 无穷大。 */
        Erf,                                 /**< 误差函数 erf(x) = (2/√π)∫₀ˣ e^(-t²) dt。 */
        Ei,                                  /**< 指数积分 Ei(x)。 */
        Si,                                  /**< 正弦积分 Si(x) = ∫₀ˣ sin(t)/t dt。 */
        Ci,                                  /**< 余弦积分 Ci(x)。 */
        Li,                                  /**< 对数积分 Li(x) = ∫₀ˣ 1/ln(t) dt。 */
        Max,                                 /**< 最大值 max(a, b)。 */
        Min,                                 /**< 最小值 min(a, b)。 */
        Sgn,                                 /**< 符号函数 sgn(x) ∈ {-1, 0, 1}。 */
        Floor,                               /**< 下取整 ⌊x⌋。 */
        Ceil,                                /**< 上取整 ⌈x⌉。 */
        Round,                               /**< 四舍五入 round(x)。 */
        RealPart,                            /**< 实部 Re(z)。 */
        ImagPart,                            /**< 虚部 Im(z)。 */
        Conjugate,                           /**< 共轭 conj(z)。 */
        ComplexAbs,                          /**< 复数模 |z|。 */
        ComplexArg                           /**< 复数辐角 arg(z)。 */
    };

private:
    LMCAS_AST_NODE_FACTORY_FRIEND;

    const FuncType type_;
    const std::vector<std::shared_ptr<const SymbolicNode>> arguments_;

    FunctionNode(FuncType t, std::vector<std::shared_ptr<const SymbolicNode>> args)
        : type_(t), arguments_(std::move(args)) {
        for (const auto& arg : arguments_) {
            if (!arg) {
                throw std::invalid_argument("FunctionNode argument cannot be null");
            }
        }
    }

public:
    FuncType type() const noexcept { return type_; }
    const std::vector<std::shared_ptr<const SymbolicNode>>& arguments() const noexcept {
        return arguments_;
    }

    int type_priority() const override { return 2; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, static_cast<size_t>(type_));
        for (const auto& arg : arguments_) {
            hash_combine(seed, arg->hash());
        }
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const FunctionNode&>(other);
        if (type_ != o.type_) {
            return static_cast<int>(type_) < static_cast<int>(o.type_) ? -1 : 1;
        }
        if (arguments_.size() != o.arguments_.size()) {
            return arguments_.size() < o.arguments_.size() ? -1 : 1;
        }
        for (size_t i = 0; i < arguments_.size(); ++i) {
            int cmp = arguments_[i]->compare(*o.arguments_[i]);
            if (cmp != 0) return cmp;
        }
        return 0;
    }

public:
    void accept(LMCAS::detail::SymbolicVisitor& visitor) const override { LMCAS::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }
    std::shared_ptr<const SymbolicNode> clone() const override {
        std::vector<std::shared_ptr<const SymbolicNode>> new_args;
        for (const auto& arg : arguments_) new_args.push_back(arg->clone());
        return LMCAS::detail::make_node<FunctionNode>(type_, std::move(new_args));
    }
};
class UninterpretedFunctionNode : public SymbolicNode {
private:
    LMCAS_AST_NODE_FACTORY_FRIEND;
    const std::string name_;
    const std::vector<std::shared_ptr<const SymbolicNode>> arguments_;

    UninterpretedFunctionNode(
        std::string name,
        std::vector<std::shared_ptr<const SymbolicNode>> arguments)
        : name_(std::move(name)), arguments_(std::move(arguments)) {
        if (name_.empty()) {
            throw std::invalid_argument("uninterpreted function name cannot be empty");
        }
        for (const auto& argument : arguments_) {
            if (!argument) {
                throw std::invalid_argument("uninterpreted function arguments cannot be null");
            }
        }
    }

public:
    const std::string& name() const noexcept { return name_; }
    const auto& arguments() const noexcept { return arguments_; }
    int type_priority() const override { return 111; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = static_cast<std::size_t>(type_priority());
        hash_combine(seed, std::hash<std::string>{}(name_));
        for (const auto& argument : arguments_) hash_combine(seed, argument->hash());
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& function = static_cast<const UninterpretedFunctionNode&>(other);
        int comparison = name_.compare(function.name_);
        if (comparison != 0) return comparison;
        const auto count = std::min(arguments_.size(), function.arguments_.size());
        for (std::size_t index = 0; index < count; ++index) {
            comparison = arguments_[index]->compare(*function.arguments_[index]);
            if (comparison != 0) return comparison;
        }
        if (arguments_.size() == function.arguments_.size()) return 0;
        return arguments_.size() < function.arguments_.size() ? -1 : 1;
    }

public:
    void accept(LMCAS::detail::SymbolicVisitor& visitor) const override {
        LMCAS::detail::SymbolicVisitor::DepthGuard guard(visitor);
        visitor.visit(*this);
    }

    std::shared_ptr<const SymbolicNode> clone() const override {
        std::vector<std::shared_ptr<const SymbolicNode>> arguments;
        arguments.reserve(arguments_.size());
        for (const auto& argument : arguments_) arguments.push_back(argument->clone());
        return LMCAS::detail::make_node<UninterpretedFunctionNode>(
            name_, std::move(arguments));
    }
};

}
