#pragma once
#include "internal/symbolic_ast/relations.hpp"

namespace LMCAS {

/**
 * @brief 分段函数节点,表示条件分支表达式.
 *
 * 存储有序的 (表达式, 条件) 对列表和可选的默认表达式.
 * 条件应为 RelationalNode 或 LogicalNode 表达式.
 */
class PiecewiseNode : public SymbolicNode {
public:
    struct Branch {
        std::shared_ptr<const SymbolicNode> expression; /**< 分支值。 */
        std::shared_ptr<const SymbolicNode> condition;  /**< 条件（RelationalNode 或 LogicalNode）。 */
    };

private:
    LMCAS_AST_NODE_FACTORY_FRIEND;

    const std::vector<Branch> branches_;
    const std::shared_ptr<const SymbolicNode> default_expr_;

    /**
     * @brief 构造分段函数节点.
     * @param br 分支列表
     * @param def 默认表达式(可为 nullptr)
     */
    PiecewiseNode(std::vector<Branch> br, std::shared_ptr<const SymbolicNode> def = nullptr)
        : branches_(std::move(br)), default_expr_(std::move(def)) {
        if (branches_.empty()) {
            throw std::invalid_argument("PiecewiseNode requires at least one branch");
        }
        for (const auto& branch : branches_) {
            if (!branch.expression || !branch.condition) {
                throw std::invalid_argument("PiecewiseNode branch fields cannot be null");
            }
        }
    }

public:
    const std::vector<Branch>& branches() const noexcept { return branches_; }
    const std::shared_ptr<const SymbolicNode>& default_expr() const noexcept {
        return default_expr_;
    }

    int type_priority() const override { return 7; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        for (const auto& b : branches_) {
            hash_combine(seed, b.expression->hash());
            hash_combine(seed, b.condition->hash());
        }
        if (default_expr_) {
            hash_combine(seed, default_expr_->hash());
        }
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const PiecewiseNode&>(other);
        if (branches_.size() != o.branches_.size()) {
            return branches_.size() < o.branches_.size() ? -1 : 1;
        }
        for (size_t i = 0; i < branches_.size(); ++i) {
            int cmp = branches_[i].expression->compare(*o.branches_[i].expression);
            if (cmp != 0) return cmp;
            cmp = branches_[i].condition->compare(*o.branches_[i].condition);
            if (cmp != 0) return cmp;
        }
        bool has_def = (default_expr_ != nullptr);
        bool o_has_def = (o.default_expr_ != nullptr);
        if (has_def != o_has_def) return has_def ? 1 : -1;
        if (has_def && o_has_def) {
            return default_expr_->compare(*o.default_expr_);
        }
        return 0;
    }

public:
    void accept(LMCAS::detail::SymbolicVisitor& visitor) const override { LMCAS::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        std::vector<Branch> new_branches;
        new_branches.reserve(branches_.size());
        for (const auto& b : branches_) {
            new_branches.push_back({b.expression->clone(), b.condition->clone()});
        }
        auto new_def = default_expr_ ? default_expr_->clone() : nullptr;
        return LMCAS::detail::make_node<PiecewiseNode>(std::move(new_branches), std::move(new_def));
    }
};

/**
 * @brief 求和节点,表示符号有限/无限求和 sum_{k=a}^{b} f(k).
 */
class SummationNode : public SymbolicNode {
private:
    LMCAS_AST_NODE_FACTORY_FRIEND;

    const std::shared_ptr<const SymbolicNode> body_;
    const std::string index_var_;
    const std::shared_ptr<const SymbolicNode> lower_bound_;
    const std::shared_ptr<const SymbolicNode> upper_bound_;

    /**
     * @brief 构造求和节点.
     * @param b 通项表达式
     * @param idx 指标变量名
     * @param lo 下界
     * @param hi 上界
     */
    SummationNode(std::shared_ptr<const SymbolicNode> b, std::string idx,
                  std::shared_ptr<const SymbolicNode> lo, std::shared_ptr<const SymbolicNode> hi)
        : body_(std::move(b)), index_var_(std::move(idx)),
          lower_bound_(std::move(lo)), upper_bound_(std::move(hi)) {
        if (!body_ || !lower_bound_ || !upper_bound_) {
            throw std::invalid_argument("SummationNode children cannot be null");
        }
        if (index_var_.empty()) {
            throw std::invalid_argument("SummationNode index variable cannot be empty");
        }
    }

public:
    const std::shared_ptr<const SymbolicNode>& body() const noexcept { return body_; }
    const std::string& index_var() const noexcept { return index_var_; }
    const std::shared_ptr<const SymbolicNode>& lower_bound() const noexcept {
        return lower_bound_;
    }
    const std::shared_ptr<const SymbolicNode>& upper_bound() const noexcept {
        return upper_bound_;
    }

    int type_priority() const override { return 8; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, body_->hash());
        hash_combine(seed, std::hash<std::string>{}(index_var_));
        hash_combine(seed, lower_bound_->hash());
        hash_combine(seed, upper_bound_->hash());
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const SummationNode&>(other);
        int cmp = index_var_.compare(o.index_var_);
        if (cmp != 0) return cmp;
        cmp = lower_bound_->compare(*o.lower_bound_);
        if (cmp != 0) return cmp;
        cmp = upper_bound_->compare(*o.upper_bound_);
        if (cmp != 0) return cmp;
        return body_->compare(*o.body_);
    }

public:
    void accept(LMCAS::detail::SymbolicVisitor& visitor) const override { LMCAS::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return LMCAS::detail::make_node<SummationNode>(
            body_->clone(), index_var_, lower_bound_->clone(), upper_bound_->clone());
    }
};

/**
 * @brief 连乘节点,表示符号有限/无限连乘 prod_{k=a}^{b} f(k).
 */
class ProductNode : public SymbolicNode {
private:
    LMCAS_AST_NODE_FACTORY_FRIEND;

    const std::shared_ptr<const SymbolicNode> body_;
    const std::string index_var_;
    const std::shared_ptr<const SymbolicNode> lower_bound_;
    const std::shared_ptr<const SymbolicNode> upper_bound_;

    /**
     * @brief 构造连乘节点.
     * @param b 通项表达式
     * @param idx 指标变量名
     * @param lo 下界
     * @param hi 上界
     */
    ProductNode(std::shared_ptr<const SymbolicNode> b, std::string idx,
                   std::shared_ptr<const SymbolicNode> lo, std::shared_ptr<const SymbolicNode> hi)
        : body_(std::move(b)), index_var_(std::move(idx)),
          lower_bound_(std::move(lo)), upper_bound_(std::move(hi)) {
        if (!body_ || !lower_bound_ || !upper_bound_) {
            throw std::invalid_argument("ProductNode children cannot be null");
        }
        if (index_var_.empty()) {
            throw std::invalid_argument("ProductNode index variable cannot be empty");
        }
    }

public:
    const std::shared_ptr<const SymbolicNode>& body() const noexcept { return body_; }
    const std::string& index_var() const noexcept { return index_var_; }
    const std::shared_ptr<const SymbolicNode>& lower_bound() const noexcept {
        return lower_bound_;
    }
    const std::shared_ptr<const SymbolicNode>& upper_bound() const noexcept {
        return upper_bound_;
    }

    int type_priority() const override { return 9; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, body_->hash());
        hash_combine(seed, std::hash<std::string>{}(index_var_));
        hash_combine(seed, lower_bound_->hash());
        hash_combine(seed, upper_bound_->hash());
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const ProductNode&>(other);
        int cmp = index_var_.compare(o.index_var_);
        if (cmp != 0) return cmp;
        cmp = lower_bound_->compare(*o.lower_bound_);
        if (cmp != 0) return cmp;
        cmp = upper_bound_->compare(*o.upper_bound_);
        if (cmp != 0) return cmp;
        return body_->compare(*o.body_);
    }

public:
    void accept(LMCAS::detail::SymbolicVisitor& visitor) const override { LMCAS::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return LMCAS::detail::make_node<ProductNode>(
            body_->clone(), index_var_, lower_bound_->clone(), upper_bound_->clone());
    }
};

/**
 * @brief 积分变换节点,表示 Laplace,Fourier,Z 变换及其逆变换.
 *
 * 统一表示各类积分变换:L{f(t)}(s),F{f(t)}(omega),Z{f[n]}(z) 等.
 */
class TransformNode : public SymbolicNode {
public:
    enum class TransformType {
        Laplace,         /**< Laplace 变换 L{f(t)}(s)。 */
        InverseLaplace,  /**< 逆 Laplace 变换 L-¹{F(s)}(t)。 */
        Fourier,         /**< Fourier 变换 F{f(t)}(omega)。 */
        InverseFourier,  /**< 逆 Fourier 变换 F-¹{F(omega)}(t)。 */
        ZTransform       /**< Z 变换 Z{f[n]}(z)。 */
    };

private:
    LMCAS_AST_NODE_FACTORY_FRIEND;

    const TransformType transform_type_;
    const std::shared_ptr<const SymbolicNode> body_;
    const std::string source_var_;
    const std::shared_ptr<const SymbolicNode> target_;

    /**
     * @brief 构造积分变换节点.
     * @param tt 变换类型
     * @param b 被变换的表达式
     * @param src 源变量名
     * @param tgt 目标变量名
     */
    TransformNode(TransformType tt, std::shared_ptr<const SymbolicNode> b,
                  std::string src, std::shared_ptr<const SymbolicNode> target)
        : transform_type_(tt), body_(std::move(b)),
          source_var_(std::move(src)), target_(std::move(target)) {
        if (!body_ || !target_) {
            throw std::invalid_argument("TransformNode children cannot be null");
        }
        if (source_var_.empty()) {
            throw std::invalid_argument("TransformNode source variable cannot be empty");
        }
    }

public:
    TransformType transform_type() const noexcept { return transform_type_; }
    const std::shared_ptr<const SymbolicNode>& body() const noexcept { return body_; }
    const std::string& source_var() const noexcept { return source_var_; }
    const std::shared_ptr<const SymbolicNode>& target() const noexcept { return target_; }

    int type_priority() const override { return 11; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, static_cast<std::size_t>(transform_type_));
        hash_combine(seed, body_->hash());
        hash_combine(seed, std::hash<std::string>{}(source_var_));
        hash_combine(seed, target_->hash());
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const TransformNode&>(other);
        if (transform_type_ != o.transform_type_) {
            return static_cast<int>(transform_type_) < static_cast<int>(o.transform_type_) ? -1 : 1;
        }
        int cmp = source_var_.compare(o.source_var_);
        if (cmp != 0) return cmp;
        cmp = target_->compare(*o.target_);
        if (cmp != 0) return cmp;
        return body_->compare(*o.body_);
    }

public:
    void accept(LMCAS::detail::SymbolicVisitor& visitor) const override { LMCAS::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return LMCAS::detail::make_node<TransformNode>(
            transform_type_, body_->clone(), source_var_, target_->clone());
    }
};

}
