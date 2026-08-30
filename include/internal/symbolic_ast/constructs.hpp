/** @file internal/symbolic_ast/constructs.hpp */
#pragma once
#include "relations.hpp"

/**
 * @brief 分段函数节点，表示条件分支表达式。
 *
 * 存储有序的 (表达式, 条件) 对列表和可选的默认表达式。
 * 条件应为 RelationalNode 或 LogicalNode 表达式。
 */
class PiecewiseNode : public SymbolicNode {
public:
    /** @brief 分支结构，包含表达式和对应条件。 */
    struct Branch {
        std::shared_ptr<const SymbolicNode> expression; ///< 分支值
        std::shared_ptr<const SymbolicNode> condition;  ///< 条件（RelationalNode 或 LogicalNode）
    };

private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::vector<Branch> branches_;
    const std::shared_ptr<const SymbolicNode> default_expr_;

    /**
     * @brief 构造分段函数节点。
     * @param br 分支列表
     * @param def 默认表达式（可为 nullptr）
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
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        std::vector<Branch> new_branches;
        new_branches.reserve(branches_.size());
        for (const auto& b : branches_) {
            new_branches.push_back({b.expression->clone(), b.condition->clone()});
        }
        auto new_def = default_expr_ ? default_expr_->clone() : nullptr;
        return lamina::detail::make_node<PiecewiseNode>(std::move(new_branches), std::move(new_def));
    }
};

/**
 * @brief 求和节点，表示符号有限/无限求和 ∑_{k=a}^{b} f(k)。
 */
class SummationNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::shared_ptr<const SymbolicNode> body_;
    const std::string index_var_;
    const std::shared_ptr<const SymbolicNode> lower_bound_;
    const std::shared_ptr<const SymbolicNode> upper_bound_;

    /**
     * @brief 构造求和节点。
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
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<SummationNode>(
            body_->clone(), index_var_, lower_bound_->clone(), upper_bound_->clone());
    }
};

/**
 * @brief 连乘节点，表示符号有限/无限连乘 ∏_{k=a}^{b} f(k)。
 */
class ProductNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::shared_ptr<const SymbolicNode> body_;
    const std::string index_var_;
    const std::shared_ptr<const SymbolicNode> lower_bound_;
    const std::shared_ptr<const SymbolicNode> upper_bound_;

    /**
     * @brief 构造连乘节点。
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
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<ProductNode>(
            body_->clone(), index_var_, lower_bound_->clone(), upper_bound_->clone());
    }
};

/**
 * @brief 积分变换节点，表示 Laplace、Fourier、Z 变换及其逆变换。
 *
 * 统一表示各类积分变换：L{f(t)}(s)、F{f(t)}(ω)、Z{f[n]}(z) 等。
 */
class TransformNode : public SymbolicNode {
public:
    /** @brief 变换类型枚举 */
    enum class TransformType {
        Laplace,         ///< Laplace 变换 L{f(t)}(s)
        InverseLaplace,  ///< 逆 Laplace 变换 L⁻¹{F(s)}(t)
        Fourier,         ///< Fourier 变换 F{f(t)}(ω)
        InverseFourier,  ///< 逆 Fourier 变换 F⁻¹{F(ω)}(t)
        ZTransform       ///< Z 变换 Z{f[n]}(z)
    };

private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const TransformType transform_type_;
    const std::shared_ptr<const SymbolicNode> body_;
    const std::string source_var_;
    const std::shared_ptr<const SymbolicNode> target_;

    /**
     * @brief 构造积分变换节点。
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
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<TransformNode>(
            transform_type_, body_->clone(), source_var_, target_->clone());
    }
};

/**
 * @brief 量词节点，表示全称量词 (∀) 或存在量词 (∃)。
 *
 * 用于表示逻辑公式中的量化表达式，如 ∀x∈D: P(x) 或 ∃x∈D: P(x)。
 */
class QuantifierNode : public SymbolicNode {
public:
    /** @brief 量词类型枚举 */
    enum class Type {
        ForAll, ///< 全称量词 ∀
        Exists  ///< 存在量词 ∃
    };

private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const Type quantifier_type_;
    const std::string bound_var_;
    const std::shared_ptr<const SymbolicNode> domain_;
    const std::shared_ptr<const SymbolicNode> predicate_;

    /**
     * @brief 构造量词节点。
     * @param qt 量词类型
     * @param var 约束变量名
     * @param dom 定义域表达式
     * @param pred 谓词表达式
     */
    QuantifierNode(Type qt, std::string var,
                   std::shared_ptr<const SymbolicNode> dom, std::shared_ptr<const SymbolicNode> pred)
        : quantifier_type_(qt), bound_var_(std::move(var)),
          domain_(std::move(dom)), predicate_(std::move(pred)) {
        if (!domain_ || !predicate_) {
            throw std::invalid_argument("QuantifierNode domain and predicate cannot be null");
        }
        if (bound_var_.empty()) {
            throw std::invalid_argument("QuantifierNode bound variable cannot be empty");
        }
    }

public:
    Type quantifier_type() const noexcept { return quantifier_type_; }
    const std::string& bound_var() const noexcept { return bound_var_; }
    const std::shared_ptr<const SymbolicNode>& domain() const noexcept { return domain_; }
    const std::shared_ptr<const SymbolicNode>& predicate() const noexcept {
        return predicate_;
    }

    int type_priority() const override { return 102; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, static_cast<std::size_t>(quantifier_type_));
        hash_combine(seed, std::hash<std::string>{}(bound_var_));
        hash_combine(seed, domain_->hash());
        hash_combine(seed, predicate_->hash());
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const QuantifierNode&>(other);
        if (quantifier_type_ != o.quantifier_type_) {
            return static_cast<int>(quantifier_type_) < static_cast<int>(o.quantifier_type_) ? -1 : 1;
        }
        int cmp = bound_var_.compare(o.bound_var_);
        if (cmp != 0) return cmp;
        cmp = domain_->compare(*o.domain_);
        if (cmp != 0) return cmp;
        return predicate_->compare(*o.predicate_);
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<QuantifierNode>(
            quantifier_type_, bound_var_, domain_->clone(), predicate_->clone());
    }
};

/**
 * @brief 集合构造器节点，表示集合构造式 {x ∈ D | P(x)}。
 *
 * 用于表示满足特定条件的元素集合。
 */
class SetBuilderNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;

    const std::string element_var_;
    const std::shared_ptr<const SymbolicNode> domain_;
    const std::shared_ptr<const SymbolicNode> predicate_;

    /**
     * @brief 构造集合构造器节点。
     * @param var 元素变量名
     * @param dom 定义域表达式
     * @param pred 成员条件表达式
     */
    SetBuilderNode(std::string var, std::shared_ptr<const SymbolicNode> dom,
                   std::shared_ptr<const SymbolicNode> pred)
        : element_var_(std::move(var)), domain_(std::move(dom)),
          predicate_(std::move(pred)) {
        if (!domain_ || !predicate_) {
            throw std::invalid_argument("SetBuilderNode domain and predicate cannot be null");
        }
        if (element_var_.empty()) {
            throw std::invalid_argument("SetBuilderNode element variable cannot be empty");
        }
    }

public:
    const std::string& element_var() const noexcept { return element_var_; }
    const std::shared_ptr<const SymbolicNode>& domain() const noexcept { return domain_; }
    const std::shared_ptr<const SymbolicNode>& predicate() const noexcept {
        return predicate_;
    }

    int type_priority() const override { return 103; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, std::hash<std::string>{}(element_var_));
        hash_combine(seed, domain_->hash());
        hash_combine(seed, predicate_->hash());
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const SetBuilderNode&>(other);
        int cmp = element_var_.compare(o.element_var_);
        if (cmp != 0) return cmp;
        cmp = domain_->compare(*o.domain_);
        if (cmp != 0) return cmp;
        return predicate_->compare(*o.predicate_);
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override { lamina::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<SetBuilderNode>(
            element_var_, domain_->clone(), predicate_->clone());
    }
};

/** A finite, unordered set of structurally distinct expressions. */
class FiniteSetNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;
    const std::vector<std::shared_ptr<const SymbolicNode>> elements_;

    explicit FiniteSetNode(std::vector<std::shared_ptr<const SymbolicNode>> elements)
        : elements_(canonicalize(std::move(elements))) {}

    static std::vector<std::shared_ptr<const SymbolicNode>> canonicalize(
        std::vector<std::shared_ptr<const SymbolicNode>> elements) {
        for (const auto& element : elements) {
            if (!element) throw std::invalid_argument("finite set elements cannot be null");
        }
        std::sort(elements.begin(), elements.end(),
                  [](const auto& lhs, const auto& rhs) { return lhs->compare(*rhs) < 0; });
        elements.erase(std::unique(elements.begin(), elements.end(),
                                   [](const auto& lhs, const auto& rhs) {
                                       return lhs->equals(*rhs);
                                   }),
                       elements.end());
        return elements;
    }

public:
    const std::vector<std::shared_ptr<const SymbolicNode>>& elements() const noexcept {
        return elements_;
    }
    bool contains(const SymbolicNode& element) const {
        return std::any_of(elements_.begin(), elements_.end(),
                           [&](const auto& candidate) { return candidate->equals(element); });
    }
    int type_priority() const override { return 104; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = static_cast<std::size_t>(type_priority());
        for (const auto& element : elements_) hash_combine(seed, element->hash());
        return seed;
    }
    int compare_same_type(const SymbolicNode& other) const override {
        const auto& set = static_cast<const FiniteSetNode&>(other);
        const auto count = std::min(elements_.size(), set.elements_.size());
        for (std::size_t index = 0; index < count; ++index) {
            const int comparison = elements_[index]->compare(*set.elements_[index]);
            if (comparison != 0) return comparison;
        }
        if (elements_.size() == set.elements_.size()) return 0;
        return elements_.size() < set.elements_.size() ? -1 : 1;
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override {
        lamina::detail::SymbolicVisitor::DepthGuard guard(visitor);
        visitor.visit(*this);
    }
    std::shared_ptr<const SymbolicNode> clone() const override {
        std::vector<std::shared_ptr<const SymbolicNode>> elements;
        elements.reserve(elements_.size());
        for (const auto& element : elements_) elements.push_back(element->clone());
        return lamina::detail::make_node<FiniteSetNode>(std::move(elements));
    }
};

/** A real interval with independently open or closed finite endpoints. */
class IntervalNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;
    const std::shared_ptr<const SymbolicNode> lower_;
    const std::shared_ptr<const SymbolicNode> upper_;
    const bool lower_closed_;
    const bool upper_closed_;

    IntervalNode(std::shared_ptr<const SymbolicNode> lower,
                 std::shared_ptr<const SymbolicNode> upper,
                 bool lower_closed, bool upper_closed)
        : lower_(std::move(lower)), upper_(std::move(upper)),
          lower_closed_(lower_closed), upper_closed_(upper_closed) {
        if (!lower_ || !upper_) throw std::invalid_argument("interval endpoints cannot be null");
    }

public:
    const auto& lower() const noexcept { return lower_; }
    const auto& upper() const noexcept { return upper_; }
    bool lower_closed() const noexcept { return lower_closed_; }
    bool upper_closed() const noexcept { return upper_closed_; }
    int type_priority() const override { return 105; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = static_cast<std::size_t>(type_priority());
        hash_combine(seed, lower_->hash());
        hash_combine(seed, upper_->hash());
        hash_combine(seed, lower_closed_);
        hash_combine(seed, upper_closed_);
        return seed;
    }
    int compare_same_type(const SymbolicNode& other) const override {
        const auto& interval = static_cast<const IntervalNode&>(other);
        int comparison = lower_->compare(*interval.lower_);
        if (comparison != 0) return comparison;
        comparison = upper_->compare(*interval.upper_);
        if (comparison != 0) return comparison;
        if (lower_closed_ != interval.lower_closed_) return lower_closed_ ? 1 : -1;
        if (upper_closed_ != interval.upper_closed_) return upper_closed_ ? 1 : -1;
        return 0;
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override {
        lamina::detail::SymbolicVisitor::DepthGuard guard(visitor);
        visitor.visit(*this);
    }
    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<IntervalNode>(
            lower_->clone(), upper_->clone(), lower_closed_, upper_closed_);
    }
};

/** A symbolic membership proposition. */
class MembershipNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;
    const std::shared_ptr<const SymbolicNode> element_;
    const std::shared_ptr<const SymbolicNode> set_;
    MembershipNode(std::shared_ptr<const SymbolicNode> element,
                   std::shared_ptr<const SymbolicNode> set)
        : element_(std::move(element)), set_(std::move(set)) {
        if (!element_ || !set_) throw std::invalid_argument("membership operands cannot be null");
    }

public:
    const auto& element() const noexcept { return element_; }
    const auto& set() const noexcept { return set_; }
    int type_priority() const override { return 106; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = static_cast<std::size_t>(type_priority());
        hash_combine(seed, element_->hash());
        hash_combine(seed, set_->hash());
        return seed;
    }
    int compare_same_type(const SymbolicNode& other) const override {
        const auto& membership = static_cast<const MembershipNode&>(other);
        const int comparison = element_->compare(*membership.element_);
        return comparison != 0 ? comparison : set_->compare(*membership.set_);
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override {
        lamina::detail::SymbolicVisitor::DepthGuard guard(visitor);
        visitor.visit(*this);
    }
    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<MembershipNode>(element_->clone(), set_->clone());
    }
};

/** A symbolic value paired with its canonical dimension and display scale. */
class QuantityNode : public SymbolicNode {
private:
    LAMINA_AST_NODE_FACTORY_FRIEND;
    const std::shared_ptr<const SymbolicNode> value_;
    const lamina::DimensionSignature dimension_;
    const Rational scale_to_base_;
    const std::string display_unit_;

    QuantityNode(std::shared_ptr<const SymbolicNode> value,
                 lamina::DimensionSignature dimension,
                 Rational scale_to_base, std::string display_unit)
        : value_(std::move(value)), dimension_(std::move(dimension)),
          scale_to_base_(std::move(scale_to_base)),
          display_unit_(std::move(display_unit)) {
        if (!value_) throw std::invalid_argument("quantity value cannot be null");
        if (scale_to_base_.is_zero()) throw std::invalid_argument("quantity scale cannot be zero");
    }

public:
    const auto& value() const noexcept { return value_; }
    const auto& dimension() const noexcept { return dimension_; }
    const Rational& scale_to_base() const noexcept { return scale_to_base_; }
    const std::string& display_unit() const noexcept { return display_unit_; }
    int type_priority() const override { return 107; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = static_cast<std::size_t>(type_priority());
        hash_combine(seed, value_->hash());
        hash_combine(seed, std::hash<std::string>{}(dimension_.to_string()));
        hash_combine(seed, scale_to_base_.hash());
        hash_combine(seed, std::hash<std::string>{}(display_unit_));
        return seed;
    }
    int compare_same_type(const SymbolicNode& other) const override {
        const auto& quantity = static_cast<const QuantityNode&>(other);
        int comparison = dimension_.to_string().compare(quantity.dimension_.to_string());
        if (comparison != 0) return comparison;
        if (scale_to_base_ < quantity.scale_to_base_) return -1;
        if (quantity.scale_to_base_ < scale_to_base_) return 1;
        comparison = display_unit_.compare(quantity.display_unit_);
        return comparison != 0 ? comparison : value_->compare(*quantity.value_);
    }

public:
    void accept(lamina::detail::SymbolicVisitor& visitor) const override {
        lamina::detail::SymbolicVisitor::DepthGuard guard(visitor);
        visitor.visit(*this);
    }
    std::shared_ptr<const SymbolicNode> clone() const override {
        return lamina::detail::make_node<QuantityNode>(
            value_->clone(), dimension_, scale_to_base_, display_unit_);
    }
};

#undef LAMINA_AST_NODE_FACTORY_FRIEND

inline std::shared_ptr<const SymbolicNode> SymbolicFactory::create_number(const ::BigInt& v) { return lamina::detail::make_node<NumberNode>(v); }
inline std::shared_ptr<const SymbolicNode> SymbolicFactory::create_number(const ::Rational& v) { return lamina::detail::make_node<NumberNode>(v); }
inline std::shared_ptr<const SymbolicNode> SymbolicFactory::create_number(lmmc_real_t v) { return lamina::detail::make_node<NumberNode>(v); }
inline std::shared_ptr<const SymbolicNode> SymbolicFactory::create_variable(const std::string& name) { return lamina::detail::make_node<VariableNode>(name); }

inline std::shared_ptr<const SymbolicNode> SymbolicFactory::create_add(std::vector<std::shared_ptr<const SymbolicNode>> ops) {
    if (ops.empty()) return create_number(::BigInt(0));
    for (const auto& op : ops) {
        if (!op) {
            throw std::invalid_argument("create_add operand cannot be null");
        }
    }
    if (ops.size() == 1) return ops[0];

    std::vector<std::shared_ptr<const SymbolicNode>> flat_ops;
    flat_ops.reserve(ops.size());
    for (const auto& op : ops) {
        if (op->is_zero()) continue;
        if (auto add = std::dynamic_pointer_cast<const AddNode>(op)) {
            flat_ops.insert(flat_ops.end(), add->operands().begin(), add->operands().end());
        } else {
            flat_ops.push_back(op);
        }
    }

    if (flat_ops.empty()) return create_number(::BigInt(0));
    if (flat_ops.size() == 1) return flat_ops[0];
    return lamina::detail::make_node<AddNode>(std::move(flat_ops));
}

inline std::shared_ptr<const SymbolicNode> SymbolicFactory::create_multiply(std::vector<std::shared_ptr<const SymbolicNode>> ops) {
    if (ops.empty()) return create_number(::BigInt(1));

    for (const auto& op : ops) {
        if (!op) {
            throw std::invalid_argument("create_multiply operand cannot be null");
        }
        if (op->is_zero()) return op;
    }

    if (ops.size() == 1) return ops[0];

    std::vector<std::shared_ptr<const SymbolicNode>> flat_ops;
    flat_ops.reserve(ops.size());
    for (const auto& op : ops) {
        if (op->is_one()) continue;
        if (auto mul = std::dynamic_pointer_cast<const MultiplyNode>(op)) {
            flat_ops.insert(flat_ops.end(), mul->operands().begin(), mul->operands().end());
        } else {
            flat_ops.push_back(op);
        }
    }

    if (flat_ops.empty()) return create_number(::BigInt(1));
    if (flat_ops.size() == 1) return flat_ops[0];
    return lamina::detail::make_node<MultiplyNode>(std::move(flat_ops));
}

inline std::shared_ptr<const SymbolicNode> SymbolicFactory::create_power(std::shared_ptr<const SymbolicNode> base, std::shared_ptr<const SymbolicNode> exponent) {
    if (!base || !exponent) {
        throw std::invalid_argument("create_power operands cannot be null");
    }
    if (exponent->is_zero()) {
        /// x^0 的化简需要 base 的定义域证明；base=0 时保留 PowerNode。
        if (base->is_zero()) {
            return lamina::detail::make_node<PowerNode>(std::move(base), std::move(exponent));
        }
        if (auto num = std::dynamic_pointer_cast<const NumberNode>(base)) {
            if (!num->is_zero()) return create_number(::BigInt(1));
        }
        if (auto complex = std::dynamic_pointer_cast<const ComplexNode>(base)) {
            auto real_num = std::dynamic_pointer_cast<const NumberNode>(complex->real());
            auto imag_num = std::dynamic_pointer_cast<const NumberNode>(complex->imag());
            if ((real_num && !real_num->is_zero()) || (imag_num && !imag_num->is_zero())) {
                return create_number(::BigInt(1));
            }
        }
        return lamina::detail::make_node<PowerNode>(std::move(base), std::move(exponent));
    }
    if (exponent->is_one()) return base;
    if (base->is_zero()) {
        /// 指数严格为正时 0^x 化简为 0，其余指数保留 PowerNode 及其定义域条件。
        if (exponent->is_positive()) return create_number(::BigInt(0));
        return lamina::detail::make_node<PowerNode>(std::move(base), std::move(exponent));
    }
    if (base->is_one()) return create_number(::BigInt(1));
    return lamina::detail::make_node<PowerNode>(std::move(base), std::move(exponent));
}

inline std::shared_ptr<const SymbolicNode> SymbolicFactory::create_complex(std::shared_ptr<const SymbolicNode> real, std::shared_ptr<const SymbolicNode> imag) {
    if (!real || !imag) {
        throw std::invalid_argument("create_complex operands cannot be null");
    }
    if (imag->is_zero()) {
        return real;
    }
    return lamina::detail::make_node<ComplexNode>(std::move(real), std::move(imag));
}
