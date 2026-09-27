#pragma once
#include "internal/symbolic_ast/bindings.hpp"

namespace LMCAS {
class QuantifierNode : public SymbolicNode {
public:
    enum class Type {
        ForAll, /**< 全称量词 forall。 */
        Exists  /**< 存在量词 exists。 */
    };

private:
    LMCAS_AST_NODE_FACTORY_FRIEND;

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
    void accept(LMCAS::detail::SymbolicVisitor& visitor) const override { LMCAS::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return LMCAS::detail::make_node<QuantifierNode>(
            quantifier_type_, bound_var_, domain_->clone(), predicate_->clone());
    }
};
class SetBuilderNode : public SymbolicNode {
private:
    LMCAS_AST_NODE_FACTORY_FRIEND;

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
    void accept(LMCAS::detail::SymbolicVisitor& visitor) const override { LMCAS::detail::SymbolicVisitor::DepthGuard guard(visitor); visitor.visit(*this); }

    std::shared_ptr<const SymbolicNode> clone() const override {
        return LMCAS::detail::make_node<SetBuilderNode>(
            element_var_, domain_->clone(), predicate_->clone());
    }
};
class FiniteSetNode : public SymbolicNode {
private:
    LMCAS_AST_NODE_FACTORY_FRIEND;
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
    void accept(LMCAS::detail::SymbolicVisitor& visitor) const override {
        LMCAS::detail::SymbolicVisitor::DepthGuard guard(visitor);
        visitor.visit(*this);
    }
    std::shared_ptr<const SymbolicNode> clone() const override {
        std::vector<std::shared_ptr<const SymbolicNode>> elements;
        elements.reserve(elements_.size());
        for (const auto& element : elements_) elements.push_back(element->clone());
        return LMCAS::detail::make_node<FiniteSetNode>(std::move(elements));
    }
};
class IntervalNode : public SymbolicNode {
private:
    LMCAS_AST_NODE_FACTORY_FRIEND;
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
    int compare_same_type(const SymbolicNode& other) const override;

public:
    void accept(LMCAS::detail::SymbolicVisitor& visitor) const override {
        LMCAS::detail::SymbolicVisitor::DepthGuard guard(visitor);
        visitor.visit(*this);
    }
    std::shared_ptr<const SymbolicNode> clone() const override {
        return LMCAS::detail::make_node<IntervalNode>(
            lower_->clone(), upper_->clone(), lower_closed_, upper_closed_);
    }
};
class MembershipNode : public SymbolicNode {
private:
    LMCAS_AST_NODE_FACTORY_FRIEND;
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
    void accept(LMCAS::detail::SymbolicVisitor& visitor) const override {
        LMCAS::detail::SymbolicVisitor::DepthGuard guard(visitor);
        visitor.visit(*this);
    }
    std::shared_ptr<const SymbolicNode> clone() const override {
        return LMCAS::detail::make_node<MembershipNode>(element_->clone(), set_->clone());
    }
};
class QuantityNode : public SymbolicNode {
private:
    LMCAS_AST_NODE_FACTORY_FRIEND;
    const std::shared_ptr<const SymbolicNode> value_;
    const LMCAS::DimensionSignature dimension_;
    const Rational scale_to_base_;
    const std::string display_unit_;

    QuantityNode(std::shared_ptr<const SymbolicNode> value,
                 LMCAS::DimensionSignature dimension,
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
    void accept(LMCAS::detail::SymbolicVisitor& visitor) const override {
        LMCAS::detail::SymbolicVisitor::DepthGuard guard(visitor);
        visitor.visit(*this);
    }
    std::shared_ptr<const SymbolicNode> clone() const override {
        return LMCAS::detail::make_node<QuantityNode>(
            value_->clone(), dimension_, scale_to_base_, display_unit_);
    }
};

}
