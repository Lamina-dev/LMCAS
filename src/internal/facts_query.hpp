#pragma once

#include "assumption.hpp"
#include "result.hpp"
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace LMCAS {

class SymbolicNode;
class SymbolicExpr;
class ComputationContext;
class FactsQuery {
public:
    virtual Result<Tribool> is_positive(const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const = 0;
    virtual Result<Tribool> is_negative(const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const = 0;
    virtual Result<Tribool> is_nonnegative(const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const = 0;
    virtual Result<Tribool> is_nonzero(const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const = 0;
    virtual Result<Tribool> is_real(const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const = 0;

protected:
    ~FactsQuery() = default;
};

namespace detail {
const FactsQuery& no_facts() noexcept;

/**
 * @brief 为绑定体提供词法作用域，隔离同名自由变量的外层事实。
 * 无关表达式继续借用外层事实及其关系推论；栈上视图串联保留嵌套作用域。
 */
class ScopedFacts final : public FactsQuery {
public:
    ScopedFacts(const FactsQuery& outer, std::string_view bound_name)
        : outer_(outer), bound_name_(bound_name) {}

    Result<Tribool> is_positive(const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const override;
    Result<Tribool> is_negative(const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const override;
    Result<Tribool> is_nonnegative(const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const override;
    Result<Tribool> is_nonzero(const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const override;
    Result<Tribool> is_real(const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const override;

private:
    Result<bool> shadows(const std::shared_ptr<const SymbolicNode>& node, ComputationContext& context) const;
    const FactsQuery& outer_;
    const std::string bound_name_;
};
Result<Tribool> query_definedness(
    const std::shared_ptr<const SymbolicNode>& node, const FactsQuery& facts, Domain domain, ComputationContext& context);
Result<Tribool> query_nonzero_value(
    const std::shared_ptr<const SymbolicNode>& node, const FactsQuery& facts, Domain domain, ComputationContext& context);
Result<Tribool> query_real_value(
    const std::shared_ptr<const SymbolicNode>& node, const FactsQuery& facts, ComputationContext& context);
Result<Tribool> query_positive_value(
    const std::shared_ptr<const SymbolicNode>& node, const FactsQuery& facts, ComputationContext& context);
Result<Tribool> query_nonnegative_value(
    const std::shared_ptr<const SymbolicNode>& node, const FactsQuery& facts, ComputationContext& context);

/**
 * @brief 返回完整定义域约束：空向量表示已证明，nullopt 表示无法完整表达。
 * @note 非空向量表示条件合取；空定义域用假谓词表示。
 */
Result<std::optional<std::vector<std::shared_ptr<SymbolicExpr>>>> domain_constraints(
    const std::shared_ptr<const SymbolicNode>& expression, const FactsQuery& facts,
    Domain domain, ComputationContext& context);

}
}
