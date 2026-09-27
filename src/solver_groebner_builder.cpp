#include "internal/solver_groebner_builder.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

namespace LMCAS::groebner_detail {
namespace {
std::optional<BigInt> integer_exponent(const PowerNode& node) {
    if (!node.exponent()->is_number()) {
        return std::nullopt;
    }
    auto number = std::dynamic_pointer_cast<const NumberNode>(node.exponent());
    if (!number) {
        return std::nullopt;
    }
    const auto& value = number->value();
    if (std::holds_alternative<BigInt>(value)) {
        return std::get<BigInt>(value);
    }
    if (std::holds_alternative<Rational>(value)) {
        const auto& rational = std::get<Rational>(value);
        if (rational.is_integer()) {
            return rational.to_bigint();
        }
    }
    if (std::holds_alternative<lmmc_real_t>(value)) {
        lmmc_real_t real = std::get<lmmc_real_t>(value);
        if (std::isfinite(real) && real == std::floor(real)) {
            return Rational::from_double(real).to_bigint();
        }
    }
    return std::nullopt;
}

bool power_degrees_fit(const Poly& base, int exponent) {
    for (const auto& [monomial, coefficient] : base.terms) {
        (void)coefficient;
        for (int degree : monomial) {
            if (degree > std::numeric_limits<int>::max() / exponent) return false;
        }
    }
    return true;
}
}

    PolyBuilder::PolyBuilder(const std::vector<std::string>& v, PolyContext& ctx, bool strict)
        : vars(v), ext_vars(ctx.ext_vars), transcendental_map(ctx.transcendental_map),
          aux_to_node(&ctx.aux_to_node), result(ctx.ext_vars.size()), strict_mode(strict) {}

        void PolyBuilder::visit(const NumberNode& node) {
            result = Poly(ext_vars.size());

            if (std::holds_alternative<Rational>(node.value())) {
                result.add_term(std::vector<int>(ext_vars.size(), 0), std::get<Rational>(node.value()));
            } else if (std::holds_alternative<BigInt>(node.value())) {
                result.add_term(std::vector<int>(ext_vars.size(), 0), Rational(std::get<BigInt>(node.value())));
            } else if (std::holds_alternative<lmmc_real_t>(node.value())) {
                result.add_term(std::vector<int>(ext_vars.size(), 0),
                    Rational::from_double(std::get<lmmc_real_t>(node.value())));
            }
        }

        void PolyBuilder::visit(const VariableNode& node) {
            result = Poly(ext_vars.size());

            auto it = std::find(ext_vars.begin(), ext_vars.end(), node.name());
            if (it != ext_vars.end()) {
                Monomial m(ext_vars.size(), 0);
                m[std::distance(ext_vars.begin(), it)] = 1;
                result.add_term(m, Rational(1));
            } else {

                if (strict_mode) {
                    failed = true;
                    return;
                }
                size_t idx = get_or_create_aux_var(LMCAS::detail::make_node<VariableNode>(node.name()));

                result = Poly(ext_vars.size());
                Monomial m(ext_vars.size(), 0);
                m[idx] = 1;
                result.add_term(m, Rational(1));
            }
        }

        void PolyBuilder::visit(const AddNode& node) {
            Poly sum(ext_vars.size());
            for (auto& op : node.operands()) {
                PolyBuilder b(vars, ext_vars, transcendental_map, aux_to_node, strict_mode);
                op->accept(b);
                if (b.failed) { failed = true; return; }

                Poly b_res = b.get_result();
                align_variables(sum, b_res);
                sum = add_poly(sum, b_res);
            }
            result = sum;
        }

        void PolyBuilder::visit(const MultiplyNode& node) {
            Poly prod(ext_vars.size());
            prod.add_term(std::vector<int>(ext_vars.size(), 0), Rational(1));

            for (auto& op : node.operands()) {
                PolyBuilder b(vars, ext_vars, transcendental_map, aux_to_node, strict_mode);
                op->accept(b);
                if (b.failed) { failed = true; return; }
                Poly b_res = b.get_result();

                align_variables(prod, b_res);
                prod = mul_poly(prod, b_res);
            }
            result = prod;
        }

void PolyBuilder::visit(const PowerNode& node) {
    PolyBuilder builder(vars, ext_vars, transcendental_map, aux_to_node, strict_mode);
    node.base()->accept(builder);
    if (builder.failed) {
        failed = true;
        return;
    }
    Poly base = builder.get_result();
    auto exponent = integer_exponent(node);
    if (!exponent || *exponent > BigInt(std::numeric_limits<int>::max()) ||
        *exponent < 0) {
        if (strict_mode) {
            failed = true;
            return;
        }
        represent_as_aux_or_fail(
            LMCAS::detail::make_node<PowerNode>(node.base(), node.exponent()));
        return;
    }
    if (*exponent == 0) {
        result = Poly(ext_vars.size());
        result.add_term(Monomial(ext_vars.size(), 0), Rational(1));
        return;
    }
    const int count = static_cast<int>(exponent->try_to_int64().value());
    if (!power_degrees_fit(base, count)) {
        represent_as_aux_or_fail(node.clone());
        return;
    }
    Poly product(ext_vars.size());
    product.add_term(Monomial(ext_vars.size(), 0), Rational(1));
    for (int index = 0; index < count; ++index) {
        extend_variables(product, base.num_vars);
        product = mul_poly(product, base);
    }
    result = product;
}
}
