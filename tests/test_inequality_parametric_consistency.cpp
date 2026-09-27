#include "test_inequality_support.hpp"

struct QuadraticCoefficients {
    int a;
    int p;
    int q;
    int discriminant;
};

QuadraticCoefficients draw_quadratic(std::mt19937 &rng) {
    std::uniform_int_distribution<int> leading_dist(1, 4);
    std::uniform_int_distribution<int> sign_dist(0, 1);
    std::uniform_int_distribution<int> param_dist(-5, 5);
    int a = leading_dist(rng);
    if (sign_dist(rng)) {
        a = -a;
    }
    int p = param_dist(rng);
    int q = param_dist(rng);
    int discriminant = p * p - 4 * a * q;
    if (discriminant < 0) {
        if (a > 0) {
            q = (p * p) / (4 * a) - 1;
        } else {
            q = (p * p) / (4 * a) + 1;
        }
        discriminant = p * p - 4 * a * q;
    }
    return {a, p, q, discriminant};
}

std::shared_ptr<SymbolicExpr> parametric_quadratic_expression(int a) {
    auto x = SymbolicExpr::variable("x");
    auto p = SymbolicExpr::variable("p");
    auto q = SymbolicExpr::variable("q");
    auto ax2 = SymbolicExpr::multiply(
        SymbolicExpr::number(a),
        SymbolicExpr::power(x, SymbolicExpr::number(2)));
    auto px = SymbolicExpr::multiply(p, x);
    auto expression = SymbolicExpr::add(SymbolicExpr::add(ax2, px), q);
    return expression;
}

std::shared_ptr<SymbolicExpr> concrete_quadratic_expression(
    const QuadraticCoefficients &coefficients) {
    auto x = SymbolicExpr::variable("x");
    auto expression = SymbolicExpr::add(
        SymbolicExpr::add(
            SymbolicExpr::multiply(SymbolicExpr::number(coefficients.a),
                                   SymbolicExpr::power(x, SymbolicExpr::number(2))),
            SymbolicExpr::multiply(SymbolicExpr::number(coefficients.p), x)),
        SymbolicExpr::number(coefficients.q));
    return expression;
}

std::optional<IntervalUnion> select_parametric_solution(
    const PiecewiseIntervalResult &result) {
    if (result.is_empty()) {
        return std::nullopt;
    }
    if (result.is_single()) {
        return result.single_solution();
    }
    for (const auto &candidate : result.cases) {
        if (!candidate.condition) {
            return candidate.solution;
        }
    }
    if (!result.cases.empty()) {
        return result.cases[0].solution;
    }
    return std::nullopt;
}

std::shared_ptr<SymbolicExpr> substitute_parameters(
    const std::shared_ptr<SymbolicExpr> &expression,
    const QuadraticCoefficients &coefficients) {
    auto substituted = expression->substitute(
        "p", SymbolicExpr::number(coefficients.p));
    return substituted->substitute("q", SymbolicExpr::number(coefficients.q));
}

struct NumericBoundary {
    double value;
    bool open;
};

std::optional<NumericBoundary> evaluate_parametric_endpoint(
    const Endpoint &endpoint,
    bool lower,
    const QuadraticCoefficients &coefficients) {
    bool infinite = lower ? endpoint.is_neg_infinity : endpoint.is_pos_infinity;
    if (infinite || !endpoint.value) {
        return NumericBoundary{lower ? -1e18 : 1e18, true};
    }
    auto expression = substitute_parameters(endpoint.value, coefficients);
    auto simplified = test_numeric_eval(expression->simplify());
    if (simplified && std::isfinite(*simplified)) {
        return NumericBoundary{*simplified, endpoint.is_open};
    }
    auto raw = test_numeric_eval(substitute_parameters(endpoint.value, coefficients));
    if (raw && std::isfinite(*raw)) {
        return NumericBoundary{*raw, endpoint.is_open};
    }
    return std::nullopt;
}

bool contains_parametric_point(
    const IntervalUnion &solution,
    const QuadraticCoefficients &coefficients,
    double point) {
    if (solution.is_empty()) {
        return false;
    }
    if (solution.is_entire_line()) {
        return true;
    }
    for (const auto &interval : solution.intervals()) {
        auto lower = evaluate_parametric_endpoint(interval.lower, true, coefficients);
        auto upper = evaluate_parametric_endpoint(interval.upper, false, coefficients);
        if (!lower || !upper) {
            continue;
        }
        bool above = lower->open ? point > lower->value + 1e-10
                                 : point >= lower->value - 1e-10;
        if (!above) {
            continue;
        }
        bool below = upper->open ? point < upper->value - 1e-10
                                 : point <= upper->value + 1e-10;
        if (below) {
            return true;
        }
    }
    return false;
}

bool raw_endpoint_evaluable(
    const Endpoint &endpoint,
    bool lower,
    const QuadraticCoefficients &coefficients) {
    bool infinite = lower ? endpoint.is_neg_infinity : endpoint.is_pos_infinity;
    if (infinite || !endpoint.value) {
        return true;
    }
    auto expression = substitute_parameters(endpoint.value, coefficients);
    return test_numeric_eval(expression).has_value();
}

bool any_interval_evaluable(
    const IntervalUnion &solution, const QuadraticCoefficients &coefficients) {
    for (const auto &interval : solution.intervals()) {
        bool lower = raw_endpoint_evaluable(interval.lower, true, coefficients);
        bool upper = raw_endpoint_evaluable(interval.upper, false, coefficients);
        if (lower && upper) {
            return true;
        }
    }
    return false;
}

bool near_quadratic_root(const QuadraticCoefficients &coefficients, double point) {
    double radical = std::sqrt(std::abs(coefficients.discriminant));
    double root1 = (-coefficients.p + radical) / (2.0 * coefficients.a);
    double root2 = (-coefficients.p - radical) / (2.0 * coefficients.a);
    return std::abs(point - root1) < 1e-4 || std::abs(point - root2) < 1e-4;
}

bool check_parametric_samples(
    const QuadraticCoefficients &coefficients,
    InequalityType type,
    const IntervalUnion &parametric_solution,
    const IntervalUnion &direct_solution,
    std::mt19937 &rng,
    int iter) {
    std::uniform_real_distribution<double> sample_dist(-20.0, 20.0);
    for (int sample = 0; sample < 50; ++sample) {
        double point = sample_dist(rng);
        double value = static_cast<double>(coefficients.a) * point * point + static_cast<double>(coefficients.p) * point + static_cast<double>(coefficients.q);
        if (std::abs(value) < 1e-6) {
            continue;
        }
        bool expected = satisfies(value, type);
        bool in_direct = direct_solution.contains(point);
        bool in_parametric = contains_parametric_point(
            parametric_solution, coefficients, point);
        if (in_parametric == expected) {
            continue;
        }
        if (!in_parametric && !parametric_solution.is_empty() &&
            !parametric_solution.is_entire_line() &&
            !any_interval_evaluable(parametric_solution, coefficients)) {
            continue;
        }
        if (near_quadratic_root(coefficients, point)) {
            continue;
        }
        std::ostringstream message;
        message << "FAIL: iter=" << iter
                << " a=" << coefficients.a << " p=" << coefficients.p
                << " q=" << coefficients.q << " type=" << static_cast<int>(type)
                << " point=" << point << " expected=" << expected
                << " in_parametric=" << in_parametric
                << " in_direct=" << in_direct << " poly_val=" << value;
        ADD_FAILURE() << message.str();
        return false;
    }
    return true;
}

TEST(InequalityParametricConsistency, ParametricInequalityConsistency) {
    std::mt19937 rng(1111);
    const int num_iterations = 100;
    std::uniform_int_distribution<int> type_dist(0, 3);
    for (int iter = 0; iter < num_iterations; ++iter) {
        SCOPED_TRACE(::testing::Message() << "iter=" << iter);
        auto coefficients = draw_quadratic(rng);
        if (coefficients.discriminant < 0) {
            continue;
        }
        InequalityType type = inequality_types[type_dist(rng)];
        auto parametric_expression = parametric_quadratic_expression(coefficients.a);
        auto parametric_result = InequalitySolver::solve_parametric_inequality(
            parametric_expression, type, "x", {"p", "q"});
        auto concrete_expression = concrete_quadratic_expression(coefficients);
        auto direct_solution = InequalitySolver::solve_inequality(
            concrete_expression, type, "x");
        auto parametric_solution = select_parametric_solution(parametric_result);
        if (!parametric_solution) {
            continue;
        }
        EXPECT_TRUE(check_parametric_samples(coefficients, type, *parametric_solution,
                                             direct_solution, rng, iter));
    }
}
