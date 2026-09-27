#include "test_inequality_support.hpp"

TEST(InequalityEndpoints, SmallQuadraticRoot) {
    auto x = SymbolicExpr::variable("x");

    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto bx = SymbolicExpr::multiply(
        SymbolicExpr::number(BigInt("10000000000000000")), x);
    auto ill_conditioned = SymbolicExpr::add(
        SymbolicExpr::add(x2, bx), SymbolicExpr::number(1));
    auto positive = InequalitySolver::solve_inequality(
        ill_conditioned, InequalityType::GreaterThan, "x");

    EXPECT_TRUE((positive.contains(0.0))) << "x^2 + 10^16*x + 1 > 0 must not turn its small negative root into zero";
    EXPECT_TRUE((positive.contains(-5e-17))) << "quadratic sign chart includes points above its small negative root";
    EXPECT_TRUE((!positive.contains(-2e-16))) << "quadratic sign chart excludes points between its two roots";
}

TEST(InequalityEndpoints, NearbyExactRoots) {
    auto x = SymbolicExpr::variable("x");

    auto x2 = SymbolicExpr::power(x, SymbolicExpr::number(2));
    auto scaled_x = SymbolicExpr::multiply(
        SymbolicExpr::number(
            Rational(BigInt(-1), BigInt("1000000000000"))),
        x);
    auto nearby_roots = SymbolicExpr::add(x2, scaled_x);
    auto negative = InequalitySolver::solve_inequality(
        nearby_roots, InequalityType::LessThan, "x");

    EXPECT_TRUE((negative.contains(5e-13))) << "x*(x-10^-12) < 0 contains the interval between both exact roots";
    EXPECT_TRUE((!negative.contains(-1e-13))) << "nearby exact roots exclude the lower exterior: " +
                                                     negative.to_string();
    EXPECT_TRUE((!negative.contains(2e-12))) << "nearby exact roots exclude the upper exterior";
}

TEST(InequalityEndpoints, CheckedFactorIsolations) {
    auto x = SymbolicExpr::variable("x");

    auto quadratic = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(-2));
    auto repeated = SymbolicExpr::power(
        SymbolicExpr::add(x, SymbolicExpr::number(-1)),
        SymbolicExpr::number(2));
    auto polynomial = SymbolicExpr::multiply(quadratic, repeated);

    for (InequalityType type : {
             InequalityType::GreaterThan, InequalityType::GreaterEqual,
             InequalityType::LessThan, InequalityType::LessEqual}) {
        auto solution = InequalitySolver::solve_inequality_checked(
            polynomial, type, "x");
        ASSERT_TRUE((solution.has_value())) << (solution ? "checked quartic sign sweep succeeds"
                                                         : "checked quartic sign sweep failed: " +
                                                               solution.error().message);
        if (!solution) {
            continue;
        }

        for (int numerator = -8; numerator <= 8; ++numerator) {
            const int value = (numerator * numerator - 32) *
                              (numerator - 4) * (numerator - 4);
            bool expected = false;
            switch (type) {
            case InequalityType::GreaterThan: {
                expected = value > 0;
                break;
            }
            case InequalityType::GreaterEqual: {
                expected = value >= 0;
                break;
            }
            case InequalityType::LessThan: {
                expected = value < 0;
                break;
            }
            case InequalityType::LessEqual: {
                expected = value <= 0;
                break;
            }
            }
            EXPECT_TRUE((solution.value().contains(numerator / 4.0) == expected)) << "quartic membership matches exact polynomial sign at " +
                                                                                         std::to_string(numerator) + "/4";
        }
    }
}

TEST(InequalityEndpoints, CheckedBoundedConjunction) {
    auto x = SymbolicExpr::variable("x");

    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, InequalityType>> bounded{
        {x, InequalityType::GreaterEqual},
        {SymbolicExpr::add(x, SymbolicExpr::number(-2)),
         InequalityType::LessThan}};
    auto bounded_result = InequalitySolver::solve_inequalities_checked(bounded, "x");
    ASSERT_TRUE((bounded_result.has_value())) << "checked conjunction of exact affine inequalities succeeds";
    if (bounded_result) {
        EXPECT_TRUE((bounded_result.value().contains(0.0) &&
                     bounded_result.value().contains(1.0) &&
                     !bounded_result.value().contains(2.0)))
            << "checked conjunction preserves closed and open endpoints";
    }
}

TEST(InequalityEndpoints, CheckedEmptyConjunction) {
    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, InequalityType>> empty;
    auto empty_result = InequalitySolver::solve_inequalities_checked(empty, "x");
    EXPECT_TRUE((empty_result && empty_result.value().is_entire_line())) << "empty checked conjunction denotes the entire real line";
}

TEST(InequalityEndpoints, CheckedSingleSurdConjunction) {
    auto x = SymbolicExpr::variable("x");
    auto x2_minus_2 = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(-2));

    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, InequalityType>> one_surd{
        {x2_minus_2, InequalityType::LessThan}};
    auto one_surd_result = InequalitySolver::solve_inequalities_checked(one_surd, "x");
    EXPECT_TRUE((one_surd_result && one_surd_result.value().contains(0.0))) << "single checked quadratic conjunction keeps its proven-order surd interval";
}

TEST(InequalityEndpoints, CheckedMixedSurdConjunction) {
    auto x = SymbolicExpr::variable("x");
    auto x2_minus_2 = SymbolicExpr::add(
        SymbolicExpr::power(x, SymbolicExpr::number(2)),
        SymbolicExpr::number(-2));

    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, InequalityType>> mixed_surd{
        {x2_minus_2, InequalityType::LessThan},
        {x, InequalityType::GreaterThan}};
    auto mixed_surd_result = InequalitySolver::solve_inequalities_checked(
        mixed_surd, "x");
    EXPECT_TRUE((mixed_surd_result &&
                 mixed_surd_result.value().contains(1.0) &&
                 !mixed_surd_result.value().contains(-1.0) &&
                 !mixed_surd_result.value().contains(2.0)))
        << (mixed_surd_result
                ? "mixed result=" + mixed_surd_result.value().to_string() +
                      " contains(1)=" +
                      std::to_string(mixed_surd_result.value().contains(1.0)) +
                      " contains(-1)=" +
                      std::to_string(mixed_surd_result.value().contains(-1.0)) +
                      " contains(2)=" +
                      std::to_string(mixed_surd_result.value().contains(2.0))
                : "checked conjunction failed: " +
                      mixed_surd_result.error().message);
}

TEST(InequalityEndpoints, CheckedInvalidConjunction) {
    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, InequalityType>> invalid{
        {nullptr, InequalityType::GreaterThan}};
    auto invalid_result = InequalitySolver::solve_inequalities_checked(invalid, "x");
    EXPECT_TRUE((!invalid_result && invalid_result.error().code == CasErrc::InvalidArgument)) << "checked conjunction propagates invalid component errors";
}

TEST(InequalityEndpoints, CheckedCancelledConjunction) {
    auto x = SymbolicExpr::variable("x");
    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, InequalityType>> bounded{
        {x, InequalityType::GreaterEqual},
        {SymbolicExpr::add(x, SymbolicExpr::number(-2)),
         InequalityType::LessThan}};

    CancellationToken cancellation;
    cancellation.cancel();
    ComputationContext cancelled_context(ResourceLimits{}, cancellation);
    auto cancelled = InequalitySolver::solve_inequalities_checked(
        bounded, "x", cancelled_context);
    EXPECT_TRUE((!cancelled && cancelled.error().code == CasErrc::Cancelled)) << "checked conjunction observes cancellation before processing components";
}

TEST(InequalityEndpoints, ConjunctionPreservesIntersectionErrorPayload) {
    auto x = SymbolicExpr::variable("x");
    std::vector<std::pair<std::shared_ptr<SymbolicExpr>, InequalityType>>
        bounded{
            {x, InequalityType::GreaterEqual},
            {SymbolicExpr::add(x, SymbolicExpr::number(-2)),
             InequalityType::LessThan}};

    ComputationContext measurement;
    ASSERT_TRUE(measurement.consume_steps(
        1, "solve_inequalities_checked"));
    ASSERT_TRUE(InequalitySolver::solve_inequality_checked(
        bounded[0].first, bounded[0].second, "x", measurement));
    ASSERT_TRUE(InequalitySolver::solve_inequality_checked(
        bounded[1].first, bounded[1].second, "x", measurement));

    ResourceLimits limits;
    limits.max_steps = measurement.steps_used();
    ComputationContext exhausted_at_intersection(limits);
    auto result = InequalitySolver::solve_inequalities_checked(
        bounded, "x", exhausted_at_intersection);

    ASSERT_FALSE(result);
    EXPECT_EQ(result.error().code, CasErrc::ResourceLimit);
    EXPECT_EQ(
        result.error().message, "computation step budget exhausted");
    EXPECT_EQ(result.error().operation, "interval_union_intersect");
}

std::shared_ptr<SymbolicExpr> polynomial_from_roots(const std::vector<int> &roots) {
    auto x = SymbolicExpr::variable("x");
    auto poly = SymbolicExpr::add(x, SymbolicExpr::number(-roots[0]));
    for (size_t i = 1; i < roots.size(); ++i) {
        auto factor = SymbolicExpr::add(x, SymbolicExpr::number(-roots[i]));
        poly = SymbolicExpr::multiply(poly, factor);
    }
    return poly->expand();
}

bool check_root_membership(
    const std::vector<int> &roots,
    const std::shared_ptr<SymbolicExpr> &poly,
    InequalityType type,
    const IntervalUnion &solution,
    int iter) {
    bool is_strict = type == InequalityType::GreaterThan ||
                     type == InequalityType::LessThan;
    for (int root : roots) {
        bool included = solution.contains(static_cast<double>(root));
        if (included != is_strict)
            continue;
        std::ostringstream message;
        message << "FAIL (" << (is_strict ? "strict" : "non-strict")
                << "): iter=" << iter << " root=" << root
                << " included=" << included
                << " poly=" << poly->to_string()
                << " type=" << static_cast<int>(type);
        ADD_FAILURE() << message.str();
        return false;
    }
    return true;
}

TEST(InequalityEndpoints, EndpointCorrectness) {
    std::mt19937 rng(777);
    const int num_iterations = 100;
    std::uniform_int_distribution<int> root_count_dist(1, 2);
    std::uniform_int_distribution<int> root_val_dist(-10, 10);
    std::uniform_int_distribution<int> type_dist(0, 3);
    for (int iter = 0; iter < num_iterations; ++iter) {
        SCOPED_TRACE(::testing::Message() << "iter=" << iter);
        int num_roots = root_count_dist(rng);
        std::set<int> root_set;
        while (static_cast<int>(root_set.size()) < num_roots) {
            root_set.insert(root_val_dist(rng));
        }
        std::vector<int> roots(root_set.begin(), root_set.end());
        auto poly = polynomial_from_roots(roots);
        InequalityType type = inequality_types[type_dist(rng)];
        auto solution = InequalitySolver::solve_inequality(poly, type, "x");
        EXPECT_TRUE(check_root_membership(roots, poly, type, solution, iter));
    }
}

std::shared_ptr<SymbolicExpr> repeated_root_polynomial(int root, int multiplicity) {
    auto x = SymbolicExpr::variable("x");
    auto factor = SymbolicExpr::add(x, SymbolicExpr::number(-root));
    auto poly = factor;
    for (int i = 1; i < multiplicity; ++i) {
        poly = SymbolicExpr::multiply(poly, factor);
    }
    return poly->expand();
}

int numeric_sign(double value) {
    const double sign_tol = 1e-10;
    if (value > sign_tol)
        return 1;
    if (value < -sign_tol)
        return -1;
    return 0;
}

bool check_multiplicity_signs(int root, int multiplicity, bool use_even, int iter) {
    auto poly = repeated_root_polynomial(root, multiplicity);
    const double epsilon = 0.5;
    auto sub_left = poly->substitute(
        "x", SymbolicExpr::number(static_cast<double>(root) - epsilon));
    double left_val = sub_left->simplify()->to_numeric();
    auto sub_right = poly->substitute(
        "x", SymbolicExpr::number(static_cast<double>(root) + epsilon));
    double right_val = sub_right->simplify()->to_numeric();
    int left_sign = numeric_sign(left_val);
    int right_sign = numeric_sign(right_val);
    if (left_sign == 0 || right_sign == 0)
        return true;
    if ((left_sign == right_sign) == use_even)
        return true;
    std::ostringstream message;
    message << "FAIL (" << (use_even ? "even" : "odd")
            << " mult): iter=" << iter
            << " root=" << root << " mult=" << multiplicity
            << " left_sign=" << left_sign << " right_sign=" << right_sign
            << " left_val=" << left_val << " right_val=" << right_val;
    ADD_FAILURE() << message.str();
    return false;
}

TEST(InequalityEndpoints, MultiplicitySignChanges) {
    std::mt19937 rng(888);
    const int num_iterations = 100;
    std::uniform_int_distribution<int> root_val_dist(-5, 5);
    std::uniform_int_distribution<int> even_mult_dist(0, 1);
    std::uniform_int_distribution<int> odd_mult_dist(0, 1);
    std::uniform_int_distribution<int> parity_dist(0, 1);
    for (int iter = 0; iter < num_iterations; ++iter) {
        SCOPED_TRACE(::testing::Message() << "iter=" << iter);
        int root = root_val_dist(rng);
        bool use_even = parity_dist(rng) == 0;
        int multiplicity;
        if (use_even) {
            multiplicity = even_mult_dist(rng) == 0 ? 2 : 4;
        } else {
            multiplicity = odd_mult_dist(rng) == 0 ? 1 : 3;
        }
        EXPECT_TRUE(check_multiplicity_signs(root, multiplicity, use_even, iter));
    }
}
