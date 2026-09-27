#include "test_newton_raphson_support.hpp"

static bool sturm_intervals_contain_known_roots(
    const std::vector<std::pair<Rational, Rational>> &intervals,
    const std::set<int> &distinct_roots_set, int trial) {
    bool intervals_valid = true;
    for (const auto &interval : intervals) {
        double lo = interval.first.to_double();
        double hi = interval.second.to_double();

        bool contains_known_root = false;
        for (int r : distinct_roots_set) {
            double rd = (double)r;
            if (rd >= lo - 1e-10 && rd <= hi + 1e-10) {
                contains_known_root = true;
                break;
            }
        }

        if (!contains_known_root) {
            intervals_valid = false;
            std::ostringstream msg;
            msg << "Trial " << trial
                << ": interval [" << lo << "," << hi
                << "] does not contain any known root. Known roots: [";
            bool first = true;
            for (int r : distinct_roots_set) {
                if (!first)
                    msg << ",";
                msg << r;
                first = false;
            }
            msg << "]";
            ADD_FAILURE() << msg.str();
            break;
        }
    }

    return intervals_valid;
}

static bool closed_form_count_agrees_with_sturm(
    const LMCAS::Polynomial<Rational> &poly, int deg, int sturm_count, int trial) {
    bool intervals_valid = true;
    std::vector<std::shared_ptr<SymbolicExpr>> symbolic_roots;
    if (deg == 2) {
        symbolic_roots = LMCAS::solve_cubic(
            num_expr(0),
            SymbolicExpr::number(poly.coeffs[2].to_double()),
            SymbolicExpr::number(poly.coeffs[1].to_double()),
            SymbolicExpr::number(poly.coeffs[0].to_double()),
            "x");

    } else if (deg == 3) {
        symbolic_roots = LMCAS::solve_cubic(
            SymbolicExpr::number(poly.coeffs[3].to_double()),
            SymbolicExpr::number(poly.coeffs[2].to_double()),
            SymbolicExpr::number(poly.coeffs[1].to_double()),
            SymbolicExpr::number(poly.coeffs[0].to_double()),
            "x");
    } else if (deg == 4) {
        symbolic_roots = LMCAS::solve_quartic(
            SymbolicExpr::number(poly.coeffs[4].to_double()),
            SymbolicExpr::number(poly.coeffs[3].to_double()),
            SymbolicExpr::number(poly.coeffs[2].to_double()),
            SymbolicExpr::number(poly.coeffs[1].to_double()),
            SymbolicExpr::number(poly.coeffs[0].to_double()),
            "x");
    }

    std::set<double> closed_form_real_roots;
    for (const auto &root : symbolic_roots) {
        double val = test_numeric_value(root);
        if (!std::isnan(val) && !std::isinf(val)) {

            double rounded = std::round(val * 1e6) / 1e6;
            closed_form_real_roots.insert(rounded);
        }
    }

    int closed_form_count = (int)closed_form_real_roots.size();

    if (closed_form_count > sturm_count) {
        std::ostringstream msg;
        msg << "Trial " << trial
            << ": closed-form found " << closed_form_count
            << " real roots but Sturm only found " << sturm_count
            << " (degree " << deg << ")";
        ADD_FAILURE() << msg.str();
        intervals_valid = false;
    }
    return intervals_valid;
}

TEST(NewtonRaphsonIsolation, SturmIsolationX22Has2RealRootsEarly) {
    LMCAS::Polynomial<Rational> poly("x");
    poly.coeffs = {Rational(-2), Rational(0), Rational(1)};

    auto intervals = LMCAS::isolate_real_roots_checked(poly).value();
    EXPECT_TRUE((intervals.size() == 2)) << "x^2-2 should have 2 isolated real roots";
}

TEST(NewtonRaphsonIsolation, SturmSequenceRootCountAccuracy) {
    const int NUM_TRIALS = 60;

    std::mt19937 rng(314159);
    std::uniform_int_distribution<int> degree_dist(2, 5);
    std::uniform_int_distribution<int> root_dist(-4, 4);

    for (int trial = 0; trial < NUM_TRIALS; ++trial) {
        int deg = degree_dist(rng);

        std::vector<int> all_roots;
        for (int i = 0; i < deg; ++i) {
            all_roots.push_back(root_dist(rng));
        }

        std::set<int> distinct_roots_set(all_roots.begin(), all_roots.end());
        int expected_distinct_real_roots = (int)distinct_roots_set.size();

        LMCAS::Polynomial<Rational> poly = poly_from_roots(all_roots);

        auto intervals = LMCAS::isolate_real_roots_checked(poly).value();
        int sturm_count = (int)intervals.size();

        bool count_matches = (sturm_count == expected_distinct_real_roots);

        if (!count_matches) {
            std::ostringstream msg;
            msg << "Trial " << trial << ": Sturm found "
                << sturm_count << " roots, expected " << expected_distinct_real_roots
                << " distinct real roots (degree " << deg << ", roots: [";
            for (size_t k = 0; k < all_roots.size(); ++k) {
                if (k > 0)
                    msg << ",";
                msg << all_roots[k];
            }
            msg << "])";
            ADD_FAILURE() << msg.str();
            continue;
        }

        bool intervals_valid = sturm_intervals_contain_known_roots(
            intervals, distinct_roots_set, trial);

        if (deg <= 4 && intervals_valid) {

            intervals_valid = closed_form_count_agrees_with_sturm(
                poly, deg, sturm_count, trial);
        }

        EXPECT_TRUE(intervals_valid && count_matches)
            << "Trial " << trial << ", degree=" << deg;
    }
}

TEST(NewtonRaphsonIsolation, SturmIsolationX22Has2RealRoots) {
    LMCAS::Polynomial<Rational> poly("x");
    poly.coeffs = {Rational(-2), Rational(0), Rational(1)};

    auto intervals = LMCAS::isolate_real_roots_checked(poly).value();
    EXPECT_TRUE((intervals.size() == 2)) << "x^2-2 should have 2 isolated real roots";
}
