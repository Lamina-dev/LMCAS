#include <iostream>
#include <chrono>
#include <vector>
#include <functional>
#include <string>
#include <iomanip>

#include "../symbolic.hpp"
#include "../cas.hpp"

// Simple Timer Class
class Timer {
    std::string name;
    std::chrono::high_resolution_clock::time_point start;
public:
    Timer(const std::string& n) : name(n), start(std::chrono::high_resolution_clock::now()) {}
    ~Timer() {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        std::cout << "[BENCHMARK] " << std::left << std::setw(40) << name 
                  << ": " << std::fixed << std::setprecision(2) << elapsed.count() << " ms" << std::endl;
    }
};

void bench_polynomial_expand() {
    Timer t("Polynomial Expand (x+y+z+1)^8");
    auto x = SymbolicExpr::variable("x");
    auto y = SymbolicExpr::variable("y");
    auto z = SymbolicExpr::variable("z");
    auto one = SymbolicExpr::number(1);
    
    // (x+y+z+1)
    auto base = SymbolicExpr::add(x, SymbolicExpr::add(y, SymbolicExpr::add(z, one)));
    // ^8
    auto expr = SymbolicExpr::power(base, SymbolicExpr::number(8));
    
    auto expanded = expr->expand();
    // Do not print result to avoid I/O bottleneck
}

void bench_matrix_determinant(int n) {
    std::string name = "Matrix Determinant " + std::to_string(n) + "x" + std::to_string(n);
    Timer t(name);
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> mat_data(n, std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for(int i=0; i<n; ++i) {
        for(int j=0; j<n; ++j) {
            // M_ij = x_ij
            auto val = SymbolicExpr::variable("x_" + std::to_string(i) + "_" + std::to_string(j));
            mat_data[i][j] = val;
        }
    }
    auto mat = SymbolicExpr::matrix(mat_data);
    auto det = SymbolicExpr::determinant(mat);
    // Expand to force computation if lazy properties exist
    det = det->expand();
}

void bench_bigint_factorial() {
    Timer t("BigInt Factorial 1000!");
    BigInt res = 1;
    for(int i=1; i<=1000; ++i) {
        res = res * BigInt(i);
    }
}

void bench_symbolic_diff() {
    Timer t("Symbolic Diff (sin(x^2 + 1) * e^x)^5");
    auto x = SymbolicExpr::variable("x");
    
    // inner = sin(x^2+1) * e^x
    auto x2_1 = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(1));
    auto term1 = SymbolicExpr::sin(x2_1);
    
    // Use e^x for exp
    auto e = SymbolicExpr::variable("e");
    auto term2 = SymbolicExpr::power(e, x);
    
    auto inner = SymbolicExpr::multiply(term1, term2);
    
    auto f = SymbolicExpr::power(inner, SymbolicExpr::number(5));
    
    auto d = f->differentiate("x")->simplify();
}

void bench_linear_system() {
    // Solve 3x3 Linear System
    // x + y + z = a
    // 2x - y + z = b
    // x + 2y - z = c
    Timer t("Linear System Solve 3x3");
    
    auto x = "x"; auto y = "y"; auto z = "z";
    auto a = SymbolicExpr::variable("a");
    auto b = SymbolicExpr::variable("b");
    auto c = SymbolicExpr::variable("c");
    
    // Eq1: x + y + z - a
    auto eq1 = SymbolicExpr::add(SymbolicExpr::add(SymbolicExpr::variable(x), SymbolicExpr::variable(y)), 
               SymbolicExpr::add(SymbolicExpr::variable(z), SymbolicExpr::multiply(SymbolicExpr::number(-1), a)));
    
    // Eq2: 2x - y + z - b
    auto eq2 = SymbolicExpr::add(SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::variable(x)), 
               SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::variable(y))),
               SymbolicExpr::add(SymbolicExpr::variable(z), SymbolicExpr::multiply(SymbolicExpr::number(-1), b)));

    // Eq3: x + 2y - z - c
    auto eq3 = SymbolicExpr::add(SymbolicExpr::add(SymbolicExpr::variable(x), 
               SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::variable(y))),
               SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::variable(z)), 
               SymbolicExpr::multiply(SymbolicExpr::number(-1), c)));

    SymbolicExpr::solve_system({eq1, eq2, eq3}, {x, y, z});
}

int main() {
    std::cout << "Starting LMCAS Benchmarks..." << std::endl;
    std::cout << "========================================" << std::endl;

    bench_bigint_factorial();
    bench_symbolic_diff();
    bench_linear_system();
    bench_polynomial_expand(); 
    
    bench_matrix_determinant(3);
    bench_matrix_determinant(4);
    // bench_matrix_determinant(5); // 5x5 might be too slow for symbolic determinant without optimization

    std::cout << "========================================" << std::endl;
    std::cout << "Benchmarks completed." << std::endl;
    return 0;
}
