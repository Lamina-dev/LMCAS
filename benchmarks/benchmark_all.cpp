#include <iostream>
#include <chrono>
#include <vector>
#include <functional>
#include <string>
#include <iomanip>

#include "symbolic.hpp"



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
    
    
    auto base = SymbolicExpr::add(x, SymbolicExpr::add(y, SymbolicExpr::add(z, one)));
    
    auto expr = SymbolicExpr::power(base, SymbolicExpr::number(8));
    
    auto expanded = expr->expand();
    
}

void bench_matrix_determinant(int n) {
    std::string name = "Matrix Determinant " + std::to_string(n) + "x" + std::to_string(n);
    Timer t(name);
    
    std::vector<std::vector<std::shared_ptr<SymbolicExpr>>> mat_data(n, std::vector<std::shared_ptr<SymbolicExpr>>(n));
    for(int i=0; i<n; ++i) {
        for(int j=0; j<n; ++j) {
            
            auto val = SymbolicExpr::variable("x_" + std::to_string(i) + "_" + std::to_string(j));
            mat_data[i][j] = val;
        }
    }
    auto mat = SymbolicExpr::matrix(mat_data);
    auto det = SymbolicExpr::determinant(mat);
    
    det = det->expand();
}

void bench_bigint_factorial() {
    Timer t("BigInt Factorial 1000!");
    BigInt res = BigInt::factorial(1000);
}

void bench_bigint_combinatorics() {
    {
        Timer t("BigInt nCr(5000, 2500)");
        BigInt res = BigInt::nCr(5000, 2500);
    }
    {
        Timer t("BigInt Multinomial(1000, {250, 250, 250, 250})");
        std::vector<unsigned int> ks = {250, 250, 250, 250};
        BigInt res = BigInt::multinomial(1000, ks);
    }
}

void bench_bigint_power() {
    {
        Timer t("BigInt 3^10000");
        BigInt base(3);
        BigInt res = base.power(10000);
    }
}

void bench_symbolic_diff() {
    Timer t("Symbolic Diff (sin(x^2 + 1) * e^x)^5");
    auto x = SymbolicExpr::variable("x");
    
    
    auto x2_1 = SymbolicExpr::add(SymbolicExpr::power(x, SymbolicExpr::number(2)), SymbolicExpr::number(1));
    auto term1 = SymbolicExpr::sin(x2_1);
    
    
    auto e = SymbolicExpr::variable("e");
    auto term2 = SymbolicExpr::power(e, x);
    
    auto inner = SymbolicExpr::multiply(term1, term2);
    
    auto f = SymbolicExpr::power(inner, SymbolicExpr::number(5));
    
    auto d = f->differentiate("x")->simplify();
}

void bench_linear_system() {
    
    
    
    
    Timer t("Linear System Solve 3x3");
    
    auto x = "x"; auto y = "y"; auto z = "z";
    auto a = SymbolicExpr::variable("a");
    auto b = SymbolicExpr::variable("b");
    auto c = SymbolicExpr::variable("c");
    
    
    auto eq1 = SymbolicExpr::add(SymbolicExpr::add(SymbolicExpr::variable(x), SymbolicExpr::variable(y)), 
               SymbolicExpr::add(SymbolicExpr::variable(z), SymbolicExpr::multiply(SymbolicExpr::number(-1), a)));
    
    
    auto eq2 = SymbolicExpr::add(SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::variable(x)), 
               SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::variable(y))),
               SymbolicExpr::add(SymbolicExpr::variable(z), SymbolicExpr::multiply(SymbolicExpr::number(-1), b)));

    
    auto eq3 = SymbolicExpr::add(SymbolicExpr::add(SymbolicExpr::variable(x), 
               SymbolicExpr::multiply(SymbolicExpr::number(2), SymbolicExpr::variable(y))),
               SymbolicExpr::add(SymbolicExpr::multiply(SymbolicExpr::number(-1), SymbolicExpr::variable(z)), 
               SymbolicExpr::multiply(SymbolicExpr::number(-1), c)));

    SymbolicExpr::solve_system({eq1, eq2, eq3}, {x, y, z});
}

void bench_bigint_gcd() {
    Timer t("BigInt GCD (Fibonacci 2000)");
    BigInt a = 0;
    BigInt b = 1;
    for(int i=0; i<2000; ++i) {
        BigInt temp = b;
        b = a + b;
        a = temp;
    }
    
    for(int i=0; i<100; ++i) {
         BigInt::gcd(a, b);
    }
}

int main() {
    std::cout << "Starting LMCAS Benchmarks..." << std::endl;
    std::cout << "========================================" << std::endl;

    bench_bigint_factorial();
    bench_bigint_combinatorics();
    bench_bigint_power();
    bench_bigint_gcd();
    bench_symbolic_diff();
    bench_linear_system();
    bench_polynomial_expand(); 
    
    bench_matrix_determinant(3);
    bench_matrix_determinant(4);
    

    std::cout << "========================================" << std::endl;
    std::cout << "Benchmarks completed." << std::endl;
    return 0;
}
