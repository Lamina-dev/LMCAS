#include "bigint.hpp"
#include "rational.hpp"
#include <iostream>
#include <string>
#include <cassert>

int main() {
    try {
        std::cout << "--- BigInt Debug Test ---" << std::endl;

    {
        BigInt two(2);
        std::string s = two.ToString();
        std::cout << "[PASS] BigInt(2).to_string(): " << s << std::endl;
        if (s != "2") std::cerr << "FAIL" << std::endl;
    }
    {
        BigInt x("123");
        std::string s = x.ToString();
        std::cout << "[PASS] BigInt('123').to_string(): " << s << std::endl;
        if (s != "123") std::cerr << "FAIL expect 123" << std::endl;
    }
    {
        BigInt a(2), b(3);
        BigInt c = a * b;
        std::cout << "[PASS] 2 * 3: " << c.ToString() << std::endl;
        if (c.ToString() != "6") std::cerr << "FAIL" << std::endl;
    }
    {
        std::string large = "123456789123456789";
        BigInt x(large);
        std::cout << "[PASS] BigInt(large).to_string(): " << x.ToString() << std::endl;
        if (x.ToString() != large) std::cerr << "FAIL large string" << std::endl;
        
        std::cout << "Calculating large * large..." << std::endl;
        BigInt y = x * x;
        std::cout << "Result computed." << std::endl;
        std::cout << "Result string: " << y.ToString() << std::endl;
        
        std::cout << "Testing 2-limb % 1-limb check..." << std::endl;
        BigInt one(1);
        BigInt zero = y % one;
        std::cout << "y % 1 = " << zero.ToString() << std::endl;
        
    }
    
    std::cout << "--- Rational Debug Test ---" << std::endl;
    {
       // Test modulo
       BigInt n1("123456789");
       BigInt n2("987654321");
       
       std::cout << "Testing BigInt % ..." << std::endl;
       BigInt rem = n2 % n1; 
       std::cout << "987654321 % 123456789 = " << rem.ToString() << std::endl;
       // 987654321 / 123456789 = 8.0000000729
       // 123456789 * 8 = 987654312
       // 987654321 - 987654312 = 9
       if (rem.ToString() != "9") std::cerr << "FAIL % operator" << std::endl;
       
       std::cout << "Testing Rational Construction ..." << std::endl;
       Rational r1(n1); // Denom=1
       Rational r2(n2);
       
       std::cout << "Testing Rational * ..." << std::endl;
       Rational r3 = r1 * r2;
       std::cout << "Rational mult result call to_string..." << std::endl;
       std::string s = r3.to_string();
       std::cout << "Rational mult result got string." << std::endl;
       std::cout << "Rational mult result: " << s << std::endl;
    }

    } catch (const std::exception& e) {
        std::cerr << "Caught exception: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Caught unknown exception" << std::endl;
        return 1;
    }
    return 0;
}
