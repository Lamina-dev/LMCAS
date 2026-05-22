#include <variant>
#include <iostream>
int main() {
    std::variant<std::nullptr_t, int> v1 = nullptr;
    std::variant<std::nullptr_t, int> v2 = nullptr;
    std::cout << (v1 < v2) << std::endl;
}
