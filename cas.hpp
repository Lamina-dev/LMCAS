#pragma once
#include "value.hpp"
#include <vector>

LAMINA_API Value cas_parse(const std::vector<Value>& args);
LAMINA_API Value cas_simplify(const std::vector<Value>& args);
LAMINA_API Value cas_differentiate(const std::vector<Value>& args);
LAMINA_API Value cas_integrate(const std::vector<Value>& args);
LAMINA_API Value cas_limit(const std::vector<Value>& args);
LAMINA_API Value cas_solve(const std::vector<Value>& args);
LAMINA_API Value cas_evaluate(const std::vector<Value>& args);
LAMINA_API Value cas_store(const std::vector<Value>& args);
LAMINA_API Value cas_load(const std::vector<Value>& args);
