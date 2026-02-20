#pragma once
#include "symbolic.hpp"
#include <unordered_map>
#include <string>
#include <vector>
#include <functional>

namespace lamina {


using MatchMap = std::unordered_map<std::string, SymbolicExpr>;

class LAMINA_API Matcher {
public:
    
    
    static bool match(const SymbolicExpr& pattern, const SymbolicExpr& target, 
                      const std::unordered_set<std::string>& wildcards,
                      MatchMap& results);

    
    static SymbolicExpr replace(const SymbolicExpr& template_expr, const MatchMap& bindings);
};

struct LAMINA_API Rule {
    SymbolicExpr pattern;
    SymbolicExpr replacement;
    std::unordered_set<std::string> wildcards; 

    
    std::function<bool(const MatchMap&)> condition; 

    Rule(SymbolicExpr p, SymbolicExpr r, std::unordered_set<std::string> w, 
         std::function<bool(const MatchMap&)> c = nullptr)
        : pattern(p), replacement(r), wildcards(w), condition(c) {}
};

class LAMINA_API RewriteEngine {
    std::vector<Rule> rules;

public:
    void add_rule(const Rule& rule);
    
    
    SymbolicExpr apply(const SymbolicExpr& expr, int max_iterations = 100);

    
    SymbolicExpr apply_step(const SymbolicExpr& expr);

    const std::vector<Rule>& get_rules() const { return rules; }
};




SymbolicExpr wildcard(const std::string& name);

} 
