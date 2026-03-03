#pragma once
#include <memory>
#include <vector>
#include <string>
#include <variant>
#include <algorithm> 
#include <typeindex> 
#include <unordered_map>
#include <unordered_set>
#include "rational.hpp"
#include "bigint.hpp"



class SymbolicVisitor;
class NumberNode;
class VariableNode;
class AddNode;
class MultiplyNode;
class PowerNode;
class FunctionNode;
class MatrixNode;
class RelationalNode;


inline void hash_combine(std::size_t& seed, std::size_t value) {
    seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}


class SymbolicNode {
protected:
    mutable std::size_t cached_hash = 0;
    mutable bool hash_computed = false;

    
    virtual std::size_t compute_hash() const = 0;

    
    
    virtual int compare_same_type(const SymbolicNode& other) const = 0;

public:
    virtual ~SymbolicNode() = default;
    virtual void accept(SymbolicVisitor& visitor) = 0;
    
    virtual std::shared_ptr<SymbolicNode> clone() const = 0;
    
    
    
    virtual int type_priority() const = 0;

    std::size_t hash() const {
        if (!hash_computed) {
            cached_hash = compute_hash();
            hash_computed = true;
        }
        return cached_hash;
    }

    
    int compare(const SymbolicNode& other) const {
        if (type_priority() != other.type_priority()) {
            return type_priority() < other.type_priority() ? -1 : 1;
        }
        return compare_same_type(other);
    }
    
    
    bool equals(const SymbolicNode& other) const {
        if (this == &other) return true;
        if (hash() != other.hash()) return false; 
        if (type_priority() != other.type_priority()) return false;
        return compare_same_type(other) == 0;
    }

    virtual bool is_number() const { return false; }
    virtual bool is_one() const { return false; }
    virtual bool is_zero() const { return false; }
};


class SymbolicVisitor {
public:
    virtual ~SymbolicVisitor() = default;
    
    virtual void visit(NumberNode& node) = 0;
    virtual void visit(VariableNode& node) = 0;
    virtual void visit(AddNode& node) = 0;
    virtual void visit(MultiplyNode& node) = 0;
    virtual void visit(PowerNode& node) = 0;
    virtual void visit(FunctionNode& node) = 0;
    virtual void visit(MatrixNode& node) = 0;
    virtual void visit(RelationalNode& node) {}
};


struct NodeHash {
    std::size_t operator()(const std::shared_ptr<SymbolicNode>& node) const {
        return node ? node->hash() : 0;
    }
};

struct NodeEqual {
    bool operator()(const std::shared_ptr<SymbolicNode>& lhs, const std::shared_ptr<SymbolicNode>& rhs) const {
        if (!lhs || !rhs) return lhs == rhs;
        return lhs->equals(*rhs);
    }
};


template<typename T>
using NodeMap = std::unordered_map<std::shared_ptr<SymbolicNode>, T, NodeHash, NodeEqual>;

using NodeSet = std::unordered_set<std::shared_ptr<SymbolicNode>, NodeHash, NodeEqual>;





class SymbolicFactory {
public:
    static std::shared_ptr<SymbolicNode> create_number(const ::BigInt& v);
    static std::shared_ptr<SymbolicNode> create_number(const ::Rational& v);
    static std::shared_ptr<SymbolicNode> create_number(double v);
    static std::shared_ptr<SymbolicNode> create_variable(const std::string& name);

    static std::shared_ptr<SymbolicNode> create_add(std::vector<std::shared_ptr<SymbolicNode>> ops);

    static std::shared_ptr<SymbolicNode> create_multiply(std::vector<std::shared_ptr<SymbolicNode>> ops);
    static std::shared_ptr<SymbolicNode> create_power(std::shared_ptr<SymbolicNode> base, std::shared_ptr<SymbolicNode> exponent);
};

class NumberNode : public SymbolicNode {
public:
    std::variant<BigInt, Rational, double> value;
    
    explicit NumberNode(const BigInt& v) : value(v) {}
    explicit NumberNode(const Rational& v) : value(v) {}
    explicit NumberNode(double v) : value(v) {}
    explicit NumberNode(std::variant<BigInt, Rational, double> v) : value(std::move(v)) {}
    
    
    int type_priority() const override { return -10; }

protected:
    std::size_t compute_hash() const override {
        
        
        
        
        
        if (std::holds_alternative<Rational>(value)) {
            const auto& r = std::get<Rational>(value);
            return std::hash<std::string>{}(r.to_string()); 
        } else if (std::holds_alternative<BigInt>(value)) {
            const auto& b = std::get<BigInt>(value);
            return std::hash<std::string>{}(b.to_string());
        } else {
             auto d = std::get<double>(value);
             
             if (d == std::floor(d)) {
                 
             }
             return std::hash<double>{}(d);
        }
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const NumberNode&>(other);
        
        
        auto to_double = [](const std::variant<BigInt, Rational, double>& v) {
            if (std::holds_alternative<double>(v)) return std::get<double>(v);
            if (std::holds_alternative<Rational>(v)) return std::get<Rational>(v).to_double();
            return std::get<BigInt>(v).to_double();
        };

        double v1 = to_double(value);
        double v2 = to_double(o.value);
        
        /* 
           TODO: Proper exact comparison logic:
           1. If both BigInt -> compare
           2. If both Rational -> compare
           3. If mixed, upgrade to Rational
           4. Double fallback only for Double variants
        */

        if (v1 < v2) return -1;
        if (v1 > v2) return 1;
        return 0;
    }

public:
    void accept(SymbolicVisitor& visitor) override { visitor.visit(*this); }
    std::shared_ptr<SymbolicNode> clone() const override {
        
        if (std::holds_alternative<BigInt>(value)) return std::make_shared<NumberNode>(std::get<BigInt>(value));
        if (std::holds_alternative<Rational>(value)) return std::make_shared<NumberNode>(std::get<Rational>(value));
        return std::make_shared<NumberNode>(std::get<double>(value));
    }
    
    bool is_number() const override { return true; }
    bool is_zero() const override {
        
        if (std::holds_alternative<BigInt>(value)) return std::get<BigInt>(value) == BigInt(0);
        if (std::holds_alternative<Rational>(value)) return std::get<Rational>(value) == Rational(0);
        return std::get<double>(value) == 0.0;
    }
    bool is_one() const override {
        if (std::holds_alternative<BigInt>(value)) return std::get<BigInt>(value) == BigInt(1);
        if (std::holds_alternative<Rational>(value)) return std::get<Rational>(value) == Rational(1);
        return std::get<double>(value) == 1.0;
    }
};

class VariableNode : public SymbolicNode {
public:
    std::string name;
    
    explicit VariableNode(std::string n) : name(std::move(n)) {}

    int type_priority() const override { return 10; }

protected:
    std::size_t compute_hash() const override {
        return std::hash<std::string>{}(name);
    }
    
    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const VariableNode&>(other);
        return name.compare(o.name);
    }
    
public:
    void accept(SymbolicVisitor& visitor) override { visitor.visit(*this); }
    std::shared_ptr<SymbolicNode> clone() const override { return std::make_shared<VariableNode>(name); }
};

class AddNode : public SymbolicNode {
public:
    std::vector<std::shared_ptr<SymbolicNode>> operands;
    
    explicit AddNode(std::vector<std::shared_ptr<SymbolicNode>> ops) {
        
        for (const auto& op : ops) {
            if (auto add = std::dynamic_pointer_cast<AddNode>(op)) {
                operands.insert(operands.end(), add->operands.begin(), add->operands.end());
            } else {
                operands.push_back(op);
            }
        }
        
        std::sort(operands.begin(), operands.end(), [](const auto& a, const auto& b) {
            bool a_num = a->is_number();
            bool b_num = b->is_number();
            if (a_num != b_num) return !a_num; // Numbers last
            return a->compare(*b) < 0;
        });
    }

    int type_priority() const override { return 5; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        for (const auto& op : operands) {
            hash_combine(seed, op->hash());
        }
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const AddNode&>(other);
        if (operands.size() != o.operands.size()) {
            return operands.size() < o.operands.size() ? -1 : 1;
        }
        for (size_t i = 0; i < operands.size(); ++i) {
            int cmp = operands[i]->compare(*o.operands[i]);
            if (cmp != 0) return cmp;
        }
        return 0;
    }
    
public:
    void accept(SymbolicVisitor& visitor) override { visitor.visit(*this); }
    std::shared_ptr<SymbolicNode> clone() const override {
        std::vector<std::shared_ptr<SymbolicNode>> new_ops;
        new_ops.reserve(operands.size());
        for (const auto& op : operands) new_ops.push_back(op->clone());
        return std::make_shared<AddNode>(std::move(new_ops));
    }
};

class MultiplyNode : public SymbolicNode {
public:
    std::vector<std::shared_ptr<SymbolicNode>> operands;
    
    explicit MultiplyNode(std::vector<std::shared_ptr<SymbolicNode>> ops) {
        
        for (const auto& op : ops) {
            if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(op)) {
                operands.insert(operands.end(), mul->operands.begin(), mul->operands.end());
            } else {
                operands.push_back(op);
            }
        }
        
        std::sort(operands.begin(), operands.end(), [](const auto& a, const auto& b) {
            return a->compare(*b) < 0;
        });
    }

    int type_priority() const override { return 4; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        for (const auto& op : operands) {
            hash_combine(seed, op->hash());
        }
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const MultiplyNode&>(other);
        if (operands.size() != o.operands.size()) {
            return operands.size() < o.operands.size() ? -1 : 1;
        }
        for (size_t i = 0; i < operands.size(); ++i) {
            int cmp = operands[i]->compare(*o.operands[i]);
            if (cmp != 0) return cmp;
        }
        return 0;
    }
    
public:
    void accept(SymbolicVisitor& visitor) override { visitor.visit(*this); }
    std::shared_ptr<SymbolicNode> clone() const override {
        std::vector<std::shared_ptr<SymbolicNode>> new_ops;
        new_ops.reserve(operands.size());
        for (const auto& op : operands) new_ops.push_back(op->clone());
        return std::make_shared<MultiplyNode>(std::move(new_ops));
    }
};

class PowerNode : public SymbolicNode {
public:
    std::shared_ptr<SymbolicNode> base;
    std::shared_ptr<SymbolicNode> exponent;
    
    PowerNode(std::shared_ptr<SymbolicNode> b, std::shared_ptr<SymbolicNode> e) 
        : base(std::move(b)), exponent(std::move(e)) {}
    
    int type_priority() const override { return 0; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, base->hash());
        hash_combine(seed, exponent->hash());
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const PowerNode&>(other);
        
        // Reverse compare exponent for higher degree first
        int cmp_exp = exponent->compare(*o.exponent);
        if (cmp_exp != 0) return cmp_exp > 0 ? -1 : 1; 

        return base->compare(*o.base);
    }
        
public:
    void accept(SymbolicVisitor& visitor) override { visitor.visit(*this); }
    std::shared_ptr<SymbolicNode> clone() const override {
        return std::make_shared<PowerNode>(base->clone(), exponent->clone());
    }
};

class FunctionNode : public SymbolicNode {
public:
    enum class FuncType {
        Sin, Cos, Tan, Cot, Sec, Csc,
        ArcSin, ArcCos, ArcTan,
        Sinh, Cosh, Tanh,
        Ln, Log, Abs, Sqrt,
        Exp,
        Atan2,
        Calculus_Integral,
        Infinity,
        Limit
    };
    
    FuncType type;
    std::vector<std::shared_ptr<SymbolicNode>> arguments;
    
    FunctionNode(FuncType t, std::vector<std::shared_ptr<SymbolicNode>> args) 
        : type(t), arguments(std::move(args)) {}

    int type_priority() const override { return 2; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, static_cast<size_t>(type));
        for (const auto& arg : arguments) {
            hash_combine(seed, arg->hash());
        }
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const FunctionNode&>(other);
        if (type != o.type) {
            return static_cast<int>(type) < static_cast<int>(o.type) ? -1 : 1;
        }
        if (arguments.size() != o.arguments.size()) {
            return arguments.size() < o.arguments.size() ? -1 : 1;
        }
        for (size_t i = 0; i < arguments.size(); ++i) {
            int cmp = arguments[i]->compare(*o.arguments[i]);
            if (cmp != 0) return cmp;
        }
        return 0;
    }
        
public:
    void accept(SymbolicVisitor& visitor) override { visitor.visit(*this); }
    std::shared_ptr<SymbolicNode> clone() const override {
        std::vector<std::shared_ptr<SymbolicNode>> new_args;
        for (const auto& arg : arguments) new_args.push_back(arg->clone());
        return std::make_shared<FunctionNode>(type, std::move(new_args));
    }
};

#include <map>



class MatrixNode : public SymbolicNode {
public:
    using DenseStorage = std::vector<std::shared_ptr<SymbolicNode>>;
    using SparseStorage = std::map<size_t, std::shared_ptr<SymbolicNode>>;

    static std::variant<DenseStorage, SparseStorage> create_storage_from_grid(
        const std::vector<std::vector<std::shared_ptr<SymbolicNode>>>& grid, size_t total_elements);

    
    const std::variant<DenseStorage, SparseStorage> storage;
    const size_t rows;
    const size_t cols;

    
    MatrixNode(const std::vector<std::vector<std::shared_ptr<SymbolicNode>>>& grid) 
        : storage(create_storage_from_grid(grid, grid.size() * (grid.empty() ? 0 : grid[0].size()))),
          rows(grid.size()),
          cols(grid.empty() ? 0 : grid[0].size()) {}

    MatrixNode(size_t r, size_t c, DenseStorage dense) 
        : storage(std::move(dense)), rows(r), cols(c) {}
        
    MatrixNode(size_t r, size_t c, SparseStorage sparse) 
        : storage(std::move(sparse)), rows(r), cols(c) {}

    int type_priority() const override { return 6; }

protected:
    std::size_t compute_hash() const override {
        std::size_t seed = 0;
        hash_combine(seed, type_priority());
        hash_combine(seed, rows);
        hash_combine(seed, cols);

        
        
        if (std::holds_alternative<DenseStorage>(storage)) {
            const auto& dense = std::get<DenseStorage>(storage);
            for (size_t i = 0; i < dense.size(); ++i) {
                if (dense[i] && !dense[i]->is_zero()) {
                    hash_combine(seed, i); 
                    hash_combine(seed, dense[i]->hash());
                }
            }
        } else {
            const auto& sparse = std::get<SparseStorage>(storage);
            for (const auto& [idx, val] : sparse) {
                
                if (val && !val->is_zero()) {
                    hash_combine(seed, idx);
                    hash_combine(seed, val->hash());
                }
            }
        }
        return seed;
    }

    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const MatrixNode&>(other);
        if (rows != o.rows) return rows < o.rows ? -1 : 1;
        if (cols != o.cols) return cols < o.cols ? -1 : 1;

        
        for (size_t r = 0; r < rows; ++r) {
            for (size_t c = 0; c < cols; ++c) {
                auto v1 = this->get(r, c);
                auto v2 = o.get(r, c);
                
                bool z1 = !v1 || v1->is_zero();
                bool z2 = !v2 || v2->is_zero();
                
                if (z1 && z2) continue;
                if (z1) return -1; 
                if (z2) return 1;  
                
                int cmp = v1->compare(*v2);
                if (cmp != 0) return cmp;
            }
        }
        return 0;
    }

public:
    void accept(SymbolicVisitor& visitor) override { visitor.visit(*this); }
    
    std::shared_ptr<SymbolicNode> clone() const override {
        if (std::holds_alternative<DenseStorage>(storage)) {
            const auto& dense = std::get<DenseStorage>(storage);
            DenseStorage new_dense;
            new_dense.reserve(dense.size());
            for(const auto& e : dense) new_dense.push_back(e->clone());
            return std::make_shared<MatrixNode>(rows, cols, std::move(new_dense));
        } else {
            const auto& sparse = std::get<SparseStorage>(storage);
            SparseStorage new_sparse;
            for(const auto& [idx, node] : sparse) {
                new_sparse[idx] = node->clone();
            }
            return std::make_shared<MatrixNode>(rows, cols, std::move(new_sparse));
        }
    }

    std::shared_ptr<SymbolicNode> get(size_t r, size_t c) const {
        if (r >= rows || c >= cols) return nullptr; 
        
        size_t idx = r * cols + c;
        if (std::holds_alternative<DenseStorage>(storage)) {
            return std::get<DenseStorage>(storage)[idx];
        } else {
            const auto& sparse = std::get<SparseStorage>(storage);
            auto it = sparse.find(idx);
            if (it != sparse.end()) return it->second;
            return std::make_shared<NumberNode>(0.0); 
        }
    }

    bool is_sparse() const { return std::holds_alternative<SparseStorage>(storage); }

private:
   
};


inline std::variant<MatrixNode::DenseStorage, MatrixNode::SparseStorage> MatrixNode::create_storage_from_grid(
    const std::vector<std::vector<std::shared_ptr<SymbolicNode>>>& grid, size_t total_elements) {
    
    size_t non_zeros = 0;
    for (const auto& row : grid) {
        for (const auto& item : row) {
            if (item && !item->is_zero()) non_zeros++;
        }
    }

    bool use_sparse = total_elements > 0 && ((double)non_zeros / total_elements < 0.2);

    if (use_sparse) {
        SparseStorage s;
        size_t r = 0;
        for (const auto& row : grid) {
            size_t c = 0;
            for (const auto& item : row) {
                if (item && !item->is_zero()) {
                    s[r * grid[0].size() + c] = item;
                }
                c++;
            }
            r++;
        }
        return s;
    } else {
        DenseStorage d;
        d.reserve(total_elements);
        for (const auto& row : grid) {
            d.insert(d.end(), row.begin(), row.end());
        }
        return d;
    }
}





class RelationalNode : public SymbolicNode {
public:
    enum class Op { EQ, NEQ, LT, GT, LEQ, GEQ };
    
    std::shared_ptr<SymbolicNode> left;
    std::shared_ptr<SymbolicNode> right;
    Op op;

    RelationalNode(std::shared_ptr<SymbolicNode> l, std::shared_ptr<SymbolicNode> r, Op o)
        : left(std::move(l)), right(std::move(r)), op(o) {}

    void accept(SymbolicVisitor& visitor) override { visitor.visit(*this); }
    
    std::shared_ptr<SymbolicNode> clone() const override {
        return std::make_shared<RelationalNode>(left->clone(), right->clone(), op);
    }
    
    int type_priority() const override { return 100; } 

    
    std::size_t compute_hash() const override {
        std::size_t h = std::hash<int>{}((int)op);
        hash_combine(h, left->hash());
        hash_combine(h, right->hash());
        return h;
    }

    
    int compare_same_type(const SymbolicNode& other) const override {
        const auto& o = static_cast<const RelationalNode&>(other);
        if (op != o.op) return (int)op < (int)o.op ? -1 : 1;
        int cmp = left->compare(*o.left);
        if (cmp != 0) return cmp;
        return right->compare(*o.right);
    }
    
    static std::string op_to_string(Op op) {
        switch(op) {
            case Op::EQ: return "=";
            case Op::NEQ: return "!=";
            case Op::LT: return "<";
            case Op::GT: return ">";
            case Op::LEQ: return "<=";
            case Op::GEQ: return ">=";
            default: return "?";
        }
    }
};





inline std::shared_ptr<SymbolicNode> SymbolicFactory::create_number(const ::BigInt& v) { return std::make_shared<NumberNode>(v); }
inline std::shared_ptr<SymbolicNode> SymbolicFactory::create_number(const ::Rational& v) { return std::make_shared<NumberNode>(v); }
inline std::shared_ptr<SymbolicNode> SymbolicFactory::create_number(double v) { return std::make_shared<NumberNode>(v); }
inline std::shared_ptr<SymbolicNode> SymbolicFactory::create_variable(const std::string& name) { return std::make_shared<VariableNode>(name); }

inline std::shared_ptr<SymbolicNode> SymbolicFactory::create_add(std::vector<std::shared_ptr<SymbolicNode>> ops) {
    if (ops.empty()) return create_number(0.0);
    if (ops.size() == 1) return ops[0];
    
    std::vector<std::shared_ptr<SymbolicNode>> flat_ops;
    flat_ops.reserve(ops.size());
    for (const auto& op : ops) {
        if (op->is_zero()) continue;
        if (auto add = std::dynamic_pointer_cast<AddNode>(op)) {
            flat_ops.insert(flat_ops.end(), add->operands.begin(), add->operands.end());
        } else {
            flat_ops.push_back(op);
        }
    }
    
    if (flat_ops.empty()) return create_number(0.0);
    if (flat_ops.size() == 1) return flat_ops[0];
    return std::make_shared<AddNode>(std::move(flat_ops));
}

inline std::shared_ptr<SymbolicNode> SymbolicFactory::create_multiply(std::vector<std::shared_ptr<SymbolicNode>> ops) {
    if (ops.empty()) return create_number(1.0);
    
    for (const auto& op : ops) {
        if (op->is_zero()) return op; 
    }

    if (ops.size() == 1) return ops[0];
    
    std::vector<std::shared_ptr<SymbolicNode>> flat_ops;
    flat_ops.reserve(ops.size());
    for (const auto& op : ops) {
        if (op->is_one()) continue;
        if (auto mul = std::dynamic_pointer_cast<MultiplyNode>(op)) {
            flat_ops.insert(flat_ops.end(), mul->operands.begin(), mul->operands.end());
        } else {
            flat_ops.push_back(op);
        }
    }

    if (flat_ops.empty()) return create_number(1.0);
    if (flat_ops.size() == 1) return flat_ops[0];
    return std::make_shared<MultiplyNode>(std::move(flat_ops));
}

inline std::shared_ptr<SymbolicNode> SymbolicFactory::create_power(std::shared_ptr<SymbolicNode> base, std::shared_ptr<SymbolicNode> exponent) {
    if (!base || !exponent) return nullptr;
    if (exponent->is_zero()) return create_number(1.0);
    if (exponent->is_one()) return base;
    if (base->is_zero()) return create_number(0.0);
    if (base->is_one()) return create_number(1.0);
    return std::make_shared<PowerNode>(std::move(base), std::move(exponent));
}

