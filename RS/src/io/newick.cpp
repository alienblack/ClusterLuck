#include "io/newick.h"

#include <cctype>
#include <stdexcept>

namespace {

class Parser {
public:
    Parser(const std::string& input, int n, bool internal_preorder)
        : s(input), leaf_count(n), preorder(internal_preorder) {}

    Tree parse() {
        Tree t;
        t.root = parse_subtree(t);
        skip_ws();
        if (peek() == ';') ++pos;
        skip_ws();
        if (pos != s.size()) throw std::runtime_error("trailing Newick text");
        t.prepare(leaf_count);
        return t;
    }

private:
    const std::string& s;
    int leaf_count;
    bool preorder;
    size_t pos = 0;

    char peek() const {
        return pos < s.size() ? s[pos] : '\0';
    }

    void skip_ws() {
        while (pos < s.size() && std::isspace(static_cast<unsigned char>(s[pos]))) ++pos;
    }

    void skip_label_or_length() {
        skip_ws();
        while (pos < s.size()) {
            char c = s[pos];
            if (c == ',' || c == ')' || c == ';') break;
            ++pos;
        }
    }

    NodeId parse_subtree(Tree& t) {
        skip_ws();
        if (peek() == '(') {
            ++pos;
            NodeId v = preorder ? t.add_node(0) : -1;
            NodeId a = parse_subtree(t);
            skip_ws();
            if (peek() != ',') throw std::runtime_error("expected comma");
            ++pos;
            NodeId b = parse_subtree(t);
            skip_ws();
            if (peek() != ')') throw std::runtime_error("expected close paren");
            ++pos;
            if (!preorder) v = t.add_node(0);
            t.left[v] = a;
            t.right[v] = b;
            t.parent[a] = v;
            t.parent[b] = v;
            skip_label_or_length();
            return v;
        }

        if (!std::isdigit(static_cast<unsigned char>(peek()))) {
            throw std::runtime_error("expected leaf label");
        }
        int value = 0;
        while (std::isdigit(static_cast<unsigned char>(peek()))) {
            value = value * 10 + (s[pos] - '0');
            ++pos;
        }
        if (value <= 0 || value > leaf_count) {
            throw std::runtime_error("leaf label out of range");
        }
        skip_label_or_length();
        return t.add_node(value);
    }
};

}  // namespace

Tree parse_newick_tree(const std::string& text, int leaf_count) {
    return Parser(text, leaf_count, false).parse();
}

Tree parse_newick_tree_h45_order(const std::string& text, int leaf_count) {
    return Parser(text, leaf_count, true).parse();
}
