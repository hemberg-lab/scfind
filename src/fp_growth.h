#pragma once

#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
// #include <deque>
#include <utility>

#include "typedef.h"

using TransformedPrefixPath = std::pair<Transaction, uint64_t>;



struct FPNode {
    const Item item;
    uint64_t frequency;
    std::shared_ptr<FPNode> node_link;
    std::weak_ptr<FPNode> parent;
    std::vector<std::shared_ptr<FPNode>> children;

    FPNode(const Item&, const std::shared_ptr<FPNode>&);
};

struct FPTree {
    std::shared_ptr<FPNode> root;
    std::map<Item, std::shared_ptr<FPNode>> header_table;
    uint64_t minimum_support_threshold;

    FPTree(const std::vector<Transaction>&, uint64_t);

    bool empty() const;
};


std::set<Pattern> fptree_growth(const FPTree&);


