#include "RootedNetwork.hpp"

#include <boost/tokenizer.hpp>

#include <stack>
#include <queue>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <ranges>
#include <tuple>
#include <utility>
#include <format>
#include <limits>


RootedNetwork::RootedNetwork() {
    root = new Node;
    insert_node(root);
}

// I. destructor
RootedNetwork::~RootedNetwork() {
    for (auto e : edges)
        delete e;
    for (auto u : nodes)
        delete u;
}

// II. copy constructor
RootedNetwork::RootedNetwork(const RootedNetwork& other) {
    std::unordered_map<Node*, Node*> node_map;
    for (auto u : other.nodes) {
        Node* ux = new Node(*u);
        node_map[u] = ux;
        insert_node(ux);
    }
    for (auto e : other.edges) {
        Edge* ex = new Edge(*e);
        ex->parent = node_map[e->parent];
        ex->child = node_map[e->child];
        insert_edge(ex);
    }
    root = node_map[other.root];
}

// III. copy assignment
RootedNetwork& RootedNetwork::operator=(const RootedNetwork& other) {
    if (this == &other)
        return *this;

    for (auto e : edges)
        delete e;
    for (auto u : nodes)
        delete u;

    return *this = RootedNetwork(other);
}

// IV. move constructor
RootedNetwork::RootedNetwork(RootedNetwork&& other) noexcept {
    root = std::move(other.root);
    nodes = std::move(other.nodes);
    edges = std::move(other.edges);
    incoming = std::move(other.incoming);
    outgoing = std::move(other.outgoing);
}

// v move assignment
RootedNetwork& RootedNetwork::operator=(RootedNetwork&& other) noexcept {
    if (this == &other)
        return *this;
  
    for (auto e : edges)
        delete e;
    for (auto u : nodes)
        delete u;

    root = std::move(other.root);
    nodes = std::move(other.nodes);
    edges = std::move(other.edges);
    incoming = std::move(other.incoming);
    outgoing = std::move(other.outgoing);

    return *this;
}

Edge* RootedNetwork::get_edge(Node* u, Node* v) const { // get the edge u->v
    for (auto e : incoming.at(v)) {
        if (e->parent == u) {
            return e;
        }
    }
    throw std::out_of_range("No edge found.");
}

std::vector<Node*> RootedNetwork::get_leaves() const {
    std::vector<Node*> leaves;
    for (auto u : nodes) {
        if (is_leaf(u)) {
            leaves.push_back(u);
        }
    }
    return leaves;
}

std::vector<Node*> RootedNetwork::get_leaves(Node* u) const {
    std::vector<Node*> leaves;
    std::queue<Node*> q({u});
    std::unordered_set<Node*> visited;
    while (!q.empty()) {
        Node* next_node = q.front();
        q.pop();
        if (visited.contains(next_node)) {
            continue;
        }
        visited.insert(next_node);
        if (is_leaf(next_node)) {
            leaves.push_back(next_node);
        } else {
            for (auto child : get_children(next_node)) {
                q.push(child);
            }
        }
    }
    return leaves;
}

std::vector<std::string> RootedNetwork::get_leaf_labels() const {
    std::vector<std::string> leaf_labels;
    for (auto u : nodes) {
        if (is_leaf(u)) {
            leaf_labels.push_back(u->label);
        }
    }
    return leaf_labels;
}

std::vector<Node*> RootedNetwork::get_reticulation_nodes() const {
    std::vector<Node*> hybrids;
    for (auto u : nodes) {
        if (is_hybrid(u)) {
            hybrids.push_back(u);
        }
    }
    return hybrids;
}

std::vector<Node*> RootedNetwork::get_tree_nodes() const {
    std::vector<Node*> tree_nodes;
    for (auto u : nodes) {
        if (!is_hybrid(u)) {
            tree_nodes.push_back(u);
        }
    }
    return tree_nodes;
}

std::vector<Node*> RootedNetwork::get_internal_nodes() const {
    std::vector<Node*> internal_nodes;
    for (auto u : nodes) {
        if (!is_leaf(u)) {
            internal_nodes.push_back(u);
        }
    }
    return internal_nodes;
}

std::vector<Node*> RootedNetwork::get_parents(Node* u) const {
    std::vector<Node*> parents;
    for (auto e : incoming.at(u)) {
        parents.push_back(e->parent);
    }
    return parents;
}

std::vector<Node*> RootedNetwork::get_children(Node* u) const {
    std::vector<Node*> children;
    for (auto e : outgoing.at(u)) {
        children.push_back(e->child);
    }
    return children;
}

std::vector<Node*> RootedNetwork::get_siblings(Node* u) const {
    std::vector<Node*> siblings;
    for (auto p : get_parents(u)) {
        for (auto child : get_children(p)) {
            if (child != u) {
                siblings.push_back(child);
            }
        }
    }
    return siblings;
}

bool RootedNetwork::is_parent(Node* u, Node* v) const { // is u the parent of v
    for (auto e : incoming.at(v)) {
        if (e->child == u) {
            return true;
        }
    }
    return false;
}

bool RootedNetwork::is_child(Node* u, Node* v) const { // is u the child of v
    for (auto e : incoming.at(u)) {
        if (e->parent == v) {
            return true;
        }
    }
    return false;
}

bool RootedNetwork::is_descendant(Node* u, Node* v) const { // is u a descendant of v
    std::queue<Node*> desc_queue({v});
    std::unordered_set<Node*> visited;
    while (!desc_queue.empty()) {
        auto current = desc_queue.front();
        desc_queue.pop();
        if (visited.contains(current))
            continue;
        if (current == u)
            return true;
        
        for (auto child : get_children(current)) {
            desc_queue.push(child);
        }
        visited.insert(current);
    }
    return false;
}

bool RootedNetwork::adjacent(Node* u, Node* v) const {
    for (auto e : incoming.at(u)) {
        if (e->parent == v) {
            return true;
        }
    }
    for (auto e : incoming.at(v)) {
        if (e->parent == u) {
            return true;
        }
    }
    return false;
}

void RootedNetwork::insert_node(Node* u) {
    nodes.push_back(u);
    incoming[u] = std::vector<Edge*>();
    outgoing[u] = std::vector<Edge*>();
}

void RootedNetwork::delete_node(Node* u) {
    auto removed_edge_pos = edges.end();
    for (auto e : incoming[u]) {
        std::erase(outgoing[e->parent],e);
        removed_edge_pos = std::remove(edges.begin(), removed_edge_pos, e);
        delete e;
    }
    for (auto e : outgoing[u]) {
        std::erase(incoming[e->child],e);
        removed_edge_pos = std::remove(edges.begin(), removed_edge_pos, e);
        delete e;
    }
    
    incoming.erase(u);
    outgoing.erase(u);
    edges.erase(removed_edge_pos, edges.end());
    std::erase(nodes, u);
    delete u;
}


Node* RootedNetwork::insert_attachment_point(Edge* e, double ratio) {
    // TODO: check if ratio is between 0 and 1
    double total_dist = e->distance;
    double admix = e->admixture;
    Node* child = e->child;
    Node* parent = e->parent;
    delete_edge(e);
    Node* attachment_point = new Node;
    insert_node(attachment_point);
    auto parent_edge = insert_edge(parent, attachment_point);
    parent_edge->admixture = admix;
    parent_edge->distance = ratio * total_dist;
    auto child_edge = insert_edge(attachment_point, child);
    child_edge->distance = (1-ratio) * total_dist;
    return attachment_point;
}

void RootedNetwork::merge_nodes(Node* n1, Node* n2) {
    for (auto x : incoming[n2]) {
        Edge* e = new Edge(*x);
        e->child = n1;
        insert_edge(e);
    }
    for (auto x : outgoing[n2]) {
        Edge* e = new Edge(*x);
        e->parent = n1;
        insert_edge(e);
    }
    delete_node(n2);
}

void RootedNetwork::insert_edge(Edge* e) {
    edges.push_back(e);
    incoming[e->child].push_back(e);
    outgoing[e->parent].push_back(e);
}

Edge* RootedNetwork::insert_edge(Node* parent, Node* child) {
    Edge* e = new Edge;
    e->child = child;
    e->parent = parent;
    insert_edge(e);
    return e;
}

void RootedNetwork::delete_edge(Edge* e) {
    std::erase(incoming[e->child], e);
    std::erase(outgoing[e->parent],e);
    std::erase(edges, e);
    delete e;
}

// bool RootedNetwork::rNNI(Edge* arc1, Edge* arc2, Edge* arc3) {
    // returns false if unmodified
    // https://doi.org/10.1371/journal.pcbi.1005611
    // s-u-v-t --> s-v-u-t
    // Only modify if all of the following will hold:
    // (1) the in- and outdegrees of s and t are not affected by the move;
    // (2) the in- and outdegrees of u and v remain at most 2 and
    // (3) the obtained network is acyclic.

//     return false;
// }

RootedNetwork parse_network(std::string s) {
    // TODO: destroy the network if parsing fails/invalid string is produced
    RootedNetwork net;
    Node* root = net.get_root();
    std::stack<Node*> node_stack({root});
    std::stack<Edge*> edge_stack;
    std::unordered_map<std::string, std::vector<Node*>> hybrids;
    
    int next_type = 0;
    bool finished = false;
    
    typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
    boost::char_separator<char> sep{"", "(),;:", boost::drop_empty_tokens};
    tokenizer tok{s, sep};
    
    for (const auto &t : tok) {
        if (finished) {
            // TODO: throw an error
        }
        switch(t[0]) {
            case '(':
                {
                    Node *new_node = new Node;
                    net.insert_node(new_node);
                    Edge *new_edge = net.insert_edge(node_stack.top(), new_node);
                    node_stack.push(new_node);
                    edge_stack.push(new_edge);
                }
                next_type = 0;
                break;
            case ')':
                node_stack.pop();
                edge_stack.pop();

                next_type = 0;
                break;
            case ',': 
                {
                    node_stack.pop();
                    edge_stack.pop();
                    
                    Node *new_node = new Node;
                    net.insert_node(new_node);
                    Edge *new_edge = net.insert_edge(node_stack.top(), new_node);
                    node_stack.push(new_node);
                    edge_stack.push(new_edge);
                }
                next_type = 0;
                break;
            case ':': // next token is either a branch length, support value, or admixture freq
                next_type += 1;
                break;
            case ';':
                finished = true;
                break;
            default: // name, branch length, etc
                if (next_type == 0) { // label
                    node_stack.top()->label = t;
                    if (t[0] == '#') { // hybrid node -- add to hybrids
                        if (hybrids.contains(t)) {
                            hybrids[t].push_back(node_stack.top());
                        } else {
                            hybrids[t] = std::vector({node_stack.top()});
                        }
                    }
                } else if (next_type == 1) { // branch length
                    edge_stack.top()->distance = std::stod(t);
                } else if (next_type == 2) { // support
                    edge_stack.top()->support = std::stod(t);
                } else if (next_type == 3) { // admixture frequency
                    edge_stack.top()->admixture = std::stod(t);
                } else {
                    // TODO: throw err
                }
                break;
        }
    }
    if (!finished) {
        // TODO: throw err
    }
    for (auto hybrid_nodes : hybrids) {
        if (hybrid_nodes.second.size() < 2) {
            // TODO: throw err
        }
        Node* merger = hybrid_nodes.second[0];
        for (auto other = (++hybrid_nodes.second.begin()); other < hybrid_nodes.second.end(); ++other) {
            net.merge_nodes(merger, *other);
        }
    } 
    return net;
}

std::pair<std::set<std::set<std::set<std::string>>>, std::set<std::set<std::set<std::string>>>> 
    rf_helper(const RootedNetwork& net1, const RootedNetwork& net2) {
    
    auto ll1 = net1.get_leaf_labels();
    auto ll2 = net2.get_leaf_labels();
    std::set<std::string> leaves_1(ll1.begin(), ll1.end());
    std::set<std::string> leaves_2(ll2.begin(), ll2.end());
    std::set<std::string> common_leaves;
    std::set_intersection(
        leaves_1.begin(), leaves_1.end(), 
        leaves_2.begin(), leaves_2.end(), 
        std::inserter(common_leaves, common_leaves.begin()));

    std::set<std::set<std::set<std::string>>> bips1, bips2;
    for (auto u : net1.get_internal_nodes()) {
        if (net1.is_root(u)) {
            continue;
        }
        
        std::set<std::string> in_cluster, out_cluster;
        for (auto leaf : net1.get_leaves(u)) {
            in_cluster.insert(leaf->label);
        }
        std::set_difference(
            in_cluster.begin(), in_cluster.end(), 
            common_leaves.begin(), common_leaves.end(), 
            std::inserter(out_cluster, out_cluster.begin()));
        bips1.insert({in_cluster,out_cluster});
    }
    for (auto u : net2.get_internal_nodes()) {
        if (net2.is_root(u)) {
            continue;
        }
        
        std::set<std::string> in_cluster, out_cluster;
        for (auto leaf : net2.get_leaves(u)) {
            in_cluster.insert(leaf->label);
        }
        std::set_difference(
            in_cluster.begin(), in_cluster.end(), 
            common_leaves.begin(), common_leaves.end(), 
            std::inserter(out_cluster, out_cluster.begin()));
        bips2.insert({in_cluster,out_cluster});
    }
    return std::make_pair(bips1,bips2);
}

unsigned int robinson_foulds(const RootedNetwork& net1, const RootedNetwork& net2) {
    auto [bips1, bips2] = rf_helper(net1, net2);
    std::vector<std::set<std::set<std::string>>> symdif;
    std::set_symmetric_difference(bips1.begin(), bips1.end(), bips2.begin(), bips2.end(), std::back_inserter(symdif));
    return symdif.size();
}

double normalized_robinson_foulds(const RootedNetwork& net1, const RootedNetwork& net2) {
    auto [bips1, bips2] = rf_helper(net1, net2);
    std::set<std::set<std::set<std::string>>> symdif, uni;
    std::set_symmetric_difference(bips1.begin(), bips1.end(), bips2.begin(), bips2.end(), std::inserter(symdif, symdif.begin()));
    return (double)symdif.size() / (bips1.size() + bips2.size()); //
}

// topological sort.
std::vector<std::vector<Node*>> transposed_kahn_ordering(const RootedNetwork &net) {
    std::vector<std::vector<Node*>> layers({net.get_leaves()});
    int remaining_visits = net.num_nodes() - 1;
    
    std::unordered_map<Node*, int> child_counts;
    for (auto u : net.get_nodes()) {
        child_counts[u] = net.outdegree(u);
    }

    while (remaining_visits) {
        std::vector<Node*> new_layer;
        for (auto u : layers.back()) {
            for (auto parent : net.get_parents(u)) {
                if (!--child_counts[parent])
                    new_layer.push_back(parent);
            }
        }
        remaining_visits -= layers.back().size();
        layers.push_back(new_layer);
    }
    
    return layers;
}

std::vector<std::pair<Node*,int>> unique_nodes(const RootedNetwork &net, const std::vector<std::vector<Node*>> layers) {
    std::unordered_map<Node*,std::set<Node*>> equivalence; 
    std::vector<std::pair<Node*,int>> unique;
    
    // map equivalence of leaves
    for (auto u = layers[0].begin(); u != layers[0].end(); u++) {
        for (auto v = u+1; v != layers[0].end(); v++) {
            if ((*u)->label == (*v)->label) {
                equivalence[*u].insert(*v);
                equivalence[*v].insert(*u);
            }
        }
    }

    for (auto layer_it = layers.begin()+1; layer_it != layers.end(); layer_it++) {
        for (auto u = layer_it->begin(); u != layer_it->end(); u++) {
            for (auto v = u+1; v != layer_it->end(); v++) {
                bool equivalent = true;
                for (auto c1 : net.get_children(*u)) {
                    bool found_equiv = false;
                    for (auto c2 : net.get_children(*v)) {
                        if (equivalence[c1].contains(c2)) {
                            found_equiv = true;
                            break;
                        }
                    }
                    if (!found_equiv) {
                        equivalent = false;
                        break;
                    }
                }
                if (equivalent) {
                    equivalence[*u].insert(*v);
                    equivalence[*v].insert(*u);
                }
            }
        }
    }

    std::set<Node*> visited;
    for (auto x : net.get_nodes()) {
        if (!visited.contains(x)) {
            unique.push_back(std::pair(x, equivalence[x].size()+1));
            // visited.insert(x); // probably unnecessary
            visited.insert(equivalence[x].begin(), equivalence[x].end());
        }
    }
    
    return unique;
}

std::pair<
    std::unordered_map<Node*, std::set<Node*>>, 
    std::unordered_map<Node*, std::set<Node*>>
> equivalent_nodes(
    const RootedNetwork& net1, const std::vector<std::vector<Node*>>& layers1, 
    const RootedNetwork& net2, const std::vector<std::vector<Node*>>& layers2) {
    std::unordered_map<Node*, std::set<Node*>> eq12,eq21;

    for (const auto& [u,v] : std::views::cartesian_product(layers1[0], layers2[0])) {
        if (u->label == v->label) {
            eq12[u].insert(v);
            eq21[v].insert(u);
        }
    }
    size_t last_layer = std::min(layers1.size(),layers2.size());
    for (size_t ix=1; ix < last_layer; ++ix) {
        for (const auto& [u,v] : std::views::cartesian_product(layers1[ix], layers2[ix])) {
            if (net1.outdegree(u) != net2.outdegree(v)) {
                continue;
            }
            auto u_children = net1.get_children(u);
            auto v_children = net2.get_children(v);
            std::set<Node*> u_child_set(u_children.begin(), u_children.end());
            std::set<Node*> v_child_set(v_children.begin(), v_children.end());
            
            bool all_children_have_equivalents = true;
            for (const auto& u_child : net1.get_children(u)) {
                std::vector<Node*> intersect;
                std::set_intersection(
                    eq12[u_child].begin(), eq12[u_child].end(), 
                    v_child_set.begin(), v_child_set.end(), 
                    std::back_inserter(intersect)
                );
                if (intersect.empty()) {
                    all_children_have_equivalents = false;
                    break;
                }
            }
            if (all_children_have_equivalents) {
                eq12[u].insert(v);
            }

            all_children_have_equivalents = true;
            for (const auto& v_child : net2.get_children(v)) {
                std::vector<Node*> intersect;
                std::set_intersection(
                    eq21[v_child].begin(), eq21[v_child].end(), 
                    u_child_set.begin(), u_child_set.end(), 
                    std::back_inserter(intersect)
                );
                if (intersect.empty()) {
                    all_children_have_equivalents = false;
                    break;
                }
            }
            if (all_children_have_equivalents) {
                eq21[v].insert(u);
            }
        }
    }
    
    return std::pair(eq12,eq21);
}

unsigned int nakhleh_distance(const RootedNetwork &net1, const RootedNetwork &net2) {
    // https://doi.org/10.1109/TCBB.2009.2
    // See also https://doi.org/10.1109/TCBB.2009.33
    // The two Nakhleh distance implementations that I've
    // encountered are incorrect in different ways.
    //
    // In keeping with this tradition, I will introduce a third,
    // incorrect implementation by omitting the division-by-2.
    //
    // If you're not a fan of that, you can divide the result by 2
    // or better yet:
    // implement your own (hopefully incorrect) version of Nakhleh distance.
    
    auto layers1 = transposed_kahn_ordering(net1);
    auto layers2 = transposed_kahn_ordering(net2);
    auto u1 = unique_nodes(net1, layers1);
    auto u2 = unique_nodes(net2, layers2);
    auto [eq1, eq2] = equivalent_nodes(net1, layers1, net2, layers2);
    int left_sum = 0;
    int right_sum = 0;
    
    for (const auto& [u, kappa] : u1) {
        left_sum += std::max(0, kappa-(int)eq1[u].size());
    }
    for (const auto& [v, kappa] : u2) {
        right_sum += std::max(0, kappa-(int)eq2[v].size());
    }

    return left_sum+right_sum;
}

std::unordered_map<Node*,std::unordered_map<std::string,unsigned int>> path_multiplicity_encoding(const RootedNetwork &net) {
    // https://doi.org/10.1109/TCBB.2007.70270
    auto leaves = net.get_leaf_labels(); // assuming the two networks have the same leaf sets
    auto layers = transposed_kahn_ordering(net);
    std::unordered_map<Node*,std::unordered_map<std::string,unsigned int>> path_multi;
    for (auto leaf : net.get_leaves()) {
        path_multi[leaf][leaf->label] += 1;
    }
    for (auto layer = layers.begin()+1; layer != layers.end(); ++layer) {
        for (auto node : *layer) {
            for (auto child : net.get_children(node)) {
                for (auto leaf : leaves) {
                    path_multi[node][leaf] += path_multi[child][leaf];
                }
            }
        }
    }

    return path_multi;
}

unsigned int path_multiplicity_distance(const RootedNetwork& net1, const RootedNetwork& net2) {
    // Distance for tree-child networks described in https://doi.org/10.1109/TCBB.2007.70270
    auto pm1 = path_multiplicity_encoding(net1);
    auto pm2 = path_multiplicity_encoding(net2);
    auto l1 = net1.get_leaf_labels();
    std::sort(l1.begin(), l1.end());
    auto l2 = net2.get_leaf_labels();
    std::sort(l2.begin(), l2.end());
    std::vector<std::string> common_leaves;
    std::set_intersection(l1.begin(), l1.end(), l2.begin(), l2.end(), std::back_inserter(common_leaves));

    std::vector<std::vector<std::pair<std::string, unsigned int>>> pmv1, pmv2, symdif;
    for (auto [node, path_multis] : pm1) {
        std::vector<std::pair<std::string,unsigned int>> path_multi;
        for (auto leaf : common_leaves) {
            path_multi.push_back(std::pair(leaf,path_multis[leaf]));
        }
        pmv1.push_back(path_multi);
    }
    for (auto [node, path_multis] : pm2) {
        std::vector<std::pair<std::string,unsigned int>> path_multi;
        for (auto leaf : common_leaves) {
            path_multi.push_back(std::pair(leaf,path_multis[leaf]));
        }
        pmv2.push_back(path_multi);
    }
    std::sort(pmv1.begin(), pmv1.end());
    std::sort(pmv2.begin(), pmv2.end());
    std::set_symmetric_difference(pmv1.begin(), pmv1.end(), pmv2.begin(), pmv2.end(), std::back_inserter(symdif));
    return symdif.size();
}

std::unordered_map<Node*,std::unordered_map<std::string, std::multiset<unsigned int>>> path_length_encoding(const RootedNetwork &net) {
    auto leaves = net.get_leaf_labels(); // assuming the two networks have the same leaf sets
    auto layers = transposed_kahn_ordering(net);
    std::unordered_map<Node*,std::unordered_map<std::string, std::multiset<unsigned int>>> path_lengths;

    for (auto leaf : net.get_leaves()) {
        path_lengths[leaf][leaf->label] = {0};
    }
    
    for (auto layer = layers.begin()+1; layer != layers.end(); ++layer) {
        for (auto node : *layer) {
            for (auto child : net.get_children(node)) {
                for (auto [leaf, pls] : path_lengths[child]) {
                    for (auto path_length : pls) {
                        path_lengths[node][leaf].insert(path_length+1);
                    }
                }
            }
        }
    }
    return path_lengths;
}

unsigned int path_length_distance(const RootedNetwork& net1, const RootedNetwork& net2) {
    // There is at least one distance metric with this name
    // However, it has more to do with pairwise distances between leaves
    // This is representing each node by the multiset of its path lengths to each leaf.
    //
    // It avoids the problem where if a node has no equivalent, the root
    // (and all descendants between the root and the non-equivalent node)
    // are not equivalent
    //
    // For example ((A,B),(C,D)); and ((A,C),(B,D));
    // The Nakhleh distance is 3, this path length distance is 2.
    // In practice, this often is almost identical to Nakhleh, with the
    // tradeoff of being less efficient.
    // 
    // This is for sure a metric on the same class of rooted binary phylogenetic networks as Nakhleh.
    // When I prove it's a metric on arbitrary-degree phylogenetic networks I'll throw it on arxiv or something
    auto ple1 = path_length_encoding(net1);
    auto ple2 = path_length_encoding(net2);
    auto l1 = net1.get_leaf_labels();
    std::sort(l1.begin(), l1.end());
    auto l2 = net2.get_leaf_labels();
    std::sort(l2.begin(), l2.end());
    std::vector<std::string> common_leaves;
    std::set_intersection(l1.begin(), l1.end(), l2.begin(), l2.end(), std::back_inserter(common_leaves));
    
    std::vector<std::vector<std::pair<std::string, std::multiset<unsigned int>>>> plv1, plv2, symdif;
    for (auto [node, path_multis] : ple1) {
        std::vector<std::pair<std::string,std::multiset<unsigned int>>> path_multi;
        for (auto leaf : common_leaves) {
            path_multi.push_back(std::pair(leaf,ple1[node][leaf]));
        }
        plv1.push_back(path_multi);
    }
    for (auto [node, path_multis] : ple2) {
        std::vector<std::pair<std::string,std::multiset<unsigned int>>> path_multi;
        for (auto leaf : common_leaves) {
            path_multi.push_back(std::pair(leaf,ple2[node][leaf]));
        }
        plv2.push_back(path_multi);
    }
    std::sort(plv1.begin(), plv1.end());
    std::sort(plv2.begin(), plv2.end());
    std::set_symmetric_difference(plv1.begin(), plv1.end(), plv2.begin(), plv2.end(), std::back_inserter(symdif));
    return symdif.size(); // possible to normalize by plv1.size() + plv2.size()? Kind of like Robinson-Foulds
}

template <class T> bool ckmin(T &a, const T &b) { return b < a ? a = b, 1 : 0; }

// Source: https://en.wikipedia.org/wiki/Hungarian_algorithm
template <class T> T hungarian(const std::vector<std::vector<T>> &C) {
        const int J = (int)size(C), W = (int)size(C[0]);
        assert(J <= W);
        std::vector<int> job(W + 1, -1);
        std::vector<T> ys(J), yt(W + 1);    // potentials
        std::vector<T> answers;
        const T inf = std::numeric_limits<T>::max();
        for (int j_cur = 0; j_cur < J; ++j_cur) {    // assign j_cur-th job
                int w_cur = W;
                job[w_cur] = j_cur;
                // min reduced cost over edges from Z to worker w
                std::vector<T> min_to(W + 1, inf);
                std::vector<int> prv(W + 1, -1);    // previous worker on alternating path
                std::vector<bool> in_Z(W + 1);        // whether worker is in Z
                while (job[w_cur] != -1) {     // runs at most j_cur + 1 times
                        in_Z[w_cur] = true;
                        const int j = job[w_cur];
                        T delta = inf;
                        int w_next;
                        for (int w = 0; w < W; ++w) {
                                if (!in_Z[w]) {
                                        if (ckmin(min_to[w], C[j][w] - ys[j] - yt[w]))
                                                prv[w] = w_cur;
                                        if (ckmin(delta, min_to[w])) w_next = w;
                                }
                        }
                        // delta will always be nonnegative,
                        // except possibly during the first time this loop runs
                        // if any entries of C[j_cur] are negative
                        for (int w = 0; w <= W; ++w) {
                                if (in_Z[w]) ys[job[w]] += delta, yt[w] -= delta;
                                else min_to[w] -= delta;
                        }
                        w_cur = w_next;
                }
                // update assignments along alternating path
                for (int w; w_cur != W; w_cur = w) job[w_cur] = job[w = prv[w_cur]];
                answers.push_back(-yt[W]);
        }
        return answers.back();
}

template <typename T> T earth_movers_distance(const std::multiset<T>& x, const std::multiset<T>& y) {
    T total = 0;
    if (x.size() < y.size()) {
        size_t sz_diff = y.size() - x.size();
        auto x_it = x.begin();
        auto y_it = y.begin();
        for (size_t ix = 0; ix < sz_diff; ++ix) {
            total += *(y_it);
            ++y_it;
        }
        for (size_t ix = 0; ix < x.size(); ++ix) {
            total += std::abs((((int)(*y_it)) - (int)(*x_it)));
            ++y_it;
            ++x_it;
        }
    } else {
        auto x_it = x.begin();
        auto y_it = y.begin();
        size_t sz_diff = x.size() - y.size();
        for (size_t ix = 0; ix < sz_diff; ++ix) {
            total += *(x_it);
            ++x_it;
        }
        for (size_t ix = 0; ix < y.size(); ++ix) {
            total += std::abs((((int)(*y_it)) - (int)(*x_it)));
            ++y_it;
            ++x_it;
        }
        
    }
    return total;
}

double path_length_hamming_distance(const RootedNetwork& net1, const RootedNetwork& net2) {
    // bad name.
    // hamming-distance of path lengths
    // Same idea as https://doi.org/10.1038/s41586-018-0043-0,
    // but applied to the path-length encoding.
    // Also, since the path length distance is a metric, this is too.
    //
    // In practice, this ends up "smoothing" the distribution of distances
    // Not really that useful.
    auto ple1 = path_length_encoding(net1);
    auto ple2 = path_length_encoding(net2);
    auto l1 = net1.get_leaf_labels();
    std::sort(l1.begin(), l1.end());
    auto l2 = net2.get_leaf_labels();
    std::sort(l2.begin(), l2.end());
    std::vector<std::string> common_leaves;
    std::set_intersection(l1.begin(), l1.end(), l2.begin(), l2.end(), std::back_inserter(common_leaves));
    
    std::vector<std::vector<std::multiset<unsigned int>>> plv1, plv2, symdif;
    for (auto [node, path_multis] : ple1) {
        if (net1.is_leaf(node))
            continue;
        std::vector<std::multiset<unsigned int>> path_multi;
        for (auto leaf : common_leaves) {
            path_multi.push_back(ple1[node][leaf]);
        }
        plv1.push_back(path_multi);
    }
    for (auto [node, path_multis] : ple2) {
        if (net2.is_leaf(node))
            continue;
        std::vector<std::multiset<unsigned int>> path_multi;
        for (auto leaf : common_leaves) {
            path_multi.push_back(ple2[node][leaf]);
        }
        plv2.push_back(path_multi);
    }
    double result = 0;
    
    // hm.
    std::vector<std::vector<double>> cost(
        plv1.size(),
        std::vector<double>(plv2.size(), 0)
    );
    std::vector<std::vector<double>> cost_T(
        plv2.size(),
        std::vector<double>(plv1.size(), 0)
    );

    size_t idx_2, idx_1 = 0;
    for (auto ple_1 : plv1) {
        idx_2 = 0;
        for (auto ple_2 : plv2) {
            for (size_t ix=0; ix < common_leaves.size(); ++ix) {
                if ((ple_1[ix] != ple_2[ix])) {
                   cost[idx_1][idx_2] += 1;
                    cost_T[idx_2][idx_1] += 1;
                }
            }
            ++idx_2;
        }
        ++idx_1;
    }

    if (plv1.size() < plv2.size()) {
        result += 2*hungarian(cost)/common_leaves.size() + (plv2.size() - plv1.size());
    } else {
        result += 2*hungarian(cost_T)/common_leaves.size() + (plv1.size() - plv2.size());
    }
    
    return result; // normalize? maximum distance formulation is going to be complicated
}

double path_adjustment_distance(const RootedNetwork& net1, const RootedNetwork& net2) {
    // same as above, but using a earth-mover's instead of hamming
    auto ple1 = path_length_encoding(net1);
    auto ple2 = path_length_encoding(net2);
    auto l1 = net1.get_leaf_labels();
    std::sort(l1.begin(), l1.end());
    auto l2 = net2.get_leaf_labels();
    std::sort(l2.begin(), l2.end());
    std::vector<std::string> common_leaves;
    std::set_intersection(l1.begin(), l1.end(), l2.begin(), l2.end(), std::back_inserter(common_leaves));
    
    std::vector<std::vector<std::multiset<unsigned int>>> plv1, plv2, symdif;
    for (auto [node, path_multis] : ple1) {
        if (net1.is_leaf(node))
            continue;
        std::vector<std::multiset<unsigned int>> path_multi;
        for (auto leaf : common_leaves) {
            path_multi.push_back(ple1[node][leaf]);
        }
        plv1.push_back(path_multi);
    }
    for (auto [node, path_multis] : ple2) {
        if (net2.is_leaf(node))
            continue;
        std::vector<std::multiset<unsigned int>> path_multi;
        for (auto leaf : common_leaves) {
            path_multi.push_back(ple2[node][leaf]);
        }
        plv2.push_back(path_multi);
    }
    double result = 0;
    
    // hm.
    std::vector<std::vector<double>> cost(
        plv1.size(),
        std::vector<double>(plv2.size(), 0)
    );
    std::vector<std::vector<double>> cost_T(
        plv2.size(),
        std::vector<double>(plv1.size(), 0)
    );

    size_t idx_2, idx_1 = 0;
    for (auto ple_1 : plv1) {
        idx_2 = 0;
        for (auto ple_2 : plv2) {
            for (size_t ix=0; ix < common_leaves.size(); ++ix) {
                if ((ple_1[ix] != ple_2[ix])) {
                    // Earth-mover's distance.
                    unsigned int emd = earth_movers_distance(ple_1[ix], ple_2[ix]);
                    unsigned int num_paths = ple_1.size() + ple_2.size();
                    cost[idx_1][idx_2] += (double)emd / num_paths; //(double)symdif.size()/uni.size();
                    cost_T[idx_2][idx_1] += (double)emd / num_paths; // (double)symdif.size()/uni.size(); //1; // 1.0/common_leaves.size(); // cost[idx_1][idx_2];
                }
            }
            ++idx_2;
        }
        ++idx_1;
    }

    if (plv1.size() < plv2.size()) {
        result += 2*hungarian(cost)/common_leaves.size() + (plv2.size() - plv1.size());
    } else {
        result += 2*hungarian(cost_T)/common_leaves.size() + (plv1.size() - plv2.size());
    }
    
    return result;
}

std::string to_newick(const RootedNetwork& net, const int format, const int precision) {
    std::unordered_map<Node*,std::string> labels;
    bool internal_labels, branch_lengths, support, admixture;
    // TODO: flexible formatting?
    switch(format) {
        case 0: // make a guess as to what the format should be
            internal_labels = true;
            branch_lengths = false;
            support = false;
            admixture = false;
            for (auto e : net.get_edges()) {
                if (e->distance != 0)
                    branch_lengths = true;
                if (e->support != 0)
                    support = true;
                if (e->admixture != 1)
                    admixture = true;
            }
            break;
        case 1: // suppress support, write branch lengths and internal labels
            internal_labels = true;
            branch_lengths = true;
            support = false;
            admixture = true;
            break;
        case 2: // 
            break;
        default:            
            throw std::out_of_range("");
    }
    
    std::set<Node*> visited_hybrids;
    auto layers = transposed_kahn_ordering(net);
    for (auto leaf : layers[0]) {
        labels[leaf] = leaf->label;
        if (branch_lengths) {
            // maybe handle the case where a leaf has in-degree 2? it shouldn't... but
            labels[leaf] += std::format(":{:.{}f}", net.get_incoming_edges(leaf)[0]->distance, precision);
        }
    }

    for (size_t ix = 1; ix < layers.size(); ++ix) {
        for (auto node : layers[ix]) {
            labels[node] = "(";
            for (auto child : net.get_children(node)) {
                if (net.is_hybrid(child)) {
                    if (!visited_hybrids.contains(child)) {
                        visited_hybrids.insert(child);
                    } else {
                        labels[node] += child->label;
                        Edge* hybrid_edge = net.get_edge(node, child);
                        if (branch_lengths)
                            labels[node] += std::format(":{:.{}f}", hybrid_edge->distance, precision);
                        if (support)
                            labels[node] += std::format(":{:.{}f}", hybrid_edge->support, precision);
                        if (admixture) {
                            if (!support)
                                labels[node] += ':';
                            labels[node] += std::format(":{:.{}f}", hybrid_edge->admixture, precision);
                        }
                        labels[node] += ',';
                        continue;
                    }
                }
                labels[node] += labels[child] + ',';
            }
            labels[node][labels[node].size()-1] = ')';
            if (internal_labels || net.is_hybrid(node)) {
                labels[node] += node->label;
            }
            if (branch_lengths && !net.is_root(node)) {
                // Need to test this, in particular for multiple sources/sinks
                Edge* parent_edge = nullptr;
                auto parent_edges = net.get_incoming_edges(node);
                if (parent_edges.size() == 1) {
                    parent_edge = parent_edges[0];
                } else {
                    // get the first node in the above layers that is the parent of this one.
                    auto parents = net.get_parents(node);
                    std::set<Node*> parent_set(parents.begin(), parents.end());
                    for (auto jx = ix+1; jx < layers.size(); ++jx) {
                        for (auto potential_parent : layers[jx]) {
                            if (parent_set.contains(potential_parent)) {
                                parent_edge = net.get_edge(potential_parent, node);
                                break;
                            }
                        }
                        if (parent_edge != nullptr) {
                            break;
                        }
                    }
                }
                if (parent_edge == nullptr) {
                    throw std::out_of_range(std::format("Could not find a parent for hybrid node labeled {}", node->label));
                }
                labels[node] += std::format(":{:.{}f}", parent_edge->distance, precision);
                if (parent_edges.size() > 1) {
                    if (!support)
                        labels[node] += ':';
                    labels[node] += std::format(":{:.{}f}", parent_edge->admixture, precision);
                }
            }
        }
    }
    return labels[net.get_root()] + ";";
}


