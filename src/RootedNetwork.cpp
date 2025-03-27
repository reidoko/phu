#include "RootedNetwork.hpp"

#include <boost/tokenizer.hpp>

#include <stack>
#include <queue>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <utility>
#include <format>


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


