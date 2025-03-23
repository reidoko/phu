#include <string>
#include <vector>
#include <unordered_map>

struct Node {
    std::string label = "";
};

struct Edge {
    double distance = 0;
    double support = 0;
    double admixture = 1;
    Node *parent;
    Node *child;
};

class RootedNetwork {
 public:
    RootedNetwork();
    // Rule of 5
    ~RootedNetwork(); // I. destructor
    RootedNetwork(const RootedNetwork& other); // II. copy constrctor
    RootedNetwork& operator=(const RootedNetwork& other); // III. copy assignment
    RootedNetwork(RootedNetwork&& other) noexcept; // IV. move constructor
    RootedNetwork& operator=(RootedNetwork&& other) noexcept; // V. move assignment

    Node* get_root() const { return root; }
    Edge* get_edge(Node*, Node*) const;
    std::vector<Node*> get_nodes() const { return nodes; }
    std::vector<Edge*> get_edges() const { return edges; }
    std::vector<Node*> get_leaves() const;
    std::vector<Node*> get_leaves(Node*) const;
    std::vector<std::string> get_leaf_labels() const;
    std::vector<Node*> get_reticulation_nodes() const;
    std::vector<Node*> get_tree_nodes() const;
    std::vector<Node*> get_internal_nodes() const;
    std::vector<Edge*> get_incoming_edges(Node* u) const { return incoming.at(u); }
    std::vector<Edge*> get_outgoing_edges(Node* u) const { return outgoing.at(u); }
    std::vector<Node*> get_parents(Node*) const;
    std::vector<Node*> get_children(Node*) const;
    std::vector<Node*> get_siblings(Node*) const;
    unsigned int num_nodes() const { return nodes.size(); }
    unsigned int num_edges() const { return edges.size(); }
    unsigned int num_hybrids() const { return get_reticulation_nodes().size(); }
    unsigned int num_leaves() const { return get_leaves().size(); }
    unsigned int indegree(Node* u) const { return incoming.at(u).size(); }
    unsigned int outdegree(Node* u) const { return outgoing.at(u).size(); }
    bool is_root(Node* u) const { return u == root; }
    bool is_leaf(Node* u) const { return outdegree(u) == 0; }
    bool is_hybrid(Node* u) const { return indegree(u) > 1; }
    bool is_parent(Node*, Node*) const;
    bool is_child(Node*, Node*) const;
    bool is_descendant(Node*, Node*) const;
    bool adjacent(Node*, Node*) const;

    void delete_node(Node*);
    void insert_node(Node*);
    Node* insert_attachment_point(Edge*, double=0.5);
    void merge_nodes(Node*, Node*);

    void insert_edge(Edge*);
    Edge* insert_edge(Node*, Node*);
    void delete_edge(Edge*);

 private:
    Node* root;
    std::vector<Node*> nodes;
    std::vector<Edge*> edges;
    std::unordered_map<Node*,std::vector<Edge*>> incoming;
    std::unordered_map<Node*,std::vector<Edge*>> outgoing;
};

RootedNetwork parse_network(std::string);
std::string to_newick(const RootedNetwork&, const int format=0, const int precision=4);

// distance measures
// trees
unsigned int robinson_foulds(const RootedNetwork&, const RootedNetwork&);
double normalized_robinson_foulds(const RootedNetwork&, const RootedNetwork&);

// networks
unsigned int nakhleh_distance(const RootedNetwork&, const RootedNetwork&);
unsigned int path_multiplicity_distance(const RootedNetwork&, const RootedNetwork&);
unsigned int path_length_distance(const RootedNetwork&, const RootedNetwork&);
