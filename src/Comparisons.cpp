#include "Comparisons.hpp"

#include <set>
#include <algorithm>
#include <tuple>
#include <limits>
#include <ranges>
#include <cassert>

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

