#include <print>
#include <iostream>
#include <vector>
#include <string>
#include <random>

#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options.hpp>

#include "RootedNetwork.hpp"
#include "Comparisons.hpp"

namespace po = boost::program_options;

int main_compare(int argc, char *argv[], po::variables_map &global_vm);
int main_reticulate(int argc, char *argv[], po::variables_map &global_vm);
std::vector<RootedNetwork> read_newicks(po::variables_map &vm);

int main(int argc, char *argv[]) {
    try {
        po::options_description global_opts("Global options");
        global_opts.add_options()
            ("help,h", "Print help")
            ("command", po::value<std::string>(), "Command (compare, reticulate)")
            ("subargs", po::value<std::vector<std::string>>(), "Arguments for command")
        ;
        po::positional_options_description pos;
        pos.add("command", 1).
                add("subargs", -1);

        po::variables_map vm;
        po::parsed_options parsed = po::command_line_parser(argc, argv).
            options(global_opts).
            positional(pos).
            allow_unregistered().
            run();
        
        po::store(parsed, vm);
        po::notify(vm);
        
        if (vm.count("command")) {
            std::string cmd = vm["command"].as<std::string>();
            if (cmd == "compare") {
                return main_compare(argc, argv, vm);
            } else if (cmd == "reticulate") {
                return main_reticulate(argc, argv, vm);
            } else {
                std::println("subcommand \"{}\" not recognized.", cmd);
                std::cout << global_opts;
                return 0;
            }
        }
        else if (vm.count("help")) {
            std::cout << global_opts;
            return 0;
        }
        
    } catch (const po::error &ex) {
        std::cerr << ex.what() << '\n';
        return 1;
    } 
}

int main_compare(int argc, char *argv[], po::variables_map &global_vm) {
    po::options_description cmp_opts("Comparison options");

    cmp_opts.add_options()
        ("newick,n", po::value<std::vector<std::string>>()->multitoken(), "Trees/networks to compare in newick format")
        ("measure,m", po::value<std::vector<std::string>>()->multitoken(), "Measure(s) for network distance/dissimilarity to use (options: mu, nakhleh, pl)")
    ;

    if (global_vm.count("help")) {
        std::cout << cmp_opts;
        return 0;
    }

    po::parsed_options parsed = po::command_line_parser(argc, argv).
        options(cmp_opts).
        run();

    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);


    auto nets = read_newicks(vm);
    if (nets.size() == 0) {
        std::println("No networks supplied");
        return 1;
    }
    if (nets.size() == 1) {
        std::println("Need at least two networks to compare, only one was supplied");
        return 1;
    }
    size_t max_ix = nets.size(); // TODO: add option to change this to only compare one network to the rest.

    std::vector<std::function<double(const RootedNetwork&, const RootedNetwork&)>> dist_fns;
    std::vector<std::string> dist_names;
    if (vm.count("measure")) {
        auto measures = vm["measure"].as<std::vector<std::string>>();
        for (std::string measure : measures) {
            std::transform(measure.begin(), measure.end(), measure.begin(), ::tolower);
            auto dist_fn = &(normalized_robinson_foulds);
            if (measure == "rf") {
                dist_fn = [](const RootedNetwork& a, const RootedNetwork& b) {return (double)robinson_foulds(a,b) ;} ;
                dist_names.push_back("Robinson-Foulds");
            } else if (measure == "nrf") {
                dist_fn = &(normalized_robinson_foulds);
                dist_names.push_back("Normalized Robinson-Foulds");
            } else if (measure == "nakhleh") {
                dist_fn = [](const RootedNetwork& a, const RootedNetwork& b) { return (double)nakhleh_distance(a,b) / 2.0; };
                dist_names.push_back("Nakhleh distance");
            } else if (measure == "pl") {
                dist_fn = [](const RootedNetwork& a, const RootedNetwork& b) { return (double)path_length_distance(a,b) / 2.0; };
                dist_names.push_back("Path length distance");
            } else if (measure == "pm" || measure == "mu") {
                dist_fn = [](const RootedNetwork& a, const RootedNetwork& b) { return (double)path_multiplicity_distance(a,b) / 2.0; };
                dist_names.push_back("mu-distance"); // μ
            } else {
                std::println("\nMeasure \"{}\" not recognized.", measure);
                return 1;
            }
            dist_fns.push_back(dist_fn);
        }
    } else {
        dist_names.push_back("Nakhleh distance");
        dist_fns.push_back([](const RootedNetwork& a, const RootedNetwork& b) { return (double)nakhleh_distance(a,b) / 2.0; });
    }

    if (nets.size() == 2) {
        for (size_t ix = 0; ix < dist_fns.size(); ++ix) {
            double dist = dist_fns[ix](nets[0], nets[1]);
            std::println("{}: {}", dist_names[ix], dist);
        }
        return 0;
    }
    
    std::print("Network 1,Network 2");
    for (auto dist_name : dist_names) {
        std::print(",{}", dist_name);
    }
    std::println();

    for (size_t ix=0; ix < max_ix; ++ix) {
        for (size_t jx=ix+1; jx < nets.size(); ++jx) {
            std::print("{},{}", ix+1, jx+1);
            for (auto dist_fn : dist_fns) {
                double dist = dist_fn(nets[ix], nets[jx]);
                std::print(",{}", dist);
            }
            std::println();
        }
    }

    return 0;
}

int main_reticulate(int argc, char *argv[], po::variables_map &global_vm) {
    po::options_description ret_opts("Reticulation options");
    ret_opts.add_options()
        ("newick,n", po::value<std::vector<std::string>>()->multitoken(), "Trees/networks to add reticulations to in extended newick format")
        ("number,r", po::value<unsigned int>(), "Number of reticulations to add")
        ("seed,s", po::value<int>(), "Seed for random number generator")
    ;

    if (global_vm.count("help")) {
        std::cout << ret_opts;
        return 0;
    }
        
    po::parsed_options parsed = po::command_line_parser(argc, argv).
        options(ret_opts).
        run();

    po::variables_map vm;
    po::store(parsed, vm);
    po::notify(vm);

    auto nets = read_newicks(vm);
    if (nets.size() == 0) {
        std::println("No networks supplied.");
        return 1;
    }
    
    if (!vm.count("number")) {
        std::cout << "Number of reticulations not specified\n";
        return 1;
    }
    size_t num_ret = vm["number"].as<unsigned int>();

    int seed;
    if (vm.count("seed")) {
        seed = vm["seed"].as<int>();
    } else {
        seed = std::random_device{}();
    }
    std::mt19937 rng(seed);
    
    for (auto net : nets) {
        size_t init_rets = net.num_hybrids();
        for (size_t ix=0; ix<num_ret; ++ix) {
            auto edges = net.get_edges();
            std::uniform_int_distribution<std::size_t> idx_dist(0, edges.size() - 1);
            size_t idx, jdx;
            idx = idx_dist(rng);
            jdx = idx_dist(rng);
            while (idx == jdx) jdx = idx_dist(rng); // pretend you don't see this.
            auto e1 = edges[idx];
            auto e2 = edges[jdx];
            // check if connecting these two edges induces a cycle
            // ie e1=(u1,v1), e2=(u2,v2)
            // is there a path from v2 to u1, then this introduces a cycle
            if (net.is_descendant(e1->parent, e2->child)) {
                ix--;
                continue;
            }
      
            Node* attach_1 = net.insert_attachment_point(e1);
            Node* attach_2 = net.insert_attachment_point(e2);
            attach_2->label = std::format("#H{}", init_rets+ix+1);
            net.insert_edge(attach_1, attach_2);
        }
        std::cout << to_newick(net) << '\n';
    }

    return 0;
}

std::vector<RootedNetwork> read_newicks(po::variables_map &vm) {
    std::vector<RootedNetwork> nets;
    if (!vm.count("newick")) {
        return nets;
    }
    
    for (auto x: vm["newick"].as<std::vector<std::string>>()) {
        // TODO: Handle file paths
        nets.push_back(parse_network(x));
    }
    return nets;
}
