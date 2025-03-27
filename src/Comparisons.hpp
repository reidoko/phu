#pragma once

#include "RootedNetwork.hpp"

// trees
unsigned int robinson_foulds(const RootedNetwork&, const RootedNetwork&);
double normalized_robinson_foulds(const RootedNetwork&, const RootedNetwork&);

// networks
unsigned int nakhleh_distance(const RootedNetwork&, const RootedNetwork&);
unsigned int path_multiplicity_distance(const RootedNetwork&, const RootedNetwork&);
unsigned int path_length_distance(const RootedNetwork&, const RootedNetwork&);
