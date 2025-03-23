# $\mu$-distance
Originally from https://doi.org/10.1109/TCBB.2007.70270
# Nakhleh distance
Originally from https://doi.org/10.1109/TCBB.2009.2

# Path-lengths distance
There has been a lot of work in nodal distances for phylogenetic trees,
but they usually look at pairwise leaf-to-leaf distances. This is
taking the ideas from that, but applying them to the approach
Nakhleh distance takes for denoting node equivalence.

In this distance, a node in a rooted phylogenetic network is encoded
as a vector of path lengths for all distinct directed paths from that node to
the set of leaves, and counts the distance as the symmetric difference
of the set of node encodings for the two networks.

As an example, consider the following tree ((A,B),(C,D));  
The encodings would be as follows: 
```
A: {A:{0}}  
B: {B:{0}}  
C: {C:{0}}  
D: {D:{0}}  
(A,B): {A:{1}, B:{1}}  
(C,D): {C:{1}, D:{1}}  
((A,B),(C,D)): {A:{2}, B:{2}, C:{2}, D:{2}}
```
The tree ((A,C),(B,D)); is encoded by
```
A: {A:{0}}  
B: {B:{0}}  
C: {C:{0}}  
D: {D:{0}}  
(A,C): {A:{1}, C:{1}}  
(B,D): {B:{1}, D:{1}}  
((A,B),(C,D)): {A:{2}, B:{2}, C:{2}, D:{2}}  
```
And the symmetric difference of these node encodings would be
{{A:{1}, C:{1}}, {A:{1}, B:{1}}, {C:{1},D:{1}}, {B:{1},D:{1}}.
Or to 2 non-equivalent nodes in each tree.
The Nakhleh distance would have a distance of 3 here,
as the root would also not be equivalent since none of its children are equivalent.

This is a distance metric on rooted binary phylogenetic networks
with the same separation power as Nakhleh distance. I have not yet found
a counter-example for arbitrary rooted phylogenetic networks, nor proven
one way or the other if this has the same separation power as Nakhleh.
