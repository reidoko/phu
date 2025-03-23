# phu
A few utilities for phylogenetic networks and trees.

## Usage

### Compare
Calculate the distance or dissimilarity of two or more networks.

Available comparisons for networks:
- `mu`:  Path multiplicity distance (also called μ-distance, from [Cardona et al., 2008](https://doi.org/10.1109/TCBB.2007.70270))
- `nakhleh`: Nested-labels distance (also called Nakhleh distance, from [Luay Nakhleh, 2009](https://doi.org/10.1109/TCBB.2009.2))
- `pl`: Path-lengths distance (see [here](docs/distances.md) for more details)

```
./src/phu compare -n "((B,(A)#H1),((C,E),(D,#H1)));" "((B,(A)#H1),((C,D),(E,#H1)));" -m mu nakhleh pl
mu-distance: 2
Nakhleh distance: 4
Path length distance: 2
```

### Reticulate
Add random reticulations to a network.
```
./src/phu reticulate -n "((A,B),(C,D));" -r 1 --seed 1234
(((A,(B)#H1),#H1),(C,D));
```
