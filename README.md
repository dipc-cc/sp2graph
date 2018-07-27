[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

![sp2graph logo](/doc/images/sp2graph_logo.png "sp2graph code")

The goal of _sp2graph_ project is to provide tools for analyzing the bond-order of *sp<sup>2</sup>*-hybridized carbon nanostructures. This is particular interesting for the interpretation of recent experiments performed with scanning probe microscopy on graphene-based nanostructures synthesized on metallic surfaces.

## Objectives ##
Given the `xyz` geometry of a planar carbon-based structure (*sp<sup>2</sup>*-hybridized), __sp2graph__ aims to provide:

   1. all possible Kekulé representations
   2. the most stable Clar structure according to Clar's &pi;-sextet theory
   3. allow for user defined initial constrains, by imposing:
      - single or double bonds at specified connections
      - radicals at specific sites
   4. estimation of the most stable Kekulé structure through simple nearest-neighbor tight-binding model

## Introduction ##

Carbon *sp<sup>2</sup>*-bonded structures, in particular polycyclic aromatic hydrocarbons (PAHs), are often represented graphically according to the Kekulé bond formulas, where the four valence electrons of a carbon can form single or double bonds with neighboring atoms (see figures below). In this representation only the carbon-carbon bonds are drawn, being carbon-hydrogen bonds typically omitted (i.e. carbons with only three bonds or less are implicitly assumed to make extra bonds with hydrogen atoms). Structures exhibiting a radical character, such as open-shell states, are usually marked with a dot next to the carbons containing unpaired &pi;-electrons.

The Kekulé representation of any given *sp<sup>2</sup>* carbon-based system is not unique. The benzene molecule C<sub>6</sub>H<sub>6</sub>, for example, has two possible representations as shown in the left side of the figure below. Therefore, the electronic structure of a system is given by a combination of all possible Kekulé bond formulas [(Wassmann et al, 2010)][Wassmann2010]. For benzene the two complementary Kekulé structures form an aromatic &pi;-sextet (right side of the figure below), that represents the delocalization of the six &pi;-electrons.

![benzene resonant structures](/doc/images/benzene.png)  
**benzene:** *Kekulé resonant structures (left) and the corresponding Clar sextet (right)*.

Aromatic sextets are also known as Clar's sextets since the formulation of the Clar's &pi;-sextet rule which states that the representation with largest number of disjoint aromatic sextets (known as Clar's formula or structure) is the most stable configuration and, therefore, the most representative for characterizing the system properties [(Solà, 2013)][Sola2013]. A given system can present different equivalent Clar's formulas, i.e. different representations having the same (largest) number of aromatic sextet's. In such case, the system is described by a superposition of all equivalent Clar structures. A typical example is the anthracene molecule shown below, where a Clar's sextet can be found in one of the three rings.

![anthracene resonant structures](/doc/images/anthracene.png)  
**anthracene:** *Kekulé resonant structures (left) with the corresponding Clar sextets (right)*.

For the phenanthrene, on the other hand, only one Clar structure can be found with two aromatic sextets (see image below). The Clar's rule can be further applied, in absence of radicals, to find the most stable structure among alternant PAHs [(Yeh and Chai, 2016)][Yeh2016], i.e. isomers containing the same number of six-membered rings. For instance, by comparing the Clar structures from anthracene and phenanthrene one would expect the later to be more stable. Indeed this is in agreement with standard density functional theory calculations (DFT), where the phenanthrene is found to be more stable by 0.2 eV.

![phenanthrene resonant structures](/doc/images/phenanthrene.png)  
**phenanthrene:** *Kekulé resonant structures (left) with the corresponding Clar sextets (right)*.
NB: There should be four Kekulé structures for the lower Clar structure with two sextets.

Interestingly, by performing structural relaxations with simple reactive empirical bond order (REBO) potentials optimized for hydrocarbons [(Brenner et al, 2002)][Brenner2002] one would predict wrongly the total energy of anthracene and phenanthrene to be approximately equal, because the local bonding configurations of the carbon atoms are nearly identical. Therefore, this simple example illustrates the strength of the Clar's &pi;-sextet rule.

Although conceptually simple, apply the Clar's sextet to larger *sp<sup>2</sup>* structures can be rather complex. Furthermore, when drawing Clar structures to maximize the number of sextets, one may encounter configurations with carbon atoms with radical character. As a rule of thumb, three additional sextets need to be added to the Clar structure to compensate for the energy cost associated with creating two radicals [(Das and Wu, 2015)][Das2015]. For extended structures these rules leave many possibilities and thus call for a systematic and automated approach.

## Methodology ##

The Kekulé diagrams shown in the figures above can be described mathematically by a graph *G(V,E)*, a fundamental combinatorial object defined by a set of vertices *V* and a set of edges *E* connecting distinct vertex pairs. In our case, the vertices are given by the carbon atoms and the edges by the chemical bonds among them.

A simple method for representing a graph is the so-called adjacency matrix *A*, with dimensions *VxV* and where the element *A<sub>ij</sub>* is nonzero only if there is an edge connecting from vertex *i* to *j*. For instance, the naphthalene molecule could be represent as:

![naphthalene molecule with graph and adjacency matrix representations](/doc/images/naphthalene_graph.png)  
**Graph representation of a *sp<sup>2</sup>* structure:** *A naphthalene molecule with respective graph and adjacency matrix representations*.

In the adjacency matrix above, a carbon-carbon bond between vertices *i* and *j* is expressed by having *A<sub>ij</sub>*=1, while for the absence of bond as *A<sub>ij</sub>*=0. The graph of any *sp<sup>2</sup>* structure in general has a number of properties, to list a few:

   - it is an undirected graph, meaning that an edge is defined by a pair of distinct vertices independently of their order (*ij* or *ji*), which is expressed by the symmetry of the adjacency matrix;
   - any vertex is connected to 2 or 3 other vertices (i.e., any vertex has 2 or 3 incident edges), which makes the adjacency matrix representation sparse;
   - any vertex belongs to at least one cycle or, in other words, a vertex will never be hanging alone from the graph since this would mean a *sp<sup>1</sup>*-hybridized carbon.

More generally, given the set of `xyz` coordinates of a *sp<sup>2</sup>* structure, the set of vertex and edges can be uniquely defined. However, the adjacency matrix defined above contain no information about where are the single and double bonds. There are several ways to include such information, and a simple one would be to define *A<sub>ij</sub>*=2 if there is a double bond between vertices *i* and *j*. The edges containing a double bond are not uniquely defined, for the same reason a *sp<sup>2</sup>* system has different Kekulé structures. According to this definition, one possibility for the naphthalene molecule would be:

![naphthalene Kekulé structure with graph and adjacency matrix representations](/doc/images/naphthalene_kekule.png)  
**Graph representation of a Kekulé structure:** *A naphthalene Kekulé structure with respective graph and adjacency matrix representations*.

Therefore, in addition to the properties listed before one can include:

   - any vertex is connected to at least another vertex through a double bond, i.e., for any row *i* from the adjacency matrix there is one and only one element with *A<sub>ij</sub>*=2.



## Installation ##
For __sp2graph__ installation the following packages are required:

   - numpy

Manual installation of  is performed with the command:

    python setup.py install --prefix=<prefix>
    # or
    python setup.py install --home=<my-python-home>

One may also wish to set the following environment variables:

    export PYTHONPATH=<my-python-home>/lib/<python-version>/site-packages
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<my-python-home>/lib

## Contributions, issues and bugs ##
Contributions are highly appreciated!

If you find any bug please form a [bug report/issue][issues].

If you have a fix please consider adding a [pull request][pulls].

## License ##
The __sp2graph__ license is [LGPL][lgpl], please see the LICENSE file.

<!---
Links to external and internal sites.
-->
[lgpl]: http://www.gnu.org/licenses/lgpl.html
[issues]: https://github.com/dipc-cc/sp2graph/issues
[pulls]: https://github.com/dipc-cc/sp2graph/pulls
[Wassmann2010]: https://doi.org/10.1021/ja909234y
[Sola2013]: https://doi.org/10.3389/fchem.2013.00022
[Yeh2016]: https://doi.org/10.1038/srep30562
[Brenner2002]: https://doi.org/10.1088/0953-8984/14/4/312
[Das2015]: https://doi.org/10.1002/9783527689545.ch1
