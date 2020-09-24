# CATH Tools [![Build Status](https://travis-ci.com/UCLOrengoGroup/cath-tools.svg?branch=master)](https://travis-ci.com/UCLOrengoGroup/cath-tools) [![Documentation Status](https://readthedocs.org/projects/cath-tools/badge/?version=latest)](https://readthedocs.org/projects/cath-tools/?badge=latest)

Protein structure comparison tools such as SSAP, as used by the [Orengo Group](https://www.ucl.ac.uk/orengo-group "Orengo Group website") in curating [CATH](http://www.cathdb.info/ "CATH website").

| | | | | |
|:-- |:-- |:-- |:-- |:-- |
| [**Executable DOWNLOADS**](https://github.com/UCLOrengoGroup/cath-tools/releases/latest "Downloads for the latest CATH Tools release") <br> <sup> (for Linux/Mac; chmod them to be executable)</sup> | [**Docs**](http://cath-tools.readthedocs.io/en/latest/ "CATH Tools user documentation at Read the Docs") <br> <sup> &nbsp; </sup> | [**Code**](https://github.com/UCLOrengoGroup/cath-tools "CATH Tools GitHub repository") <br> <sup> &nbsp; </sup> | [Builds](https://travis-ci.com/UCLOrengoGroup/cath-tools "CATH Tools Travis-CI builds") <br> <sup> &nbsp; </sup> | [Extras repo](https://github.com/UCLOrengoGroup/cath-tools-supplementary "CATH Tools Supplementary GitHub repository") <br> <sup> &nbsp; </sup> |


### Tools

| | |
|:-- |:-- |
| <img class="tool-thumb" style="border-style: solid; border-width: 1px;" src="https://raw.githubusercontent.com/UCLOrengoGroup/cath-tools/master/docs/tools/img/cath-cluster.jpg"         width="400" /> | [ **cath-cluster** ]( http://cath-tools.readthedocs.io/en/latest/tools/cath-cluster ) <br> Complete-linkage cluster arbitrary data. |
| <img class="tool-thumb" style="border-style: solid; border-width: 1px;" src="https://raw.githubusercontent.com/UCLOrengoGroup/cath-tools/master/docs/tools/img/cath-map-clusters.jpg"         width="400" /> | [ **cath-map-clusters** ]( http://cath-tools.readthedocs.io/en/latest/tools/cath-map-clusters ) <br> Map names from previous clusters to new clusters based on (the overlaps between) their members (which may be specified as regions within a parent sequence). Renumber any clusters with no equivalents. |
| <img class="tool-thumb" style="border-style: solid; border-width: 1px;" src="https://raw.githubusercontent.com/UCLOrengoGroup/cath-tools/master/docs/tools/img/cath-resolve-hits.example.jpg" width="400" /> | [ **cath-resolve-hits** ]( http://cath-tools.readthedocs.io/en/latest/tools/cath-resolve-hits ) <br> Collapse a list of domain matches to your query sequence(s) down to the non-overlapping subset (ie domain architecture) that maximises the sum of the hits' scores. |
| <img class="tool-thumb" style="border-style: solid; border-width: 1px;" src="https://raw.githubusercontent.com/UCLOrengoGroup/cath-tools/master/docs/tools/img/ssap_alignment.jpg"            width="400" /> | [ **cath-ssap**         ]( http://cath-tools.readthedocs.io/en/latest/tools/cath-ssap         ) <br> Structurally align a pair of proteins. |
| <img class="tool-thumb" style="border-style: solid; border-width: 1px;" src="https://raw.githubusercontent.com/UCLOrengoGroup/cath-tools/master/docs/tools/img/1fi2A00_1j58A01.good.jpg"      width="400" /> | [ **cath-superpose**    ]( http://cath-tools.readthedocs.io/en/latest/tools/cath-superpose    ) <br> Superpose two or more protein structures using an existing alignment. |

### Extra Tools

 * `build-test`          Perform the cath-tools tests (which should all pass, albeit with a few warnings)
 * `cath-assign-domains` Use an SVM model on SSAP+PRC data to form a plan for assigning the domains to CATH superfamilies/folds
 * `cath-refine-align`   Iteratively refine an existing alignment by attempting to optimise SSAP score
 * `cath-score-align`    Score an existing alignment using structural data

<!--
| |
|:-- |:-- |
| `cath-assign-domains` | Use an SVM model on SSAP+PRC data to form a plan for assigning the domains to CATH superfamilies/folds |
| `cath-refine-align`   | Iteratively refine an existing alignment by attempting to optimise SSAP score                          |
| `cath-score-align`    | Score an existing alignment using structural data                                                      |
| `check-pdb`           | Check a PDB file for some potential problems                                                           |
-->

Authors
-------

The SSAP algorithm (`cath-ssap`) was devised by Christine A Orengo and William R Taylor.

Please cite: *Protein Structure Alignment*, Taylor and Orengo, Journal of Molecular Biology 208, 1-22, PMID: 2769748. ([PubMed](https://www.ncbi.nlm.nih.gov/pubmed/2769748), [Elsevier](http://www.sciencedirect.com/science/article/pii/0022283689900843))

Since then, many people have contributed to this code, most notably:

  * [Tony E Lewis](https://github.com/tonyelewis)             (2011&ndash;&hellip;)
  * Oliver C Redfern                                          (~2003&ndash;2011)
  * James E Bray, [Ian Sillitoe](https://github.com/sillitoe) (~2000&ndash;2003)
  * [Andrew C R Martin](https://github.com/AndrewCRMartin)    (considerable edits around 2001)

Acknowledgements
----------------

`cath-ssap` typically uses [DSSP](https://github.com/cmbi/xssp), either by reading DSSP files or via its own implementation of the DSSP algorithms.

`cath-cluster` uses Fionn Murtagh's reciprocal-nearest-neighbour algorithm (see [Multidimensional clustering algorithms, volume 4 of Compstat Lectures.
Physica-Verlag, Würzburg/ Wien, 1985. ISBN 3-7051-0008-4](http://www.multiresolutions.com/strule/MClA/)) as described and refined in Daniel Müllner's [Modern hierarchical, agglomerative clustering algorithms](https://arxiv.org/abs/1109.2378) (2011, arXiv:1109.2378).

Feedback
--------

Please tell us about your cath-tools bugs/suggestions [here](https://github.com/UCLOrengoGroup/cath-tools/issues/new).

If you find this software useful, please spread the word and star the [GitHub repo](https://github.com/UCLOrengoGroup/cath-tools).
