# CATH Tools [![Build Status](https://travis-ci.org/UCLOrengoGroup/cath-tools.svg?branch=master)](https://travis-ci.org/UCLOrengoGroup/cath-tools) [![Documentation Status](https://readthedocs.org/projects/cath-tools/badge/?version=latest)](https://readthedocs.org/projects/cath-tools/?badge=latest) <iframe src="https://ghbtns.com/github-btn.html?user=UCLOrengoGroup&repo=cath-tools&type=star&count=true" frameborder="0" scrolling="0" width="170px" height="20px" style="vertical-align: middle;"></iframe>

Protein structure comparison tools such as SSAP, as used by the [Orengo Group](https://www.ucl.ac.uk/orengo-group "Orengo Group website") in curating [CATH](http://www.cathdb.info/ "CATH website").

| Resources | | | |
|:-- |:-- |:-- |:-- |
| [**Executable DOWNLOADS**](https://github.com/UCLOrengoGroup/cath-tools/releases/latest "Downloads for the latest CATH Tools release") <br> <sup> (for 64-bit Linux; chmod them to be executable)</sup> | [**Docs**](http://cath-tools.readthedocs.io/en/latest/ "CATH Tools user documentation at Read the Docs") <br> <sup> &nbsp; </sup> | [**Code**](https://github.com/UCLOrengoGroup/cath-tools "CATH Tools GitHub repository") <br> <sup> &nbsp; </sup> | |
| [secmake](http://github.com/UCLOrengoGroup/secmake) <br> <sup>(makes sec files, required by cath-ssap) </sup> | [Builds](https://travis-ci.org/UCLOrengoGroup/cath-tools "CATH Tools Travis-CI builds") <br> <sup> &nbsp; </sup> | [Extras repo](https://github.com/UCLOrengoGroup/cath-tools-supplementary "CATH Tools Supplementary GitHub repository") <br> <sup> &nbsp; </sup> | [bioplib](https://github.com/ACRMGroup/bioplib "Bioplib's GitHub Homepage")<br> <sup>(used to execute the superpositions)</sup> |


### Tools

| | |
|:-- |:-- |
| <img class="tool-thumb" style="border-style: solid; border-width: 1px;" src="https://raw.githubusercontent.com/UCLOrengoGroup/cath-tools/master/docs/tools/img/cath-resolve-hits.example.jpg" width="250" /> | [ **cath-resolve-hits** ]( http://cath-tools.readthedocs.io/en/latest/tools/cath-resolve-hits ) <br> Collapse a list of domain matches to your query sequence(s) down to the non-overlapping subset (ie domain architecture) that maximises the sum of the hits' scores. |
| <img class="tool-thumb" style="border-style: solid; border-width: 1px;" src="https://raw.githubusercontent.com/UCLOrengoGroup/cath-tools/master/docs/tools/img/ssap_alignment.jpg"            width="250" /> | [ **cath-ssap**         ]( http://cath-tools.readthedocs.io/en/latest/tools/cath-ssap         ) <br> Structurally align a pair of proteins. |
| <img class="tool-thumb" style="border-style: solid; border-width: 1px;" src="https://raw.githubusercontent.com/UCLOrengoGroup/cath-tools/master/docs/tools/img/1fi2A00_1j58A01.good.jpg"      width="250" /> | [ **cath-superpose**    ]( http://cath-tools.readthedocs.io/en/latest/tools/cath-superpose    ) <br> Superpose two or more protein structures using an existing alignment. |

### Extra Tools

 * `build-test`          Perform the cath-tools tests (which should all pass, abeit with a few warnings)
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

  * [Tony E Lewis](https://github.com/tonyelewis)             (2011  - ....)
  * Oliver C Redfern                                          (~2003 - 2011)
  * James E Bray, [Ian Sillitoe](https://github.com/sillitoe) (~2000 - 2003)
  * [Andrew C R Martin](https://github.com/AndrewCRMartin)    (considerable edits around 2001)

Feedback
--------

Please tell us about your cath-tools bugs/suggestions [here](https://github.com/UCLOrengoGroup/cath-tools/issues/new).

If you find this software useful, please spread the word and star the [GitHub repo](https://github.com/UCLOrengoGroup/cath-tools). <iframe src="https://ghbtns.com/github-btn.html?user=UCLOrengoGroup&repo=cath-tools&type=star&count=true" frameborder="0" scrolling="0" width="170px" height="20px" style="vertical-align: middle;"></iframe>
