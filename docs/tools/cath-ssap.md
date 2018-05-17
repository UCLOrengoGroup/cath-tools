# cath-ssap

[**Downloads**](https://github.com/UCLOrengoGroup/cath-tools/releases/latest)

The binary for SSAP is called `cath-ssap`. For usage, run `cath-ssap --help` or see [below](#usage).

## Preparing to run SSAP

`cath-ssap` needs to know where to find PDB files. You can tell it using the `--pdb-path` option, which allows you to specify multiple directories in the order they should be searched by separating them with colons (`:`), eg:

`cath-ssap --pdb-path .:/global_pdbs 1cukA01 1bvsA01`

Similarly, you can specify the style of prefix/suffix you use for your PDB files, with the `--pdb-prefix` and `--pdb-suffix` options.

Since these options' values will typically be fairly stable, it can be tedious to specify them every time. Fortunately, the cath-tools allow their all their command-line options to also be specified via environment variables, eg:

~~~~~no-highlight
export CATH_TOOLS_PDB_PATH=.:/global/data/directories/pdb
export CATH_TOOLS_PDB_PREFIX=pdb
export CATH_TOOLS_PDB_SUFFIX=.ent
~~~~~

The defaults values for these parameters are included in the usage information (run `cath-ssap --help` or see [below](#usage)).

## Running SSAP

Once you've set up these environment variables, you can use a command like:

`cath-ssap 1cukA 1bvsA`

This prints a short summary of the resulting scores (see `cath-ssap --scores-help` for format) and writes an alignment file to the current directory (see `cath-ssap --alignment-help` for format).

Once you've aligned structures with `cath-ssap`, you can make better-looking superpositions with [`cath-superpose`](cath-superpose).


## Usage

The current usage information is:

~~~~~no-highlight
Usage: cath-ssap [options] <protein1> <protein2>

Run a SSAP pairwise structural alignment
[algorithm devised by C A Orengo and W R Taylor, see --citation-help]

cath-ssap uses two types of structural comparison:
  1. Fast SSAP: a quick secondary-structure based SSAP alignment
  2. Slow SSAP: residue alignment only

If both structures have more than one SS element, a fast SSAP is run first. If the fast SSAP score isn't good, another fast SSAP is run with looser cutoffs. If the (best) fast SSAP score isn't good, a slow SSAP is run. Only the best of these scores is output. These behaviours can be configured using the parameters below.)

Miscellaneous:
  -h [ --help ]                            Output help message
  -v [ --version ]                         Output version information

Standard SSAP options:
  --debug                                  Output debugging information
  -o [ --outfile ] <file>                  [DEPRECATED] Output scores to <file> rather than to stdout
  --clique-file <file>                     Read clique from <file>
  --domin-file <file>                      Read domin from <file>
  --max-score-to-fast-rerun <score> (=65)  Run a second fast SSAP with looser cutoffs if the first fast SSAP's score falls below <score>
  --max-score-to-slow-rerun <score> (=75)  Perform a slow SSAP if the (best) fast SSAP score falls below <score>
  --slow-ssap-only                         Don't try any fast SSAPs; only use slow SSAP
  --local-ssap-score                       [DEPRECATED] Normalise the SSAP score over the length of the smallest domain rather than the largest
  --all-scores                             [DEPRECATED] Output all SSAP scores from fast and slow runs, not just the highest
  --prot-src-files <set> (=PDB)            Read the protein data from the set of files <set>, of available sets:
                                           PDB, PDB_DSSP, PDB_DSSP_SEC, WOLF_SEC
  --supdir <dir>                           [DEPRECATED] Output a superposition to directory <dir>
  --aligndir <dir> (=".")                  Write alignment to directory <dir>
  --min-score-for-files <score> (=0)       Only output alignment/superposition files if the SSAP score exceeds <score>
  --min-sup-score <score> (=-0.25)         [DEPRECATED] Calculate superposition based on the residue-pairs with scores greater than <score>
  --rasmol-script                          [DEPRECATED] Write a rasmol superposition script to load and colour the superposed structures
  --xmlsup                                 [DEPRECATED] Write a small xml superposition file, from which a larger superposition file can be reconstructed

Conversion between a protein's name and its data files:
  --pdb-path <path> (=.)                   Search for PDB files using the path <path>
  --dssp-path <path> (=.)                  Search for DSSP files using the path <path>
  --wolf-path <path> (=.)                  Search for wolf files using the path <path>
  --sec-path <path> (=.)                   Search for sec files using the path <path>
  --pdb-prefix <pre>                       Prepend the prefix <pre> to a protein's name to form its PDB filename
  --dssp-prefix <pre>                      Prepend the prefix <pre> to a protein's name to form its DSSP filename
  --wolf-prefix <pre>                      Prepend the prefix <pre> to a protein's name to form its wolf filename
  --sec-prefix <pre>                       Prepend the prefix <pre> to a protein's name to form its sec filename
  --pdb-suffix <suf>                       Append the suffix <suf> to a protein's name to form its PDB filename
  --dssp-suffix <suf> (=.dssp)             Append the suffix <suf> to a protein's name to form its DSSP filename
  --wolf-suffix <suf> (=.wolf)             Append the suffix <suf> to a protein's name to form its wolf filename
  --sec-suffix <suf> (=.sec)               Append the suffix <suf> to a protein's name to form its sec filename

Regions:
  --align-regions <regions>                Handle region(s) <regions> as the alignment part of the structure.
                                           May be specified multiple times, in correspondence with the structures.
                                           Format is: D[5inwB02]251-348:B,408-416A:B
                                           (Put <regions> in quotes to prevent the square brackets confusing your shell ("No match"))

Detailed help:
  --alignment-help                         Help on alignment format
  --citation-help                          Help on SSAP authorship & how to cite it
  --scores-help                            Help on scores format

Please tell us your cath-tools bugs/suggestions : https://github.com/UCLOrengoGroup/cath-tools/issues/new
~~~~~

## DSSP, WOLF and sec

We recommend you run `cath-ssap` from PDB files only. However it's also capable of reading its data from other combinations of input files:

* PDB and DSSP
* PDB, DSSP and SEC
* WOLF and SEC

The reason for this is that `cath-ssap` used to require a WOLF/DSSP file for per-residue secondary structure information and a SEC file for information on each secondary structure; it can now calculate all the information it needs itself. For this reason, the PDB-only mode typically runs slightly slower, however it is simpler to use and avoids various problems with those formats. Our benchmarking (comparing ROC curves to assess ability to discriminate whether pairs are homologous) shows that the different modes can give slightly different results in some cases but perform equally well overall.

To generate DSSP and sec files from PDB files, you can use:

* [`dssp`](http://swift.cmbi.ru.nl/gv/dssp), the CMBI tool for generating DSSP files from PDB files
* [`secmake`](http://github.com/UCLOrengoGroup/secmake), a tool for generating sec files from PDB + DSSP files

Once you've prepared these files, you need to tell `cath-ssap` where to find them. This can be done in the same way as for PDBs (see [above](#preparing-to-run-ssap)) using the environment variables `CATH_TOOLS_DSSP_PATH` and `CATH_TOOLS_SEC_PATH` (and for non-standard prefixes/suffixes `CATH_TOOLS_DSSP_PREFIX`, `CATH_TOOLS_SEC_PREFIX`, `CATH_TOOLS_DSSP_SUFFIX` and `CATH_TOOLS_SEC_SUFFIX`).

## Feedback

Please tell us about your cath-tools bugs/suggestions [here](https://github.com/UCLOrengoGroup/cath-tools/issues/new).
