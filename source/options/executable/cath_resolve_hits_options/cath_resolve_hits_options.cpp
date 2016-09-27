/// \file
/// \brief The cath_resolve_hits_options class definitions

/// \copyright
/// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
/// Copyright (C) 2011, Orengo Group, University College London
///
/// This program is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "cath_resolve_hits_options.h"

#include <boost/program_options.hpp>
#include <boost/shared_array.hpp>

#include "acquirer/alignment_acquirer/alignment_acquirer.h"
#include "acquirer/pdbs_acquirer/file_list_pdbs_acquirer.h"
#include "acquirer/pdbs_acquirer/istream_pdbs_acquirer.h"
#include "acquirer/selection_policy_acquirer/selection_policy_acquirer.h"
#include "acquirer/superposition_acquirer/align_based_superposition_acquirer.h"
#include "alignment/alignment.h"
#include "alignment/common_atom_selection_policy/common_atom_select_ca_policy.h"
#include "alignment/common_residue_selection_policy/common_residue_select_all_policy.h"
#include "alignment/common_residue_selection_policy/common_residue_select_best_score_percent_policy.h"
#include "common/argc_argv_faker.h"
#include "common/type_aliases.h"
#include "exception/invalid_argument_exception.h"
#include "exception/not_implemented_exception.h"
#include "exception/runtime_error_exception.h"
#include "file/pdb/pdb.h"
#include "file/pdb/pdb_atom.h"
#include "file/pdb/pdb_residue.h"
#include "outputter/alignment_outputter/alignment_outputter.h"
#include "outputter/alignment_outputter/alignment_outputter_list.h"
#include "outputter/superposition_outputter/superposition_outputter.h"
#include "outputter/superposition_outputter/superposition_outputter_list.h"
#include "superposition/superposition_context.h"

#include <iostream>

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::opts;

using boost::lexical_cast;
using boost::none;
using boost::program_options::options_description;
using boost::program_options::positional_options_description;
using boost::ptr_vector;
using std::string;

/// The name of the program that uses this executable_options
const string cath_resolve_hits_options::PROGRAM_NAME("cath-resolve-hits");

/// \brief Get the name of the program that uses this executable_options
string cath_resolve_hits_options::do_get_program_name() const {
	return PROGRAM_NAME;
}

/// \brief TODOCUMENT
positional_options_description cath_resolve_hits_options::get_positional_options() {
	positional_options_description positionals;
	positionals.add( resolve_hits_input_options_block::PO_INPUT_FILE_OR_STDIN.c_str(), 1 );
	return positionals;
}

/// \brief Review all specified options and return a string containing any errors or a help string
///        (possibly using a description of all visible options)
///
/// This is a concrete definition of a virtual method that's pure in executable_options
///
/// This should only be called by executable_options, as the last step of the parse_options()
/// method, after all real parsing has completed.
///
/// \pre The options should have been parsed
///
/// \returns Any error/help string arising from the newly specified options
///          or an empty string if there aren't any
opt_str cath_resolve_hits_options::do_get_error_or_help_string() const {
	// // Grab the objects from the options blocks
	// const ptr_vector<alignment_acquirer> alignment_acquirers = the_alignment_input_options_block.get_alignment_acquirers();
	// const ptr_vector<pdbs_acquirer>      pdb_acquirers       = the_pdb_input_options_block.get_pdbs_acquirers();
	// const display_spec                   the_display_spec    = the_display_options_block.get_display_spec();
	// const alignment_outputter_list       aln_outputters      = the_alignment_output_options_block.get_alignment_outputters( the_display_spec );
	// const superposition_outputter_list   sup_outputters      = the_superposition_output_options_block.get_superposition_outputters( the_display_spec );

	// // If there are no objects then no options were specified so just output the standard usage error string
	// if (alignment_acquirers.empty() && pdb_acquirers.empty() && sup_outputters.empty()) {
	// 	return ""s;
	// }

	// // Check that there is exactly one source of alignment
	// if (alignment_acquirers.size() != 1) {
	// 	return "Please specify one source of an alignment or superposition ("
	// 	       + lexical_cast<string>(alignment_acquirers.size())
	// 	       + " specified)";
	// }

	// // Check that there is exactly one source of PDBs
	// if (pdb_acquirers.size() != 1) {
	// 	return "Please specify one source of PDBs ("
	// 	       + lexical_cast<string>(pdb_acquirers.size())
	// 	       + " specified)";
	// }

	// const variables_map &vm = get_variables_map();

	// const bool using_display_spec =    any_alignment_outputters_involve_display_spec    ( aln_outputters )
	//                                 || any_superposition_outputters_involve_display_spec( sup_outputters );
	// if ( the_display_options_block.has_specified_options( vm ) && ! using_display_spec ) {
	// 	return "Cannot specify display options because no display is being used in superposition/alignment outputters";
	// }

	return none;
}

/// \brief Get a string to prepend to the standard help
string cath_resolve_hits_options::do_get_help_prefix_string() const {
	return "Usage: " + PROGRAM_NAME + R"( [options] input_file

)" + get_overview_string() + R"(

The input data may contain unsorted hits for different query protein sequences.

However, if your input data is already grouped by query protein sequence, then
specify the --)" + resolve_hits_input_options_block::PO_INPUT_HITS_ARE_GROUPED + R"( flag for faster runs that use less memory.)";
}

/// \brief Get a string to append to the standard help
string cath_resolve_hits_options::do_get_help_suffix_string() const {
	// Consider adding usage examples:
	//
	// Some vague notes from various places
	//  * cath-resolve_hits --)" + resolve_hits_input_options_block::PO_DOMTBLOUT + R"( my
	//  * cat pdb1 end_file pdb2 end_file pdb3 | cath-resolve_hits --pdbs-from-stdin --sup-to-stdout --res-name-align
	//      (Superpose the structures from stdin based on matching residue names and then write them to stdout [common Genome3D use case]))";

	return  R"(

Format for "raw" input data
---------------------------

One hit per line, using the following space-separated fields:

 1. query_protein_id : An identifier for the query protein sequence
 2. match_id         : An identifier for the match
 3. score            : A (strictly positive) score indicating how good it would be to have that hit in the final result
 4. starts_stops     : The starts/stops on the query sequence, given in the format: 37-124,239-331

Example lines:

qyikaz 1mkfA01/12-210-i5_1,2.9e-20 2983.29780221221 3-103
qyikaz 1mkfA01/12-210-i5_2,4.9e-19 3510.41568607646 101-224
qyikaz 1mkfA01/12-210-i5_3,7e-25 3552.10980383852 825-928
qyikaz 1mkfA01/12-210-i5_4,3.5e-15 2470.04912752062 953-1053
)";
}

/// \brief Get an overview of the job that these options are for
///
/// This can be used in the --help and --version outputs
string cath_resolve_hits_options::do_get_overview_string() const {
	return R"(Collapse a list of domain matches to your query sequence(s) down to the
non-overlapping subset (ie domain architecture) that maximises the sum of the
hits' scores.)";
}

/// TODOCUMENT
// void cath_resolve_hits_options::check_ok_to_use() const {
// 	if ( get_error_or_help_string() ) {
// 		BOOST_THROW_EXCEPTION(invalid_argument_exception("Attempt to use invalid cath_resolve_hits_options"));
// 	}
// }


/// \brief TODOCUMENT
cath_resolve_hits_options::cath_resolve_hits_options() {
	super::add_options_block( the_resolve_hits_input_options_block );
}

/// TODOCUMENT
const resolve_hits_input_spec & cath_resolve_hits_options::get_resolve_hits_input_spec() const {
	return the_resolve_hits_input_options_block.get_resolve_hits_input_spec();
}
