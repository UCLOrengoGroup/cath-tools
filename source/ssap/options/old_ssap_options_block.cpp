/// \file
/// \brief The old_ssap_options_block class definitions

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

#include "old_ssap_options_block.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include "alignment/common_residue_selection_policy/common_residue_select_min_score_policy.hpp"
#include "common/clone/make_uptr_clone.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "exception/out_of_range_exception.hpp"
#include "structure/protein/protein_source_file_set/protein_source_file_set.hpp"

#include <iostream>

using namespace boost::algorithm;
using namespace boost::program_options;
using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::opts;
using namespace cath::sup;
using namespace std;

using boost::algorithm::join;
using boost::filesystem::path;
using boost::lexical_cast;
using boost::none;

constexpr bool               old_ssap_options_block::DEF_BOOL;
constexpr protein_file_combn old_ssap_options_block::DEF_PROT_SRCS;
constexpr double             old_ssap_options_block::DEF_REFAST;
constexpr double             old_ssap_options_block::DEF_RESLOW;
constexpr double             old_ssap_options_block::DEF_FILE_SC;
constexpr double             old_ssap_options_block::DEF_SUP;

const string old_ssap_options_block::PO_NAME                 = { "name"                    }; ///< The option name for the names option

const string old_ssap_options_block::PO_DEBUG                = { "debug"                   }; ///< The option name for the debug option
const string old_ssap_options_block::PO_OUT_FILE             = { "outfile"                 }; ///< The option name for the output_filename option

const string old_ssap_options_block::PO_CLIQUE_FILE          = { "clique-file"             }; ///< The option name for the clique_file option
const string old_ssap_options_block::PO_DOMIN_FILE           = { "domin-file"              }; ///< The option name for the domin_file option

const string old_ssap_options_block::PO_MAX_SCORE_TO_REFAST  = { "max-score-to-fast-rerun" }; ///< The option name for the max_score_to_fast_ssap_rerun option
const string old_ssap_options_block::PO_MAX_SCORE_TO_RESLOW  = { "max-score-to-slow-rerun" }; ///< The option name for the max_score_to_slow_ssap_rerun option
const string old_ssap_options_block::PO_SLOW_SSAP_ONLY       = { "slow-ssap-only"          }; ///< The option name for the slow_ssap_only option

const string old_ssap_options_block::PO_LOC_SSAP_SCORE       = { "local-ssap-score"        }; ///< The option name for the use_local_ssap_score option
const string old_ssap_options_block::PO_ALL_SCORES           = { "all-scores"              }; ///< The option name for the write_all_scores option
const string old_ssap_options_block::PO_PROTEIN_SOURCE_FILES = { "prot-src-files"          }; ///< The option name for the protein_source_files option

const string old_ssap_options_block::PO_SUPN_DIR             = { "supdir"                  }; ///< The option name for the superposition_dir option
const string old_ssap_options_block::PO_ALIGN_DIR            = { "aligndir"                }; ///< The option name for the alignment_dir option
const string old_ssap_options_block::PO_MIN_OUT_SCORE        = { "min-score-for-files"     }; ///< The option name for the min_score_for_writing_files option
const string old_ssap_options_block::PO_MIN_SUP_SCORE        = { "min-sup-score"           }; ///< The option name for the min_score_for_superposition option
const string old_ssap_options_block::PO_RASMOL_SCRIPT        = { "rasmol-script"           }; ///< The option name for the write_rasmol_script option
const string old_ssap_options_block::PO_XML_SUP              = { "xmlsup"                  }; ///< The option name for write_xml_sup option

/// \brief The single-character for the output file option
constexpr char old_ssap_options_block::PO_CHAR_OUT_FILE;

/// \brief A standard do_clone() method to act as a virtual copy-ctor
///
/// This is a concrete definition of a virtual method that's pure in options_block
unique_ptr<options_block> old_ssap_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Provide a name for the options block, as used in the options description text
///
/// This is a concrete definition of a virtual method that's pure in options_block
string old_ssap_options_block::do_get_block_name() const {
	return "Standard SSAP options";
}

/// \brief Add the block's non-hidden options to the provided options_description
///
/// This is a concrete definition of a virtual method that's pure in options_block
void old_ssap_options_block::do_add_visible_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                   ) {
	const string dir_varname  { "<dir>"   };
	const string file_varname { "<file>"  };
	const string score_varname{ "<score>" };
	const string set_varname  { "<set>"   };

	const auto write_rasmol_script_notifier = [&] (const bool &x) { set_write_rasmol_script( x ? sup_pdbs_script_policy::WRITE_RASMOL_SCRIPT : sup_pdbs_script_policy::LEAVE_RAW_PDBS ); };

	const string PO_OUT_FILE_W_CHAR = PO_OUT_FILE + ',' + PO_CHAR_OUT_FILE;

	const string protein_file_combns_str = join( get_ssap_ready_protein_file_combn_strings(), ", " );
	arg_desc.add_options()
		( PO_DEBUG.c_str(),                bool_switch              ( &debug                        )                           ->default_value(DEF_BOOL      ),   "Output debugging information"                                                                                          )
		( PO_OUT_FILE_W_CHAR.c_str(),      value<path>              ( &output_filename              )->value_name(file_varname ),                                ( "[DEPRECATED] Output scores to " + file_varname + " rather than to stdout" ).c_str()                                    )

		( PO_CLIQUE_FILE.c_str(),          value<path>              ( &clique_file                  )->value_name(file_varname ),                                ( "Read clique from " + file_varname ).c_str()                                                                            )
		( PO_DOMIN_FILE.c_str(),           value<path>              ( &domin_file                   )->value_name(file_varname ),                                ( "Read domin from "  + file_varname ).c_str()                                                                            )

		( PO_MAX_SCORE_TO_REFAST.c_str(),  value<double>            ( &max_score_to_fast_ssap_rerun )->value_name(score_varname)->default_value(DEF_REFAST    ), ( "Run a second fast SSAP with looser cutoffs if the first fast SSAP's score falls below " + score_varname ).c_str()      )
		( PO_MAX_SCORE_TO_RESLOW.c_str(),  value<double>            ( &max_score_to_slow_ssap_rerun )->value_name(score_varname)->default_value(DEF_RESLOW    ), ( "Perform a slow SSAP if the (best) fast SSAP score falls below " + score_varname ).c_str()                              )
		( PO_SLOW_SSAP_ONLY.c_str(),       bool_switch              ( &slow_ssap_only               )                           ->default_value(DEF_BOOL      ),   "Don't try any fast SSAPs; only use slow SSAP"                                                                          )

		( PO_LOC_SSAP_SCORE.c_str(),       bool_switch              ( &use_local_ssap_score         )                           ->default_value(DEF_BOOL      ),   "[DEPRECATED] Normalise the SSAP score over the length of the smallest domain rather than the largest"                  )
		( PO_ALL_SCORES.c_str(),           bool_switch              ( &write_all_scores             )                           ->default_value(DEF_BOOL      ),   "[DEPRECATED] Output all SSAP scores from fast and slow runs, not just the highest"                                     )
		( PO_PROTEIN_SOURCE_FILES.c_str(), value<protein_file_combn>( &protein_source_files         )->value_name( set_varname )->default_value(DEF_PROT_SRCS ), ( "Read the protein data from the set of files "+set_varname+", of available sets:\n" + protein_file_combns_str ).c_str() )

		( PO_SUPN_DIR.c_str(),             value<path>              ( &superposition_dir            )->value_name( dir_varname ),                                ( "[DEPRECATED] Output a superposition to directory " + dir_varname ).c_str()                                             )
		( PO_ALIGN_DIR.c_str(),            value<path>              ( &alignment_dir                )->value_name( dir_varname )->default_value( path(".")    ), ( "Write alignment to directory " + dir_varname ).c_str()                                                                 )
		( PO_MIN_OUT_SCORE.c_str(),        value<double>            ( &min_score_for_writing_files  )->value_name(score_varname)->default_value(DEF_FILE_SC   ), ( "Only output alignment/superposition files if the SSAP score exceeds " + score_varname ).c_str()                        )
		( PO_MIN_SUP_SCORE.c_str(),        value<double>            ( &min_score_for_superposition  )->value_name(score_varname)->default_value(DEF_SUP       ), ( "[DEPRECATED] Calculate superposition based on the residue-pairs with scores greater than " + score_varname ).c_str()   )
		( PO_RASMOL_SCRIPT.c_str(),        bool_switch()->notifier  ( write_rasmol_script_notifier  ),                                                             "[DEPRECATED] Write a rasmol superposition script to load and colour the superposed structures"                         )
		( PO_XML_SUP.c_str(),              bool_switch              ( &write_xml_sup                )                           ->default_value(DEF_BOOL      ),   "[DEPRECATED] Write a small xml superposition file, from which a larger superposition file can be reconstructed"        );
}

/// \brief Add any hidden options to the provided options_description
void old_ssap_options_block::do_add_hidden_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                  ) {
	arg_desc.add_options()
		( PO_NAME.c_str(), value<str_vec>( &names ), "Structure names" );
}

/// \brief Identify any conflicts that make the currently stored options invalid
///
/// This is a concrete definition of a virtual method that's pure in options_block
///
/// These methods should only catch absolute conflicts that could never be acceptable.
/// In particular, this permits no names being specified because the client code may
/// want to accept that if the user has requested help.
///
/// Current checks:
///  * Always accept if no names have been specified
///  * Otherwise, reject if there aren't exactly two names
///  * Reject if the min_score_for_superposition is below common_residue_select_min_score_policy::MIN_CUTOFF
///  * Reject if a specified clique file isn't a valid input file
///  * Reject if a specified domin file isn't a valid input file
///  * Reject if a specified superposition output directory isn't a valid output directory
///  * Reject if the alignment output directory isn't a valid output directory
///  * Otherwise accept
///
/// \returns A string describing the conflict in the options or an empty string if there's none
str_opt old_ssap_options_block::do_invalid_string(const variables_map &/*arg_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                  ) const {
	// Always accept if no names have been specified
	if ( names.empty() ) {
		return none;
	}

	// Otherwise, reject if there aren't exactly two names
	if ( names.size() != 2 ) {
		return "Requires exactly two protein structures (but " + lexical_cast<string>(names.size()) + (names.size() == 1 ? " has" : " have") + " been specified)";
	}

	// Reject if the min_score_for_superposition is below common_residue_select_min_score_policy::MIN_CUTOFF
	if (get_min_score_for_superposition() < common_residue_select_min_score_policy::MIN_CUTOFF) {
		return "The sup-score value must be at least as large as " + lexical_cast<string>(common_residue_select_min_score_policy::MIN_CUTOFF);
	}

	// Reject if a specified clique file isn't a valid input file
	if ( has_clique_file( *this ) && ! is_acceptable_input_file( get_clique_file( *this ) ) ) {
		return "Clique file " + get_clique_file( *this ).string() + " is not a valid input file";
	}

	// Reject if a specified domin file isn't a valid input file
	if ( has_domin_file( *this ) && ! is_acceptable_input_file( get_domin_file( *this ) ) ) {
		return "Domin file " + get_domin_file( *this ).string() + " is not a valid input file";
	}

	// Reject if a specified superposition output directory isn't a valid output directory
	if ( has_superposition_dir( *this ) && ! is_acceptable_output_dir( get_superposition_dir( *this ) ) ) {
		return "Superposition directory " + get_superposition_dir( *this ).string() + " is not a valid output directory";
	}

	// Reject if the alignment output directory isn't a valid output directory
	if (!is_acceptable_output_dir(get_alignment_dir())) {
		return "Alignment directory " + get_alignment_dir().string() + "";
	}

	// Otherwise accept
	return none;
}

/// \brief Return all options names for this block
str_vec old_ssap_options_block::do_get_all_options_names() const {
	return {
		old_ssap_options_block::PO_NAME,
		old_ssap_options_block::PO_DEBUG,
		old_ssap_options_block::PO_OUT_FILE,
		old_ssap_options_block::PO_CLIQUE_FILE,
		old_ssap_options_block::PO_DOMIN_FILE,
		old_ssap_options_block::PO_MAX_SCORE_TO_REFAST,
		old_ssap_options_block::PO_MAX_SCORE_TO_RESLOW,
		old_ssap_options_block::PO_SLOW_SSAP_ONLY,
		old_ssap_options_block::PO_LOC_SSAP_SCORE,
		old_ssap_options_block::PO_ALL_SCORES,
		old_ssap_options_block::PO_PROTEIN_SOURCE_FILES,
		old_ssap_options_block::PO_SUPN_DIR,
		old_ssap_options_block::PO_ALIGN_DIR,
		old_ssap_options_block::PO_MIN_OUT_SCORE,
		old_ssap_options_block::PO_MIN_SUP_SCORE,
		old_ssap_options_block::PO_RASMOL_SCRIPT,
		old_ssap_options_block::PO_XML_SUP,
	};
}

/// \brief Whether any protein names have been specified
bool old_ssap_options_block::protein_names_specified() const {
	return ! names.empty();
}

/// \brief Getter for the name of the first protein structure to be compared
string old_ssap_options_block::get_protein_name_a() const {
	if (names.size() < 2) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("Cannot get protein name a because two names have not (yet) been specified"));
	}
	return names[ 0 ];
}

/// \brief Getter for the name of the second protein structure to be compared
string old_ssap_options_block::get_protein_name_b() const {
	if (names.size() < 2) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("Cannot get protein name a because two names have not (yet) been specified"));
	}
	return names[ 1 ];
}

/// \brief Getter for debug
bool old_ssap_options_block::get_debug() const {
	return debug;
}

/// \brief Getter for whether an output_filename has been specified
bool old_ssap_options_block::get_output_to_file() const {
	return !output_filename.empty();
}

/// \brief Getter for output_filename
path old_ssap_options_block::get_output_filename() const {
	return output_filename;
}

/// \brief TODOCUMENT
path_opt old_ssap_options_block::get_opt_clique_file() const {
	return ( ! clique_file.empty() ) ? path_opt( clique_file ) : none;
}

/// \brief TODOCUMENT
path_opt old_ssap_options_block::get_opt_domin_file() const {
	return ( ! domin_file.empty() ) ? path_opt( domin_file ) : none;
}

/// \brief Getter for max_score_considered_distant
double old_ssap_options_block::get_max_score_to_fast_ssap_rerun() const {
	return max_score_to_fast_ssap_rerun;
}

/// \brief Getter for max_score_for_slow_ssap_rerun
double old_ssap_options_block::get_max_score_to_slow_ssap_rerun() const {
	return max_score_to_slow_ssap_rerun;
}

/// \brief Getter for no_fast_ssap
bool old_ssap_options_block::get_slow_ssap_only() const {
	return slow_ssap_only;
}

/// \brief Getter for use_local_score
bool old_ssap_options_block::get_use_local_ssap_score() const {
	return use_local_ssap_score;
}

/// \brief Getter for write_all_scores
bool old_ssap_options_block::get_write_all_scores() const {
	return write_all_scores;
}

/// \brief Getter for protein_source_files
unique_ptr<const protein_source_file_set> old_ssap_options_block::get_protein_source_files() const {
	return get_protein_source_file_set( protein_source_files );
}

/// \brief TODOCUMENT
path_opt old_ssap_options_block::get_opt_superposition_dir() const {
	return ( ! superposition_dir.empty() ) ? path_opt( superposition_dir ) : none;
}

/// \brief Getter for alignment_dir
path old_ssap_options_block::get_alignment_dir() const {
	return alignment_dir;
}

/// \brief Getter for min_score_for_writing_files
double old_ssap_options_block::get_min_score_for_writing_files() const {
	return min_score_for_writing_files;
}

/// \brief Getter for min_score_for_superposition
double old_ssap_options_block::get_min_score_for_superposition() const {
	return min_score_for_superposition;
}

/// \brief Getter for write_script
sup_pdbs_script_policy old_ssap_options_block::get_write_rasmol_script() const {
	return write_rasmol_script;
}

/// \brief Getter for write_xml_sup
bool old_ssap_options_block::get_write_xml_sup() const {
	return write_xml_sup;
}

/// \brief Setter for write_script
old_ssap_options_block & old_ssap_options_block::set_write_rasmol_script(const sup_pdbs_script_policy &arg_write_rasmol_script ///< The new policy for writing a script for superposition PDBs
                                                                         ) {
	write_rasmol_script = arg_write_rasmol_script;
	return *this;
}

/// \brief Getter for whether a clique_file has been specified
bool cath::opts::has_clique_file(const old_ssap_options_block &arg_old_ssap_options_block ///< TODOCUMENT
                                 ) {
	return static_cast<bool>( arg_old_ssap_options_block.get_opt_clique_file() );
}

/// \brief Getter for clique_file
path cath::opts::get_clique_file(const old_ssap_options_block &arg_old_ssap_options_block ///< TODOCUMENT
                                 ) {
	return *arg_old_ssap_options_block.get_opt_clique_file();
}

/// \brief Getter for whether a domin_file has been specified
bool cath::opts::has_domin_file(const old_ssap_options_block &arg_old_ssap_options_block ///< TODOCUMENT
                                ) {
	return static_cast<bool>( arg_old_ssap_options_block.get_opt_domin_file() );
}

/// \brief Getter for domin_file
path cath::opts::get_domin_file(const old_ssap_options_block &arg_old_ssap_options_block ///< TODOCUMENT
                                ) {
	return *arg_old_ssap_options_block.get_opt_domin_file();
}

/// \brief Getter for whether a superposition_dir has been specified
bool cath::opts::has_superposition_dir(const old_ssap_options_block &arg_old_ssap_options_block ///< TODOCUMENT
                                       ) {
	return static_cast<bool>( arg_old_ssap_options_block.get_opt_superposition_dir() );
}

/// \brief Getter for superposition_dir
path cath::opts::get_superposition_dir(const old_ssap_options_block &arg_old_ssap_options_block ///< TODOCUMENT
                                       ) {
	return *arg_old_ssap_options_block.get_opt_superposition_dir();
}
