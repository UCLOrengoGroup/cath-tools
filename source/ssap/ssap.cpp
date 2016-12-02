// \copyright
// CATH Tools - Protein structure comparison tools such as SSAP and SNAP
// Copyright (C) 1989, Orengo Group, University College London
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

/// \file
/// \brief Definitions of core SSAP functions

/// \mainpage
///
/// About
/// =====
///
/// Algorithm devised by Christine A Orengo and William R Taylor
///
/// Please cite: "Protein Structure Alignment", Taylor and Orengo [1989]
///              Journal of Molecular Biology 208, 1-22
///              PMID: 2769748
///
/// Many people have contributed to this code, most notably:
///   * Tony E Lewis               (  2011 - ....)
///   * Oliver C Redfern           (~ 2003 - 2011)
///   * James E Bray, Ian Sillitoe (~ 2000 - 2003)
///   * Andrew C R Martin          (considerable edits around 2001)
///
/// Code Design Principles
/// ======================
///
/// Disclaimer
/// ----------
///
/// For backward consistency, this project is not being written from scratch (which would be simpler)
/// but is instead evolving out of the old SSAP code base. This is necessary but it means that the code
/// standards are inconsistent between the older code (mostly in ssap/ssap.cpp) and the newer code.
/// It also means the new code is often constrained to mirror the old code in various ways.
/// Over time, the constraints are being removed and the inconsistencies eradicated but it's unclear
/// whether enough time will be devoted to this project to warrant removing this disclaimer.
///
/// Principles
/// ----------
///
/// An overview of some of the code's common standards/acronyms etc (please web-search and read up on
/// stuff you need to understand) :
/// * The code aims towards the SOLID design principles
///   (Single responsibility, Open-closed, Liskov substitution, Interface segregation and Dependency Inversion).
/// * An ABC is an Abstract Base Class, which is a bit like an interface in several other languages.
///   It's a class that contains at least one pure-virtual method: a method that can be overridden by
///   derived classes ("virtual") but that provides no definition of the method itself ("pure").
///   This means the class itself cannot be instantiated but its derived types can.
/// * An NVI is a Non-Virtual Interface. This refers to the widely acknowledged (C++) design principle
///   that a class's public interface should not contain virtual methods (basically because the
///   virtual methods perform one role (providing an interface to derived classes) and the public interface
///   performs another (providing an interface to clients) so these should be kept separate.
///   So a common design is to have a (private) pure-virtual method that's accessed by a simple, NVI pass-through
///   (public, non-virtual) method
/// * Angles are generally in radians unless explicitly stated otherwise
/// * In principle, arrays should always be indexed from 0 (rather than from an offset of 1) unless explicitly stated otherwise
///   although this goal may not have been achieved yet.
///
/// Glossary
/// --------
///
///  * ctor is short for constructor
///
/// Functionality Possible with Different Input File Formats
/// ========================================================
///
/// Quite a bit of complexity arises from different functionality being possible or not depending on the format(s)
/// from which the structures are been read. In the long run, this project's user documentation should provide comprehensive
/// explanations that make it very clear which functions are available on structures read from which formats.
///
/// This is not that; it's an incomplete and possibly partially inaccurate place to store notes, so that it can evolve
/// towards something more complete and useful.
///
/// The basic idea is: to get a particular function the structure must be read from at least one of the files that offers it.
///
/// There are some rules about the combinations of files that are permitted:
///  * A structure requires residue-level data and so cannot be read from a sec file alone
///  * A structure cannot be read from both a DSSP file and a Wolf file (because these roles are very similar)
///
/// Meanings of table entries:
///  * Yes : Currently implemented
///  * NYI : Not Yet Implemented (ie the plan is to implement this but it isn't there at the moment)
///  * No  : Not possible or no plans to implement
///
/// Function                     | PDB    | DSSP    | Wolf | Sec
/// ---------------------------- | ------ | ------- | ---- | ---
/// Residue-level data           | Yes    | Yes     | Yes  | No
/// PHI/PSI angles               | Yes    | Yes     | Yes  | No
/// Secondary structure          | No     | Yes     | Yes  | Yes
/// CB atom                      | Yes    | No      | Yes  | No
/// Residue Orientation          | Yes    | No      | Yes  | No
/// Non CA/CB atoms (ie N/C/O)   | NYI    | No      | No   | No
/// Handling incomplete residues | Yes    | Partial | No   | No
/// Robust residue indexing      | Yes    | Yes     | Yes  | No

/// \todo Remove every snprintf() call
///       (have at least already changed them from sprintf() calls!)

#include "ssap.hpp"

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/range/algorithm/stable_sort.hpp>

#include "alignment/alignment_coord_extractor.hpp"
#include "alignment/common_residue_selection_policy/common_residue_select_min_score_policy.hpp"
#include "alignment/dyn_prog_align/dyn_prog_score_source/entry_querier_dyn_prog_score_source.hpp"
#include "alignment/dyn_prog_align/dyn_prog_score_source/mask_dyn_prog_score_source.hpp"
#include "alignment/dyn_prog_align/dyn_prog_score_source/old_matrix_dyn_prog_score_source.hpp"
#include "alignment/dyn_prog_align/ssap_code_dyn_prog_aligner.hpp"
#include "alignment/gap/gap_penalty.hpp"
#include "alignment/io/alignment_io.hpp"
#include "alignment/pair_alignment.hpp"
#include "common/container/vector_of_vector.hpp"
#include "common/difference.hpp"
#include "common/file/open_fstream.hpp"
#include "common/size_t_literal.hpp"
#include "common/temp_check_offset_1.hpp"
#include "common/type_aliases.hpp"
#include "exception/invalid_argument_exception.hpp"
#include "exception/not_implemented_exception.hpp"
#include "exception/out_of_range_exception.hpp"
#include "ssap/clique.hpp"
#include "ssap/options/cath_ssap_options.hpp"
#include "ssap/options/old_ssap_options_block.hpp"
#include "ssap/selected_pair.hpp"
#include "ssap/ssap_scores.hpp"
#include "ssap/windowed_matrix.hpp"
#include "structure/entry_querier/residue_querier.hpp"
#include "structure/entry_querier/sec_struc_querier.hpp"
#include "structure/geometry/coord.hpp"
#include "structure/geometry/coord_list.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/protein_io.hpp"
#include "structure/protein/protein_source_file_set/protein_source_file_set.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"
#include "superposition/io/superposition_io.hpp"
#include "superposition/superposition.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>

using namespace boost::algorithm;
using namespace boost::filesystem;
using namespace boost::log;
using namespace boost::log::trivial;
using namespace cath;
using namespace cath::align;
using namespace cath::align::gap;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::sup;
using namespace cath::opts;

using boost::irange;
using boost::lexical_cast;
using boost::none;
using boost::numeric_cast;
using boost::range::stable_sort;
using std::abs;
using std::boolalpha;
using std::deque;
using std::endl;
using std::fill_n;
using std::fixed;
using std::make_pair;
using std::max;
using std::min;
using std::ofstream;
using std::ostream;
using std::pair;
using std::setprecision;
using std::string;
using std::vector;

/// \brief The number of top-scoring residue pairs to select
constexpr size_t     NUM_SELECTIONS_TO_SAVE   =  20;

/// \brief The fixed length of the string into which SSAP output lines are written
///
/// \todo Make those output lines strings and hence eradicate the need for this
constexpr size_t     SSAP_LINE_LENGTH         = 200;

/// \brief The minimum score that a lower matrix dynamic programming must achieve before its resulting alignment's scores
///        get added to the upper matrix
///
/// Note: It is not yet quite clear why the normalisation factor (that's used to normalise the total alignment score)
///       is as it is.
constexpr score_type MIN_LOWER_MAT_RES_SCORE  =  10;

constexpr size_t     SEC_STRUC_PLANAR_W_ANGLE =  10;
constexpr size_t     SEC_STRUC_PLANAR_A_ANGLE =  60;
constexpr size_t     SEC_STRUC_PLANAR_B_ANGLE =   6;
constexpr size_t     SEC_STRUC_PLANAR_C_ANGLE =  10;

// \todo Put these matrices in classes. This will help with:
//        - decent memory management
//        - passing the matrices to subroutines rather than using them globally
//        - allocating the required amount of memory rather than just using 5000x5000 ints for each!
//       Also, I suspect that quite a bit of complication throughout the file is just indexing these matrices,
//       which should be encapsulated.

/// \brief Matrix of upper scores
score_vec_of_vec global_upper_score_matrix;

/// \brief Matrix to mask out comparisons that should be skipped whilst performing upper-matrix residue comparisons
bool_vec_of_vec  global_upper_res_mask_matrix;

/// \brief Matrix to mask out comparisons that should be skipped whilst performing upper-matrix, secondary-structure comparisons
bool_vec_of_vec  global_upper_ss_mask_matrix;

/// \brief Matrix to mask out comparisons that should be skipped whilst performing lower-matrix (residue or secondary structure) comparisons
bool_vec_of_vec  global_lower_mask_matrix;

size_size_pair_vec global_selections;              ///< Selected region within matrix

size_t             global_num_selections  =     0; ///< The number of selected top-scoring residue pairs
size_t             global_window          =     0; ///< The size of the window to
size_t             global_window_add      =    70; ///< The amount that should be added to the difference in lengths to calculate window size
size_t             global_res_sim_cutoff  =   150; ///<

ptrdiff_t          global_run_counter     =     0; ///<

score_type         global_gap_penalty     =    50; ///< The gap penalty to be used in dynamic programming

bool               global_debug           = false; ///< Whether to output debug messages
bool               global_align_pass      = false; ///< Whether the pass is a later, refining alignment pass
bool               global_supaln          =  true; ///<
bool               global_doing_fast_ssap =  true; ///< Whether currently performing a fast SSAP
bool               global_res_score       = false; ///<

double             global_frac_selected   =   0.0; ///<

double             global_score_run1      =   0.0; ///<
double             global_score_run2      =   0.0; ///<
double             global_ssap_score1     =   0.0; ///<
double             global_ssap_score2     =   0.0; ///<

char               global_ssap_line1[SSAP_LINE_LENGTH]; ///<
char               global_ssap_line2[SSAP_LINE_LENGTH]; ///<

/// \brief Reset all the global variable that are used by SSAP
///
/// This is only a temporary solution because the long-term solution should be to eradicate these global variables.
///
/// However it is important to be able to test some of the SSAP functionality in the meantime so it is
/// essential to be able to ensure these variables are reset before each test
///
/// For now, this is just called in the constructor of global_test_constants so test suites wishing to test
/// behaviour of code that uses these globals should ensure that its test suite fixture inherits from global_test_constants.
void cath::reset_ssap_global_variables() {
	global_upper_score_matrix.assign   ( 0, 0, 0     );
	global_upper_res_mask_matrix.assign( 0, 0, false );
	global_upper_ss_mask_matrix.assign ( 0, 0, false );
	global_lower_mask_matrix.assign    ( 0, 0, false );
	global_selections.clear();
	global_num_selections  =     0;
	global_window          =     0;
	global_window_add      =    70;
	global_res_sim_cutoff  =   150;
	global_run_counter     =     0;
	global_gap_penalty     =    50;
	global_debug           = false;
	global_align_pass      = false;
	global_supaln          =  true;
	global_doing_fast_ssap =  true;
	global_res_score       = false;
	global_frac_selected   =   0.0;
	global_score_run1      =   0.0;
	global_score_run2      =   0.0;
	global_ssap_score1     =   0.0;
	global_ssap_score2     =   0.0;
	fill_n(global_ssap_line1, SSAP_LINE_LENGTH, 0);
	fill_n(global_ssap_line2, SSAP_LINE_LENGTH, 0);
}

/// \brief Temporary setter for global_run_counter to allow tests to check their fixtures are
///        correctly calling reset_ssap_global_variables().
///
/// \todo Eradicate these global variables and then remove this function
void cath::temp_set_global_run_counter(const ptrdiff_t &arg_global_run_counter ///< Value to which global_run_counter should be set
                                       ) {
	global_run_counter = arg_global_run_counter;
}

/// \brief Temporary getter for global_run_counter to allow tests to check their fixtures are
///        correctly calling reset_ssap_global_variables().
///
/// \returns The current value of global_run_counter
///
/// \todo Eradicate these global variables and then remove this function
ptrdiff_t cath::temp_get_global_run_counter() {
	return global_run_counter;
}

/// \brief Read a pair of proteins following the specification in arg_cath_ssap_options
prot_prot_pair cath::read_protein_pair(const cath_ssap_options &arg_cath_ssap_options, ///< The cath_ssap options
                                       ostream                 &arg_stderr             ///< TODOCUMENT
                                       ) {
	const old_ssap_options_block the_ssap_options     = arg_cath_ssap_options.get_old_ssap_options();
	const data_dirs_spec         the_data_dirs        = arg_cath_ssap_options.get_data_dirs_spec();
	const string                 protein_name_a       = the_ssap_options.get_protein_name_a();
	const string                 protein_name_b       = the_ssap_options.get_protein_name_b();
	const auto                   protein_sources_ptr  = the_ssap_options.get_protein_source_files();
	const path_opt               domin_file           = the_ssap_options.get_opt_domin_file();
	return read_protein_pair( protein_name_a, protein_name_b, the_data_dirs, *protein_sources_ptr, domin_file, arg_stderr );
}

/// \brief Read a pair of proteins following the specification in arg_ssap_options
prot_prot_pair cath::read_protein_pair(const string                  &arg_protein_name_a,          ///< TODOCUMENT
                                       const string                  &arg_protein_name_b,          ///< TODOCUMENT
                                       const data_dirs_spec          &arg_data_dirs_spec,          ///< TODOCUMENT
                                       const protein_source_file_set &arg_protein_source_file_set, ///< TODOCUMENT
                                       const path_opt                &arg_domin_file,              ///< TODOCUMENT
                                       ostream                       &arg_stderr                   ///< TODOCUMENT
                                       ) {
	const protein protein_a      = read_protein_data_from_ssap_options_files(arg_data_dirs_spec, arg_protein_name_a, arg_protein_source_file_set, arg_domin_file, arg_stderr );
	const protein protein_b      = read_protein_data_from_ssap_options_files(arg_data_dirs_spec, arg_protein_name_b, arg_protein_source_file_set, none,           arg_stderr );
	return make_pair(protein_a, protein_b);
}

/// \brief SSAP a pair of structures as directed by a cath_ssap_options object
void cath::run_ssap(const cath_ssap_options &arg_cath_ssap_options, ///< The cath_ssap options
                    ostream                 &arg_stdout,            ///< The ostream to which any stdout-like output should be written
                    ostream                 &arg_stderr             ///< The ostream to which any stdout-like output should be written
                    ) {
	// Start be resetting the SSAP global variables
	reset_ssap_global_variables();

	// If the options are invalid or specify to do_nothing, then just return
	const auto error_or_help_string = arg_cath_ssap_options.get_error_or_help_string();
	if ( error_or_help_string ) {
		arg_stderr << error_or_help_string << endl;
		return;
	}

	global_debug = arg_cath_ssap_options.get_old_ssap_options().get_debug();
	const prot_prot_pair proteins = read_protein_pair( arg_cath_ssap_options, arg_stderr );

	global_run_counter = 0;

//	const protein &protein_a = proteins.first;
//	const protein &protein_b = proteins.second;
//	const size_t num_res_a = protein_a.get_length();
//	const size_t num_res_b = protein_b.get_length();
//	for (size_t ctr_a = 0; ctr_a < num_res_a; ++ctr_a) {
//		for (size_t ctr_b = 0; ctr_b < num_res_b; ++ctr_b) {
////			if ( ctr_a + 1 != 20 || ctr_b + 1 != 15 ) {
////				continue;
////			}
//			const residue &residue_a = protein_a.get_residue_ref_of_index( ctr_a );
//			const residue &residue_b = protein_b.get_residue_ref_of_index( ctr_b );
//			const bool result = residues_have_similar_area_angle_props(residue_a, residue_b);
//			cerr << "angle_props :\t" << ctr_a + 1;
//			cerr << "\t" << ctr_b + 1;
//			cerr << "\t" << result;
//			cerr << endl;
//		}
//	}
//	exit(0);


//	const size_t num_res_a = protein_a.get_length();
//	const size_t num_res_b = protein_b.get_length();
//	for (size_t a_view_from_index = 1; a_view_from_index <= num_res_a; ++a_view_from_index) {
//		for (size_t b_view_from_index = 1; b_view_from_index <= num_res_b; ++b_view_from_index) {
//			const residue &residue_a_view_from = get_residue_ref_of_index__offset_1( protein_a, a_view_from_index );
//			const residue &residue_b_view_from = protein_b.get_residue_ref_of_index__offset_1(b_view_from_index );
//			for (size_t a_dest_to_index = 1; a_dest_to_index <= num_res_a; ++a_dest_to_index) {
//				for (size_t b_dest_to_index = 1; b_dest_to_index <= num_res_b; ++b_dest_to_index) {
////					if ( a_view_from_index != 38 || b_view_from_index != 30 || a_dest_to_index != 15 || b_dest_to_index != 2 ) {
////					  continue;
////					}
//					const residue &residue_a_dest_to = protein_a.get_residue_ref_of_index__offset_1( a_dest_to_index );
//					const residue &residue_b_dest_to = protein_b.get_residue_ref_of_index__offset_1( b_dest_to_index );
//					cerr << "res";
//					cerr << "\t" << a_view_from_index;
//					cerr << "\t" << b_view_from_index;
//					cerr << "\t" << a_dest_to_index;
//					cerr << "\t" << b_dest_to_index;
//					cerr << "\t" << context_res(
//						residue_a_view_from, residue_b_view_from,
//						residue_a_dest_to,   residue_b_dest_to
//					);
//					cerr << endl;
//				}
//			}
//		}
//	}

//	const size_t num_ss_a = protein_a.get_num_sec_strucs();
//	const size_t num_ss_b = protein_b.get_num_sec_strucs();
//	for (size_t ctr_a = 0; ctr_a < num_ss_a; ++ctr_a) {
//		for (size_t ctr_b = 0; ctr_b < num_ss_b; ++ctr_b) {
//			const sec_struc &sec_struc_a = protein_a.get_sec_struc_ref_of_index( ctr_a );
//			const sec_struc &sec_struc_b = protein_b.get_sec_struc_ref_of_index( ctr_b );
//			for (size_t ctr_i = 1; ctr_i <= num_ss_a; ++ctr_i) {
//				for (size_t ctr_j = 1; ctr_j <= num_ss_b; ++ctr_j) {
//					cerr << "sec";
//					cerr << "\t" << ctr_a;
//					cerr << "\t" << ctr_b;
//					cerr << "\t" << ctr_i;
//					cerr << "\t" << ctr_j;
//					cerr << "\t" << context_sec(protein_a, protein_b, sec_struc_a, sec_struc_b, ctr_i, ctr_j);
//					cerr << endl;
//				}
//			}
//		}
//	}

	const old_ssap_options_block the_ssap_options = arg_cath_ssap_options.get_old_ssap_options();
	const data_dirs_spec         the_data_dirs    = arg_cath_ssap_options.get_data_dirs_spec();

	// Run SSAP
	align_proteins( proteins.first, proteins.second, the_ssap_options, the_data_dirs );

	// Choose the stream to which to output the results
	//
	// (it's a bit ugly to have this code here but:
	//   - it should eventually be superseded by a bunch of outputter classes, all inheriting from a suitable ABC and
	//   - it's quite difficult to put the functionality in old_ssap_options_block, as would make sense, because
	//     the ostream/ofstream must be returned by ostream reference (or pointer) to avoid slicing, which means the
	//     old_ssap_options_block must own the ofstream but that makes old_ssap_options_block non-copyable, which prevents an option-parsing
	//     function from returning a old_ssap_options_block object (although this would presumably be fine come C++11's move operators).
	ofstream file_out_stream;
	if (the_ssap_options.get_output_to_file()) {
		open_ofstream(file_out_stream, the_ssap_options.get_output_filename());
	}
	ostream &output_stream = the_ssap_options.get_output_to_file() ? file_out_stream : arg_stdout;

	// Print the results
	print_ssap_scores(
		output_stream,
		global_ssap_score1,
		global_ssap_score2,
		global_ssap_line1,
		global_ssap_line2,
		global_run_counter,
		the_ssap_options.get_write_all_scores()
	);
}


/// \brief Align structures
///
/// JEB v1.12 12.09.2002
/// Rewrote this function to separate out running FAST SSAP and SLOW SSAP
/// FAST SSAP performs a comparison of secondary structures first
void cath::align_proteins(const protein                 &arg_protein_a,    ///< The first protein
                          const protein                 &arg_protein_b,    ///< The second protein
                          const old_ssap_options_block  &arg_ssap_options, ///< The old_ssap_options_block to specify how things should be done
                          const data_dirs_spec          &arg_data_dirs     ///< The data directories from which data should be read
                          ) {
	BOOST_LOG_TRIVIAL( info ) << "Function: alnseq";

	// Set alignment options
	global_res_score   = false;
	global_align_pass  = false;
	global_gap_penalty =     5;
	global_window      = max( arg_protein_a.get_num_sec_strucs(), arg_protein_b.get_num_sec_strucs() );

	BOOST_LOG_TRIVIAL( info ) << "Function: alnseq:  seqa->nsec=" << arg_protein_a.get_num_sec_strucs();
	BOOST_LOG_TRIVIAL( info ) << "Function: alnseq:  seqb->nsec=" << arg_protein_b.get_num_sec_strucs();

	ssap_scores fast_ssap_scores;
	if ( !arg_ssap_options.get_slow_ssap_only() ) {
		// Check for minimum number of secondary structures
		if (arg_protein_a.get_num_sec_strucs() > 1 && arg_protein_b.get_num_sec_strucs() > 1) {
			fast_ssap_scores         = fast_ssap(arg_protein_a, arg_protein_b, arg_ssap_options, arg_data_dirs);
			const double first_score = fast_ssap_scores.get_ssap_score_over_larger();

//			if (DEBUG) {
//				cerr << "global2    is " << fast_ssap_scores.get_ssap_score_over_larger()   << endl;
//				cerr << "dist       is " << global_ssap_options.dist        << endl;
//				cerr << "use_clique is " << boolalpha << global_ssap_options.use_clique << endl;
//			}

			// Re-run with increased window and torsional angle cutoffs, if not a great alignment
			if ( first_score < arg_ssap_options.get_max_score_to_fast_ssap_rerun() && ! has_clique_file( arg_ssap_options ) ) {
				BOOST_LOG_TRIVIAL( info ) << "Dist is: " << arg_ssap_options.get_max_score_to_fast_ssap_rerun() << " Removing cutoffs....";

				--global_run_counter;

				// Set alignment options
				// \todo These shouldn't be global variables, they should be parameters to fast_ssap
				global_res_score      = false;
				global_align_pass     = false;
				global_gap_penalty    =     5;
				global_res_sim_cutoff =  1000;
				global_window_add     =  1000;
				global_window         = max( arg_protein_a.get_num_sec_strucs(), arg_protein_b.get_num_sec_strucs() );

				fast_ssap_scores          = fast_ssap(arg_protein_a, arg_protein_b, arg_ssap_options, arg_data_dirs);
				const double second_score = fast_ssap_scores.get_ssap_score_over_larger();

				// Re-run original alignment if it doesn't give a better score
				BOOST_LOG_TRIVIAL( info ) << "Comparing score " << setprecision(30) << second_score << " with " << first_score;
				if (second_score <= first_score) {
					BOOST_LOG_TRIVIAL( info ) << "Reverting back to original Fast SSAP....";

					--global_run_counter;

					// Set alignment options
					// \todo These shouldn't be global variables, they should be parameters to fast_ssap
					global_res_score      = false;
					global_align_pass     = false;
					global_gap_penalty    =     5;
					global_res_sim_cutoff =   150;
					global_window_add     =    70;
					global_window         = max( arg_protein_a.get_num_sec_strucs(), arg_protein_b.get_num_sec_strucs() );

					fast_ssap_scores = fast_ssap(arg_protein_a, arg_protein_b, arg_ssap_options, arg_data_dirs);
				}
			}
		}
	}
	// RUN FAST SSAP - END

	// Check whether the previous fast SSAP result was good
	const bool fast_ssap_result_is_close = ( fast_ssap_scores.get_ssap_score_over_larger() > arg_ssap_options.get_max_score_to_slow_ssap_rerun() );

	// Don't run slow ssap if clique information is used
	const bool run_slow_ssap = ( ! has_clique_file( arg_ssap_options ) && ! fast_ssap_result_is_close );
	if ( fast_ssap_result_is_close ) {
		BOOST_LOG_TRIVIAL( info ) << "Not running slow SSAP. Cutoff: " << arg_ssap_options.get_max_score_to_slow_ssap_rerun()
		                          << " Score: " << fast_ssap_scores.get_ssap_score_over_smaller();
	}

	// RUN SLOW SSAP - START
	if ( run_slow_ssap ) {
		BOOST_LOG_TRIVIAL( info ) << "Function: alnseq:  slow_ssap";

		// v1.14 JEB
		++global_run_counter;

		const size_t max_protein_length = max( arg_protein_a.get_length(), arg_protein_b.get_length() );
		const size_t min_protein_length = min( arg_protein_a.get_length(), arg_protein_b.get_length() );

		// Set variables for SLOW SSAP
		// \todo These shouldn't be global variables, they should be parameters to compare()
		global_res_score       = false;
		global_gap_penalty     =    50;
		global_res_sim_cutoff  =   150;
		global_window_add      =    70;
		global_window          = max_protein_length - min_protein_length + global_window_add;
		global_doing_fast_ssap = false;
		global_num_selections  =     0;

		// Perform two residue alignment passes
		for (size_t pass_ctr = 1; pass_ctr <= 2; ++pass_ctr) {
			BOOST_LOG_TRIVIAL( info ) << "Function: alnseq:  pass=" << pass_ctr;

			global_align_pass = ( pass_ctr > 1 );
			if (pass_ctr == 1 || (pass_ctr == 2 && global_res_score))  {
				compare( arg_protein_a, arg_protein_b, pass_ctr, residue_querier(), arg_ssap_options, arg_data_dirs, none );
			}
		}
	}
	// RUN SLOW SSAP - END

	fflush(stdout);
}				


/// \brief Function to run fast SSAP
ssap_scores cath::fast_ssap(const protein                 &arg_protein_a,    ///< The first protein
                            const protein                 &arg_protein_b,    ///< The second protein
                            const old_ssap_options_block  &arg_ssap_options, ///< The old_ssap_options_block to specify how things should be done
                            const data_dirs_spec          &arg_data_dirs     ///< The data directories from which data should be read
                            ) {
	ssap_scores new_ssap_scores;

	BOOST_LOG_TRIVIAL( info ) << "Fast SSAP: dtot=" << global_res_sim_cutoff << " window_add=" << global_window_add;
	BOOST_LOG_TRIVIAL( info ) << "Function: fast_ssap:  fast_ssap";

	// Perform secondary structure alignment
	++global_run_counter;
	const pair<ssap_scores, alignment> scores_and_alignment = compare( arg_protein_a, arg_protein_b, 1, sec_struc_querier(), arg_ssap_options, arg_data_dirs, none );
	new_ssap_scores = scores_and_alignment.first;
	const alignment &sec_struc_alignment = scores_and_alignment.second;
	fflush(stdout);

	// Check window setting
	const size_t max_protein_length = max( arg_protein_a.get_length(), arg_protein_b.get_length() );
	const size_t min_protein_length = min( arg_protein_a.get_length(), arg_protein_b.get_length() );

	// Align structures using subsets of residue comparisons
	// Set variables for FAST SSAP
	global_align_pass      = false;
	global_gap_penalty     =    50;
	global_window          = max_protein_length - min_protein_length + global_window_add;
	global_doing_fast_ssap =  true;
	global_num_selections  =     0;

	// Perform two residue alignment passes
	for (size_t pass_ctr = 1; pass_ctr <= 2; ++pass_ctr) {
		BOOST_LOG_TRIVIAL( info ) << "Function: fast_ssap:  pass=" << pass_ctr;
		global_align_pass = ( pass_ctr > 1 );
		if ( pass_ctr == 1 || ( pass_ctr == 2 && global_res_score ) ) {
			const pair<ssap_scores, alignment> scores_and_alignment = compare( arg_protein_a, arg_protein_b, pass_ctr, residue_querier(), arg_ssap_options, arg_data_dirs, sec_struc_alignment );
			new_ssap_scores = scores_and_alignment.first;
		}
	}

	return new_ssap_scores;
}


/// \brief Compare structures
pair<ssap_scores, alignment> cath::compare(const protein                 &arg_protein_a,            ///< The first protein
                                           const protein                 &arg_protein_b,            ///< The second protein
                                           const size_t                  &arg_pass_ctr,             ///< The pass of this comparison (where the second typically refines the alignment generated by the first)
                                           const entry_querier           &arg_entry_querier,        ///< The entry_querier to query either residues or secondary structures
                                           const old_ssap_options_block  &arg_ssap_options,         ///< The old_ssap_options_block to specify how things should be done
                                           const data_dirs_spec          &arg_data_dirs,            ///< The data directories from which data should be read
                                           const alignment_opt           &arg_previous_ss_alignment ///< An optional parameter specifying a previous secondary structure alignment
                                           ) {
	const bool   res_not_ss__hacky = arg_entry_querier.temp_hacky_is_residue();
	const string entry_plural_name = get_plural_name(arg_entry_querier);

	const size_t length_a = arg_entry_querier.get_length(arg_protein_a);
	const size_t length_b = arg_entry_querier.get_length(arg_protein_b);

	// Each of these matrices is currently indexed with offset-1
	//
	// \todo Shift each of these matrices to not use offset-1 and remove the extra " + 1"
	//       from these lines
	global_upper_score_matrix.resize   ( length_b + 1, length_a + global_window + 1, 0     );
	global_upper_res_mask_matrix.resize( length_b + 1, length_a + global_window + 1, false );
	global_upper_ss_mask_matrix.resize ( length_b + 1, length_a + global_window + 1, false );
	global_lower_mask_matrix.resize    ( length_b + 1, length_a + global_window + 1, false );

	BOOST_LOG_TRIVIAL( info ) << "Function: compare";
	BOOST_LOG_TRIVIAL( info ) << "Function: compare: [aligning " << entry_plural_name << "]";
	BOOST_LOG_TRIVIAL( info ) << "Function: compare: pass=" << arg_pass_ctr;

	if ( ! res_not_ss__hacky || arg_pass_ctr == 1 ) {
		BOOST_LOG_TRIVIAL( info ) << "Function: compare: [aligning " << entry_plural_name << "] Initialise global_lower_mask_matrix and global_upper_ss_mask_matrix";
		// Each of these matrices is currently indexed with offset-1
		//
		// \todo Shift each of these matrices to not use offset-1 and remove the extra " + 1"
		//       from these lines
		global_upper_ss_mask_matrix.assign( length_b + 1, length_a + global_window + 1, false );
		global_lower_mask_matrix.assign   ( length_b + 1, length_a + global_window + 1, false );
	}

	// Select allowed pairs
	const path_opt clique_file = arg_ssap_options.get_opt_clique_file();
	if ( res_not_ss__hacky && arg_pass_ctr == 1 ) {
		set_mask_matrix(
			arg_protein_a,
			arg_protein_b,
			arg_previous_ss_alignment,
			arg_ssap_options.get_opt_clique_file()
		);
	}

	select_pairs(arg_protein_a, arg_protein_b, arg_pass_ctr, arg_entry_querier);

	// Initialise score matrix to zeros
	//
	// Each of these matrices is currently indexed with offset-1
	//
	// \todo Shift each of these matrices to not use offset-1 and remove the extra " + 1"
	//       from these lines
	BOOST_LOG_TRIVIAL( info ) << "Function: compare: [aligning " << entry_plural_name << "] Initialise global_lower_mask_matrix and global_upper_ss_mask_matrix";
	global_upper_score_matrix.assign   ( length_b + 1, length_a + global_window + 1, 0     );
	// global_upper_res_mask_matrix.assign( length_b + 1, length_a + global_window + 1, false );

	BOOST_LOG_TRIVIAL( info ) << "Function: compare: [aligning " << entry_plural_name << "] score_matrix twice";

	// Call score_matrix() to populate
	populate_upper_score_matrix(arg_protein_a, arg_protein_b, arg_entry_querier, global_align_pass);

	// Construct a source of scores to be used for aligning using dynamic-programming
	// based on the global_upper_score_matrix
	const old_matrix_dyn_prog_score_source upper_score_matrix_score_source(
		global_upper_score_matrix,
		arg_entry_querier.get_length(arg_protein_a),
		arg_entry_querier.get_length(arg_protein_b),
		global_window
	);

	// Align the upper matrix using dynamic-programming
	score_alignment_pair score_and_alignment = ssap_code_dyn_prog_aligner().align(
		upper_score_matrix_score_source,
		gap_penalty( global_gap_penalty, 0 ),
		global_window
	);
	const score_type &score         = score_and_alignment.first;
	alignment        &new_alignment = score_and_alignment.second;

	// Save scores to alignment
	score_opt_vec scores;
	scores.reserve( new_alignment.length() );
	for (size_t alignment_ctr = 0; alignment_ctr < new_alignment.length(); ++alignment_ctr) {
		if ( has_both_positions_of_index( new_alignment, alignment_ctr  )) {
			const aln_posn_type a_position   = get_a_offset_1_position_of_index( new_alignment, alignment_ctr );
			const aln_posn_type b_position   = get_b_offset_1_position_of_index( new_alignment, alignment_ctr );
			const int           a_matrix_idx = get_window_matrix_a_index__offset_1(length_a, length_b, global_window, a_position, b_position);
			const double        score        = numeric_cast<double>( global_upper_score_matrix.get( b_position, numeric_cast<size_t>( a_matrix_idx ) ) );
			scores.push_back( score / 10.0 + 0.5 );
//			cerr << "Retrieved score:\t" << score << ",\twhich normalises to: " << ( score / 10.0 + 0.5 ) << endl;
		}
		else {
			scores.push_back( none );
		}
	}
	set_pair_alignment_duplicate_scores( new_alignment, scores );

	ssap_scores new_ssap_scores;
	if (score) {
		new_ssap_scores = plot_aln(
			arg_protein_a,
			arg_protein_b,
			arg_pass_ctr,
			arg_entry_querier,
			new_alignment,
			arg_ssap_options,
			arg_data_dirs
		);
	}

	if (res_not_ss__hacky) {
		if (score) {
			global_res_score = true;
		}
		else {
			// BOOST_LOG_TRIVIAL( warning ) << "Saving zero scores after an attempted alignment."
			//                                 " This likely indicates a problem. If you think it does"
			//                                 " and none of the previous messages has referenced an existing GitHub Issue"
			//                                 " (or if you think there isn't any problem and this message is spurious)"
			//                                 " please consider raising a new issue at https://github.com/UCLOrengoGroup/cath-tools/issues";

			// v1.14 JEB - Save zero scores
			save_zero_scores(arg_protein_a, arg_protein_b);
			global_res_score = false;
		}
	}

	return make_pair(new_ssap_scores, new_alignment);
}

/// \brief Read data for a protein based on its name and a old_ssap_options_block object
protein cath::read_protein_data_from_ssap_options_files(const data_dirs_spec          &arg_data_dirs,               ///< The old_ssap_options_block to specify how things should be done
                                                        const string                  &arg_protein_name,            ///< The name of the protein that is to be read from files
                                                        const protein_source_file_set &arg_protein_source_file_set, ///< TODOCUMENT
                                                        const path_opt                &arg_domin_file,              ///< Optional domin file
                                                        ostream                       &arg_stderr                   ///< TODOCUMENT
                                                        ) {
	// Report which files are being used
	const data_file_path_map filename_of_data_file = get_filename_of_data_file(
		arg_protein_source_file_set,
		arg_data_dirs,
		arg_protein_name
	);
	for (const data_file_path_pair &filename_and_data_file : filename_of_data_file) {
		const string file_str              = to_lower_copy( lexical_cast<string>( filename_and_data_file.first ) );
		const string right_padded_file_str = string( max_data_file_str_length() - file_str.length(), ' ' );
		BOOST_LOG_TRIVIAL( info ) << "Loading " << file_str << right_padded_file_str << " from " << filename_and_data_file.second;
	}

	// Create a protein object from the name, wolf file and sec file
	protein new_protein_to_populate = arg_protein_source_file_set.read_files(
		arg_data_dirs,
		arg_protein_name,
		arg_stderr
	);

	// Re-calculate if there is a domin file
	if ( arg_domin_file ) {
		remove_domin_res( new_protein_to_populate, *arg_domin_file, arg_stderr);
	}

	//	// Return the newly created protein object
	return new_protein_to_populate;
}

/// \brief Read Clique file
clique cath::read_clique_file(const path &arg_filename ///< The clique file to read
                              ) {

	const size_t MAX_CLIQUE_FILE_BUFFER_LENGTH(10000);
	char_vec buffer(MAX_CLIQUE_FILE_BUFFER_LENGTH, 0);

	FILE *in;
	if ( ( in = fopen( arg_filename.string().c_str(),"r" ) ) == nullptr ) {
		BOOST_LOG_TRIVIAL( error ) << "Clique file (" << arg_filename << ") not found!";
		BOOST_THROW_EXCEPTION(invalid_argument_exception("**** TEMPORARY ***** TEMPORARY ***** TEMPORARY ***** TEMPORARY *****"));
		exit(1);
	}

	clique new_clique_file;
	
	// Read clique size
	const char * const buffer_fgets_return_1 = fgets(&buffer.front(), MAX_CLIQUE_FILE_BUFFER_LENGTH-1, in);
	if (buffer_fgets_return_1 == nullptr) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Parsing error in reading clique file"));
	}
	sscanf(&buffer.front(), "%zu", &new_clique_file.cliquesize);
	for (size_t clique_ctr = 0; clique_ctr < new_clique_file.cliquesize; ++clique_ctr) {
		const char * const buffer_fgets_return_2 = fgets(&buffer.front(), MAX_CLIQUE_FILE_BUFFER_LENGTH-1, in);
		if (buffer_fgets_return_2 != nullptr) {
			BOOST_THROW_EXCEPTION(runtime_error_exception("Parsing error in reading clique file"));
		}
		sec_struc_equivalency new_equivalency;
		sscanf(
			&buffer.front(),
			"%zu %s %s %zu %s %s\n",
			&new_equivalency.prota_ssnum,
			 new_equivalency.prota_start,
			 new_equivalency.prota_end,
			&new_equivalency.protb_ssnum,
			 new_equivalency.protb_start,
			 new_equivalency.protb_end
		);
		new_clique_file.equivs.push_back(new_equivalency);
	}
	fclose(in);

	return new_clique_file;
}


/// \brief Prepare the mask matrices for the next comparison
///
/// This currently only gets called from one location, which is when performing the first
/// pass of a residue comparison
void cath::set_mask_matrix(const protein       &arg_protein_a,        ///< The first protein
                           const protein       &arg_protein_b,        ///< The second protein
                           const alignment_opt &arg_opt_ss_alignment, ///< A secondary structure alignment that is required in some modes so that it can be transferred to a residue mask matrix
                           const path_opt      &arg_clique_file       ///< An optional clique file to use
                           ) {
	const size_t length_a = arg_protein_a.get_length();
	const size_t length_b = arg_protein_b.get_length();

	// Initialise arrays
	for (size_t residue_ctr_b = 0; residue_ctr_b <= length_b; ++residue_ctr_b ) {
		for (size_t residue_ctr_a = 0; residue_ctr_a <= length_a; ++residue_ctr_a ) {
			global_upper_res_mask_matrix.set( residue_ctr_b, residue_ctr_a, false );
			global_upper_ss_mask_matrix.set ( residue_ctr_b, residue_ctr_a, false );
			global_lower_mask_matrix.set    ( residue_ctr_b, residue_ctr_a, false );
		}
	}

	// If using clique file
	if ( arg_clique_file ) {
		// Read clique file
		const clique clique_data = read_clique_file( *arg_clique_file );
		const size_t clique_size = clique_data.equivs.size();

		// Ollie's variables
		const int boundary = 5; // Amount to add to boundary

		// Set equivalent residues in secondary structures for protein A and protein b
		for (size_t ctr_b = length_b; ctr_b > 0; --ctr_b) {
			const residue &residue_b = get_residue_ref_of_index__offset_1(arg_protein_b, ctr_b);
			for (size_t ctr_a = length_a; ctr_a > 0; --ctr_a) {
				const residue residue_a = get_residue_ref_of_index__offset_1(arg_protein_a, ctr_a);

				// Look to see if they match any secondary structures
				for (size_t k = 0; k < clique_size; ++k) {
					int bstart = atoi( clique_data.equivs[k].protb_start ) - boundary;
					int bend   = atoi( clique_data.equivs[k].protb_end   ) + boundary;
					int astart = atoi( clique_data.equivs[k].prota_start ) - boundary;
					int aend   = atoi( clique_data.equivs[k].prota_end   ) + boundary;

					if (get_pdb_name_number( residue_a )           && get_pdb_name_number( residue_b )         &&
						get_pdb_name_number( residue_b ) >= bstart && get_pdb_name_number( residue_b ) <= bend &&
						get_pdb_name_number( residue_a ) >= astart && get_pdb_name_number( residue_a ) <= aend) {
						global_lower_mask_matrix.set( ctr_b, ctr_a, true );
						break;
					}
				}

				// Match regions between secondary structures
				for (size_t k = 0; k < clique_size - 1; ++k) {
					int bstart = atoi( clique_data.equivs[ k + 1 ].protb_start ) + boundary;
					int bend   = atoi( clique_data.equivs[ k     ].protb_end   ) - boundary;
					int astart = atoi( clique_data.equivs[ k + 1 ].prota_start ) + boundary;
					int aend   = atoi( clique_data.equivs[ k     ].prota_end   ) - boundary;

					if (get_pdb_name_number( residue_a )          && get_pdb_name_number( residue_b )        &&
						get_pdb_name_number( residue_b ) < bstart && get_pdb_name_number( residue_b ) > bend &&
						get_pdb_name_number( residue_a ) < astart && get_pdb_name_number( residue_a ) > aend) {
						global_lower_mask_matrix.set( ctr_b, ctr_a, true );
						break;
					}
				}
			}
		}

		// Find first and last residues in the clique
		const int firsta  = atoi( clique_data.equivs[ 0               ].prota_start ) + boundary;
		const int firstb  = atoi( clique_data.equivs[ 0               ].protb_start ) + boundary;
		const int lasta   = atoi( clique_data.equivs[ clique_size - 1 ].prota_end   ) - boundary;
		const int lastb   = atoi( clique_data.equivs[ clique_size - 1 ].protb_end   ) - boundary;

		// Set equivalent pairs at beginning and end of alignment
		for (size_t ctr_b = length_b; ctr_b > 0; --ctr_b) {
			const residue &residue_b = get_residue_ref_of_index__offset_1(arg_protein_b, ctr_b);
			for (size_t ctr_a = length_a; ctr_a > 0; --ctr_a) {
				const residue &residue_a = get_residue_ref_of_index__offset_1(arg_protein_a, ctr_a);

				// Tail end of alignment
				if ( get_pdb_name_number( residue_a ) > lasta  && get_pdb_name_number( residue_b ) > lastb ) {
					global_lower_mask_matrix.set( ctr_b, ctr_a, true );
				}
				// Start of alignment
				if ( get_pdb_name_number( residue_a ) < firsta && get_pdb_name_number( residue_b ) < firstb ) {
					global_lower_mask_matrix.set( ctr_b, ctr_a, true );
				}
			}
		}
	}

	// Save aligned secondary structure pairs in sec_struc_match_matrix
	//
	// Note that vector<bool> is not like other vectors because the std infamously made
	// the mistake of specialising vector<bool> to be a bitset-style class that can't
	// return actually a reference to a bool. However none of that matters here, so
	// it's OK to use vector<bool>
	bool_vec_of_vec sec_struc_match_matrix;
	if ( arg_opt_ss_alignment ) {
		const auto last_present_a_opt = get_last_present_a_position( *arg_opt_ss_alignment );
		const auto last_present_b_opt = get_last_present_b_position( *arg_opt_ss_alignment );
		if ( last_present_a_opt && last_present_b_opt ) {
			sec_struc_match_matrix.assign(
				*last_present_b_opt + 2, // ( + 1 to go from position to size and +1 because it will be indexed offset_1)
				*last_present_a_opt + 2, // ( + 1 to go from position to size and +1 because it will be indexed offset_1)
				false
			);
			BOOST_LOG_TRIVIAL( debug ) << "Setting secondary structure alignment : " << *arg_opt_ss_alignment;;
			for (const size_t &alignment_ctr : irange( 0_z, arg_opt_ss_alignment->length() ) ) {
				if ( has_both_positions_of_index( *arg_opt_ss_alignment, alignment_ctr ) ) {
					sec_struc_match_matrix.set(
						get_b_offset_1_position_of_index( *arg_opt_ss_alignment, alignment_ctr ),
						get_a_offset_1_position_of_index( *arg_opt_ss_alignment, alignment_ctr ),
						true
					);
				}
			}
		}
	}

	global_num_selections = 0;
	size_t total_num_residues_considered = 0;
	size_t num_residues_selected         = 0;
	for (size_t residue_ctr_b = length_b; residue_ctr_b > 0; --residue_ctr_b) {
		const residue &residue_b    = get_residue_ref_of_index__offset_1(arg_protein_b, residue_ctr_b);
		const size_t   window_start = get_window_start_a_for_b__offset_1( length_a, length_b, global_window, residue_ctr_b );
		const size_t   window_stop  = get_window_stop_a_for_b__offset_1 ( length_a, length_b, global_window, residue_ctr_b );

		for (size_t residue_ctr_a = window_stop; residue_ctr_a >= window_start; --residue_ctr_a ) {
			const residue &residue_a    = get_residue_ref_of_index__offset_1(arg_protein_a, residue_ctr_a);
			const int      a_matrix_idx = get_window_matrix_a_index__offset_1(length_a, length_b, global_window, residue_ctr_a, residue_ctr_b);

			++total_num_residues_considered;

			// IF USING SEC STR. ALIGNMENT TO GUIDE RESIDUE SELECTION
			if ( global_doing_fast_ssap ) {
				// Use clique method
				if ( arg_clique_file ) {
					if ( global_lower_mask_matrix.get( residue_ctr_b, residue_ctr_a ) && residues_have_similar_area_angle_props( residue_a, residue_b ) ) {
						++num_residues_selected;
						global_upper_res_mask_matrix.set( residue_ctr_b, numeric_cast<size_t>( a_matrix_idx ), true );
					}
				}
				// If no clique data is present, use built-in secondary structure method
				else if (residue_a.get_sec_struc_number()
				         && residue_b.get_sec_struc_number()
				         && arg_opt_ss_alignment
				         && sec_struc_match_matrix.get( residue_b.get_sec_struc_number(), residue_a.get_sec_struc_number() )
				         && residues_have_similar_area_angle_props(residue_a, residue_b) ) {
					++num_residues_selected;
					global_upper_res_mask_matrix.set( residue_ctr_b, numeric_cast<size_t>( a_matrix_idx ), true );
				}
			}
			else {
				if (residues_have_similar_area_angle_props(residue_a, residue_b)) {
					++num_residues_selected;
					global_upper_res_mask_matrix.set( residue_ctr_b, numeric_cast<size_t>( a_matrix_idx ), true );
				}
			}
		}
	}
	global_frac_selected = numeric_cast<double>( num_residues_selected ) / numeric_cast<double>( total_num_residues_considered );
}


/// \brief Selects residue pairs in similar structural locations or secondary structures of same type
///
/// This sets global_lower_mask_matrix and possibly global_upper_ss_mask_matrix with the selections
void cath::select_pairs(const protein       &arg_protein_a,    ///< The first protein
                        const protein       &arg_protein_b,    ///< The second protein
                        const size_t        &arg_pass,         ///< The pass of this comparison (where the second typically refines the alignment generated by the first)
                        const entry_querier &arg_entry_querier ///< The entry_querier to query either residues or secondary structures
                        ) {
	const size_t length_a = arg_entry_querier.get_length(arg_protein_a);
	const size_t length_b = arg_entry_querier.get_length(arg_protein_b);

	deque<selected_pair> selected_pairs;

	// Reset variables/arrays for selected residue pairs
	size_t num_entries_selected         = 0;
	size_t total_num_entries_considered = 0;

	// Compare properties of residue/SS pairs for each cell in matrix window
	for (size_t ctr_b = length_b; ctr_b > 0; --ctr_b) {
		const size_t window_start = get_window_start_a_for_b__offset_1( length_a, length_b, global_window, ctr_b );
		const size_t window_stop  = get_window_stop_a_for_b__offset_1(  length_a, length_b, global_window, ctr_b );

		for (size_t ctr_a = window_stop; ctr_a >= window_start; --ctr_a ) {
			++total_num_entries_considered;
			const int a_matrix_idx = get_window_matrix_a_index__offset_1(length_a, length_b, global_window, ctr_a, ctr_b);
	
			// First pass:
			//   for residues:             select if areas/angles similar
			//   for secondary structures: select if both are of same type
			if ( arg_pass == 1 ) {
				if ( arg_entry_querier.are_similar__offset_1(arg_protein_a, arg_protein_b, ctr_a, ctr_b) ) {
					++num_entries_selected;
					global_lower_mask_matrix.set   ( ctr_b, ctr_a, true );
					global_upper_ss_mask_matrix.set( ctr_b, ctr_a, true );
				}
				else {
					global_lower_mask_matrix.set   ( ctr_b, ctr_a, false );
					global_upper_ss_mask_matrix.set( ctr_b, ctr_a, false );
				}
			}
			// Subsequent passes (must be residues):
			//   select 20 highest scoring residue pairs from first pass
			else {
				const score_type score = global_upper_score_matrix.get( ctr_b, numeric_cast<size_t>( a_matrix_idx ) );
				update_best_pair_selections( selected_pairs, selected_pair(ctr_a, ctr_b, score), NUM_SELECTIONS_TO_SAVE );
			}
		}
	}

	// For second pass and residue comparisons, copy selected residues into select structure
	if ( global_align_pass && arg_entry_querier.temp_hacky_is_residue() ) {
		global_selections.assign( NUM_SELECTIONS_TO_SAVE + 1, make_pair( 0_z, 0_z ) );
		for (size_t selected_ctr = 0; selected_ctr < selected_pairs.size(); ++selected_ctr) {
			// Index is calculated to put the selection at the end of the positions with indices 1..NUM_TO_SAVE
			const size_t index_in_global_selections = NUM_SELECTIONS_TO_SAVE + 1 - ( selected_pairs.size() - selected_ctr );
			global_selections[ index_in_global_selections ] = make_pair(
				selected_pairs[ selected_ctr ].get_index_a(),
				selected_pairs[ selected_ctr ].get_index_b()
			);
		}
		global_num_selections = NUM_SELECTIONS_TO_SAVE;
	}

	// Calculate fraction of total residue pairs selected
	if ( arg_pass > 1 ) {
		num_entries_selected = NUM_SELECTIONS_TO_SAVE;
	}
	if ( global_align_pass && arg_entry_querier.temp_hacky_is_residue()) {
		global_frac_selected = numeric_cast<double>( num_entries_selected ) / numeric_cast<double>( total_num_entries_considered );
	}
}


/// \brief Potentially update a limited list of best seen pairs with a new entry
///        (ie replace the worst if the list's already full or just add otherwise)
///
/// \todo Move the lines that set global_lower_mask_matrix out of this subroutine
void cath::update_best_pair_selections(deque<selected_pair> &arg_selected_pairs,    ///< The best scoring pairs so far, in ascending order by score
                                       const selected_pair  &arg_potential_pair,    ///< A potential new pair
                                       const size_t         &arg_max_num_selections ///< The maximum number of selections to store
                                       ) {
	const size_t index_a = arg_potential_pair.get_index_a();
	const size_t index_b = arg_potential_pair.get_index_b();
	global_lower_mask_matrix.set( index_b, index_a, false );

	// If arg_selected_pairs isn't yet full or if the new score is better than the lowest score
	// (which comes first because arg_selected_pairs is sorted) then...
	const bool full            =  arg_selected_pairs.size()  >= arg_max_num_selections;
	const bool new_beats_first = !arg_selected_pairs.empty() && (arg_potential_pair.get_score() > arg_selected_pairs.front().get_score());
	if (!full || new_beats_first) {
		// Add the new entry to the front and then re-sort
		arg_selected_pairs.push_front( arg_potential_pair );
		stable_sort( arg_selected_pairs );

		// If arg_max_num_selections was already full, remove the first entry
		if (full) {
			const size_t first_index_a = arg_selected_pairs.front().get_index_a();
			const size_t first_index_b = arg_selected_pairs.front().get_index_b();
			global_lower_mask_matrix.set( first_index_b, first_index_a, false );

			arg_selected_pairs.pop_front();
		}

		global_lower_mask_matrix.set( index_b, index_a, true );
	}
}


/// \brief Check whether residue pair have similar area/angle properties.
///
/// Current globals used:
///  - global_residue_similarity_cutoff
///
/// \todo Consider potential problems in this code:
///       -# the code checks the sum of accessibilities rather than the difference which makes little sense
///          (although the difference is implied in buried_difference)
///       -# Natalie has identified a pair of very similar structures for which cath-ssap gives all zeroes
///          that are handled correctly if the difference is used here rather than the sum.
///          See: https://github.com/UCLOrengoGroup/cath-tools/issues/8
///       -# the code doesn't allow for wrapping of phi and psi angles
///       -# the code doesn't do anything to handle undetermined phi/psi angles at breaks in the chain
///          (which, at present, get set to 360.0)
bool cath::residues_have_similar_area_angle_props(const residue &arg_residue_i, ///< The first  residue to compare
                                                  const residue &arg_residue_j  ///< The second residue to compare
                                                  ) {
	const int    buried_i                   = get_accessi_of_residue( arg_residue_i );
	const int    buried_j                   = get_accessi_of_residue( arg_residue_j );
	const size_t buried_difference          = numeric_cast<size_t>( difference( buried_i, buried_j ) );
	const size_t phi_angle_diff_in_degrees  = numeric_cast<size_t>( round( difference(
		angle_in_degrees( arg_residue_i.get_phi_angle() ),
		angle_in_degrees( arg_residue_j.get_phi_angle() )
	) ) );
	const size_t psi_angle_diff_in_degrees  = numeric_cast<size_t>( round( difference(
		angle_in_degrees( arg_residue_i.get_psi_angle() ),
		angle_in_degrees( arg_residue_j.get_psi_angle() )
	) ) );

	const size_t mean_angle_diff_in_degrees = ( phi_angle_diff_in_degrees + psi_angle_diff_in_degrees ) / 2;
	const size_t accessibility_sum          = arg_residue_i.get_access() + arg_residue_j.get_access();
//	const size_t accessibility_difference   = difference( arg_residue_i.get_access(), arg_residue_j.get_access() );

//	cerr << "Phi     a                  : " << arg_residue_i.get_phi_angle()               << endl;
//	cerr << "Phi     b                  : " << arg_residue_j.get_phi_angle()               << endl;
//	cerr << "Psi     a                  : " << arg_residue_i.get_psi_angle()               << endl;
//	cerr << "Psi     b                  : " << arg_residue_j.get_psi_angle()               << endl;
//	cerr << "Average phi/psi difference : " << mean_angle_diff_in_degrees                  << endl;
//	cerr << "Access  a                  : " << arg_residue_i.get_access()                  << endl;
//	cerr << "Access  b                  : " << arg_residue_j.get_access()                  << endl;
//	cerr << "Residue a                  : " << arg_residue_i.get_amino_acid().get_letter() << endl;
//	cerr << "Residue b                  : " << arg_residue_j.get_amino_acid().get_letter() << endl;
//	cerr << "Buried  a                  : " << buried_i                                    << endl;
//	cerr << "Buried  b                  : " << buried_j                                    << endl;
//	cerr << "Buried difference          : " << buried_difference                           << endl;

	// Combined areas and angles
	return ( buried_difference + accessibility_sum        + mean_angle_diff_in_degrees < global_res_sim_cutoff );
//	return ( buried_difference + accessibility_difference + mean_angle_diff_in_degrees < global_res_sim_cutoff );
}

/// \brief Populate the scores for the upper (ie major, whole) matrix
///
/// This iterates over certain cells in the upper matrix and calls compare_upper_cell()
/// on them, which does Dynamic Programming (DP) on the views from that pair and then
/// adds the individual scores along that alignment to the upper matrix.
///
/// \pre Presumably global_upper_score_matrix must be zeroed
///
/// \post global_upper_score_matrix will have appropriate scores added to it
///
/// This code used to be incorporated into score_matrix and has been separated out,
/// making both quite a bit easier to understand.
///
/// For an align_pass of residues, only the top-scoring selections are considered.
///
/// For other cases, a mask (global_upper_res_mask_matrix, global_upper_ss_mask_matrix
/// or global_lower_mask_matrix) is used to determine which cells are considered.
///
/// \todo In general, abstract matrix iteration into a class so that:
///         - different matrix-iterating pieces of code don't need to repeat
///           calculations and double loops
///         - the selection-based iteration used in this function can be injected
///           via an iterator argument (ie dependency injection) so that this code
///           doesn't need nasty conditional code to handle it
///       Ensure that the standard matrix iteration matches the sweep that's required
///       by the dynamic-programming code in score_matrix().
///
/// \todo For this function, ensure that the particular masking behaviour is also dependency-injected
void cath::populate_upper_score_matrix(const protein       &arg_protein_a,     ///< The first protein
                                       const protein       &arg_protein_b,     ///< The second protein
                                       const entry_querier &arg_entry_querier, ///< The entry_querier to query either residues or secondary structures
                                       const bool          &arg_align_pass     ///< Whether this is a later, alignment-refining pass
                                       ) {
	// If this is a later, alignment-refining pass of residue this use the selected
	// set of top-scoring residue pairs
	const bool res_not_ss__hacky = arg_entry_querier.temp_hacky_is_residue();
	const bool using_selections  = (res_not_ss__hacky && arg_align_pass);

	// Set number of elements in protein A and B to compare

	const size_t full_length_a = arg_entry_querier.get_length(arg_protein_a);
	const size_t full_length_b = arg_entry_querier.get_length(arg_protein_b);
	const size_t length_a      =                                            full_length_a;
	const size_t length_b      = using_selections ? global_num_selections : full_length_b;

	// Set normalisation constant
	//
	// Note: it's not very clear where these two constants come from and the value
	//       only seems to get used in compare_upper_cell() if comparing residues
	//       (not secondary structures) anyway.
	const double normalisation_num = res_not_ss__hacky ? 200.0 : 25.0;
	const double normalisation     = global_frac_selected * sqrt( normalisation_num * numeric_cast<double>( min( length_a, length_b ) ) );

	size_t num_potential_upper_cell_comps = 0;
	size_t num_actual_upper_cell_comps    = 0;
	bool   found_non_zero_cell            = false;
	bool   found_threshold_cell           = false;

	// Reverse-iterate over the elements in arg_protein_b
	// (or over the selections if using them)
	for (size_t ctr_b = length_b; ctr_b > 0; --ctr_b) {
		// Calculate the arg_protein_a window start/stop for this arg_protein_b entry
		// (or just set them both from the selected pair if using selections)
		const size_t window_start = using_selections ? global_selections[ctr_b].first
		                                             : get_window_start_a_for_b__offset_1( length_a, length_b, global_window, ctr_b );
		const size_t window_stop  = using_selections ? global_selections[ctr_b].first
		                                             : get_window_stop_a_for_b__offset_1(  length_a, length_b, global_window, ctr_b );
		const size_t jval         = using_selections ? global_selections[ctr_b].second
		                                             : ctr_b;

		// Iterate over the window that's been calculated
		for (size_t ctr_a = window_stop; ctr_a >= window_start; --ctr_a ) {
			// Determine whether this pair should be compared:
			//  - If using selections,           then true, else
			//  - If using residues,             then consult global_upper_res_mask_matrix, else
			//  -    Using secondary structures, so   consult global_upper_ss_mask_matrix
			bool should_compare_pair = true;
			if ( ! using_selections ) {
				if ( res_not_ss__hacky ) {
					const int a_matrix_idx = get_window_matrix_a_index__offset_1( length_a, length_b, global_window, ctr_a, ctr_b );
					should_compare_pair = global_upper_res_mask_matrix.get( ctr_b, numeric_cast<size_t>( a_matrix_idx ) );
				}
				else {
					should_compare_pair = global_upper_ss_mask_matrix.get( ctr_b, ctr_a );
				}
			}

			// Compare environments of allowed pairs
			++num_potential_upper_cell_comps;
			if ( should_compare_pair ) {
				++num_actual_upper_cell_comps;
				const auto compare_result = compare_upper_cell(
					arg_protein_a,
					arg_protein_b,
					ctr_a,
					jval,
					arg_entry_querier,
					normalisation
				);
				found_non_zero_cell  = found_non_zero_cell  || ( compare_result != compare_upper_cell_result::ZERO   );
				found_threshold_cell = found_threshold_cell || ( compare_result == compare_upper_cell_result::SCORED );
			}
		}
	}


	const string msg_context_prfx = "When populating upper_score_matrix ("
	                                + arg_entry_querier.get_entry_name()
	                                + "; pass "
	                                + ::std::to_string( arg_align_pass )
	                                + "), ";
	BOOST_LOG_TRIVIAL( debug ) << msg_context_prfx
	                           << "compared "
	                           << num_actual_upper_cell_comps
	                           << " residue pairs out of a possible "
	                           << num_potential_upper_cell_comps;
	if ( res_not_ss__hacky && ! arg_align_pass ) {
		if ( num_actual_upper_cell_comps == 0 ) {
			BOOST_LOG_TRIVIAL( warning ) << msg_context_prfx
			                             << "chose no residue pairs out of a possible "
			                             << num_potential_upper_cell_comps
			                             << " to compare."
			                                " This may relate to https://github.com/UCLOrengoGroup/cath-tools/issues/8"
			                                " - please see that issue for more information and please add a comment"
			                                " if it's causing you problems (or open a new issue if this message is spurious).";
		}
		else if ( ! found_threshold_cell ) {
			if ( found_non_zero_cell ) {
				BOOST_LOG_TRIVIAL( warning ) << msg_context_prfx
				                             << "attempted alignment for "
				                             << num_potential_upper_cell_comps
				                             << " cells in the upper matrix and though some achieved non-zero scores,"
				                                " none of them reached the threshold after their normalisation";
			}
			else {
				BOOST_LOG_TRIVIAL( warning ) << msg_context_prfx
				                             << "attempted alignment for "
				                             << num_potential_upper_cell_comps
				                             << " cells in the upper matrix but none of them achieved non-zero scores";
			}
		}
	}
}

/// \brief Compares residue environments in lower level matrix, if score above threshold,
///        adds alignment path to upper level matrix
///
/// \todo Figure out what's going on
compare_upper_cell_result cath::compare_upper_cell(const protein       &arg_protein_a,                   ///< The first  protein
                                                   const protein       &arg_protein_b,                   ///< The second protein
                                                   const size_t        &arg_a_view_from_index__offset_1, ///< The index of the residue/secondary-structure in the first  protein on which this should be performed
                                                   const size_t        &arg_b_view_from_index__offset_1, ///< The index of the residue/secondary-structure in the second protein on which this should be performed
                                                   const entry_querier &arg_entry_querier,               ///< The entry_querier to query either residues or secondary structures
                                                   const double        &arg_normalisation                ///< The value that should be used to normalise the score for residues before comparison against MIN_LOWER_MAT_RES_SCORE
                                                   ) {
	const bool   res_not_ss__hacky = arg_entry_querier.temp_hacky_is_residue();
	const size_t length_a          = arg_entry_querier.get_length(arg_protein_a);
	const size_t length_b          = arg_entry_querier.get_length(arg_protein_b);

	// Construct two sources of scores to be used for aligning using dynamic-programming:
	//  * the first just uses arg_entry_querier, arg_a_view_from_index and arg_b_view_from_index
	//  * the second is a masked version of the first, using global_lower_mask_matrix
	check_offset_1(arg_a_view_from_index__offset_1);
	check_offset_1(arg_b_view_from_index__offset_1);
	const entry_querier_dyn_prog_score_source entry_querier_score_source(
		arg_entry_querier,
		arg_protein_a,
		arg_protein_b,
		arg_a_view_from_index__offset_1 - 1,
		arg_b_view_from_index__offset_1 - 1
	);
	const mask_dyn_prog_score_source mask_score_source(
		global_lower_mask_matrix,
		entry_querier_score_source
	);

	// Choose between the two score sources:
	//  * if this is an aligning pass, then use entry_querier_score_source;
	//  * otherwise, use mask_score_source, which is like entry_querier_score_source but masked
	const dyn_prog_score_source &the_score_source = global_align_pass ? static_cast<const dyn_prog_score_source &>(entry_querier_score_source)
	                                                                  : static_cast<const dyn_prog_score_source &>(mask_score_source);

	// Align the lower matrix using dynamic-programming
	score_alignment_pair score_and_alignment = ssap_code_dyn_prog_aligner().align(
		the_score_source,
		gap_penalty(global_gap_penalty, 0),
		global_window
	);
	score_type       score        = score_and_alignment.first;
	const alignment &my_alignment = score_and_alignment.second;

//	cerr << "Comparing\t" << arg_a_view_from_index << "\t" << arg_b_view_from_index << "\t" << arg_entry_querier.get_entry_name();
//	cerr << "\tscore is " << score << "\twith alignment length " << my_alignment.length() << endl;

	// Check whether normalised score is above threshold
	if ( res_not_ss__hacky ) {
		if ( arg_normalisation != 0.0 ) {
			score = numeric_cast<score_type>( numeric_cast<double>( score ) / arg_normalisation );
		}
		else {
			score = 0;
		}
	}

	if ( score == 0 ) {
		return compare_upper_cell_result::ZERO;
	}
	if ( res_not_ss__hacky && score < MIN_LOWER_MAT_RES_SCORE ) {
		return compare_upper_cell_result::NON_ZERO_BELOW_THRESHOLD;
	}

	// If yes, trace distance (lower) level alignment path, onto residue (upper) level matrix
	for (size_t alignment_ctr = 0; alignment_ctr < my_alignment.length(); ++alignment_ctr) {
		if (has_both_positions_of_index(my_alignment, alignment_ctr)) {
			const aln_posn_type a_dest_to_index__offset_1 = get_a_offset_1_position_of_index( my_alignment, alignment_ctr );
			const aln_posn_type b_dest_to_index__offset_1 = get_b_offset_1_position_of_index( my_alignment, alignment_ctr );
			const int           a_matrix_idx              = get_window_matrix_a_index__offset_1(length_a, length_b, global_window, a_dest_to_index__offset_1, b_dest_to_index__offset_1);
			const score_type    score_addend              = arg_entry_querier.distance_score__offset_1(
				arg_protein_a,                   arg_protein_b,
				arg_a_view_from_index__offset_1, arg_b_view_from_index__offset_1,
				a_dest_to_index__offset_1,       b_dest_to_index__offset_1
			);
			global_upper_score_matrix.get( b_dest_to_index__offset_1, numeric_cast<size_t>( a_matrix_idx ) ) += score_addend;
//			cerr << "At\t" << ( arg_a_view_from_index__offset_1 - 1 );
//			cerr << "\t"   << ( arg_b_view_from_index__offset_1 - 1 );
//			cerr << "\t"   << ( a_dest_to_index__offset_1       - 1 );
//			cerr << "\t"   << ( b_dest_to_index__offset_1       - 1 );
//			cerr << "\tadding score:\t" << score_addend;
//			cerr << "\tto get:\t" << global_upper_score_matrix[b_dest_to_index__offset_1][ numeric_cast<size_t>( a_matrix_idx ) ];
//			cerr <<"\t["   << get_plural_name(arg_entry_querier) << "]" << endl;
		}
	}
	return compare_upper_cell_result::SCORED;
}



/// \brief Compares vectors/scalars/overlap/packing between secondary structures in two proteins
score_type cath::context_sec(const protein &arg_protein_a,         ///< The first  protein
                             const protein &arg_protein_b,         ///< The second protein
                             const size_t  &arg_a_view_from_index, ///< The "from" secondary structure in the first  protein
                             const size_t  &arg_b_view_from_index, ///< The "from" secondary structure in the second protein
                             const size_t  &arg_to_ss_index_a,     ///< The index of the "to" secondary structure in the first  protein
                             const size_t  &arg_to_ss_index_b      ///< The index of the "to" secondary structure in the second protein
                             ) {
	const sec_struc               &from_sec_struc_a = arg_protein_a.get_sec_struc_ref_of_index( arg_a_view_from_index );
	const sec_struc               &from_sec_struc_b = arg_protein_b.get_sec_struc_ref_of_index( arg_b_view_from_index );
	const sec_struc               &to_sec_struc_a   = arg_protein_a.get_sec_struc_ref_of_index( arg_to_ss_index_a     );
	const sec_struc               &to_sec_struc_b   = arg_protein_b.get_sec_struc_ref_of_index( arg_to_ss_index_b     );
	const sec_struc_planar_angles &planar_angles_a  = from_sec_struc_a.get_planar_angles_of_index( arg_to_ss_index_a );
	const sec_struc_planar_angles &planar_angles_b  = from_sec_struc_b.get_planar_angles_of_index( arg_to_ss_index_b );

	// If types of beta strands are different, return
	if (from_sec_struc_a.get_type() != from_sec_struc_b.get_type() || to_sec_struc_a.get_type() != to_sec_struc_b.get_type() ) {
		return 0;
	}

	const coord orig_from_to_vec_a       = calculate_inter_sec_struc_vector( arg_protein_a, arg_a_view_from_index, arg_to_ss_index_a );
	const coord orig_from_to_vec_b       = calculate_inter_sec_struc_vector( arg_protein_b, arg_b_view_from_index, arg_to_ss_index_b );
	const coord scaled_from_to_vec_a     = numeric_cast<double>(entry_querier::INTEGER_SCALING) * orig_from_to_vec_a;
	const coord scaled_from_to_vec_b     = numeric_cast<double>(entry_querier::INTEGER_SCALING) * orig_from_to_vec_b;
	const coord int_scaled_from_to_vec_a = int_cast_copy( scaled_from_to_vec_a );
	const coord int_scaled_from_to_vec_b = int_cast_copy( scaled_from_to_vec_b );

	// If types of helices are different, return
	const size_t  a_dist       = numeric_cast<size_t>(length(scaled_from_to_vec_a));
	const size_t  b_dist       = numeric_cast<size_t>(length(scaled_from_to_vec_b));
	const size_t  d_dist       = difference( a_dist,                                     b_dist                                     );
	const double  d_angle1     = difference( planar_angles_a.get_planar_angle_x(),       planar_angles_b.get_planar_angle_x()       );
	const double  d_angle2     = difference( planar_angles_a.get_planar_angle_minus_y(), planar_angles_b.get_planar_angle_minus_y() );
	const double  d_angle3     = difference( planar_angles_a.get_planar_angle_z(),       planar_angles_b.get_planar_angle_z()       );
	const size_t  mean_d_angle = numeric_cast<size_t>( ( d_angle1 + d_angle2 + d_angle3 ) / 3.0 );

	if ( ( from_sec_struc_a.get_type() == sec_struc_type::ALPHA_HELIX ) && ( from_sec_struc_b.get_type() == sec_struc_type::ALPHA_HELIX ) &&
	     (   to_sec_struc_a.get_type() == sec_struc_type::ALPHA_HELIX ) && (   to_sec_struc_b.get_type() == sec_struc_type::ALPHA_HELIX ) &&
	     ( d_dist < 15 ) &&
	     ( d_angle1 > 90 || d_angle2 > 90 || d_angle3 > 90 ) ) {
		return 0;
	}

	size_t s_vect = 0;
	size_t squared_distance = 0;
	if ( ( abs( int_scaled_from_to_vec_a.get_x() ) + abs( int_scaled_from_to_vec_a.get_y() ) + abs( int_scaled_from_to_vec_a.get_z() ) != 0.0 ) &&
	     ( abs( int_scaled_from_to_vec_b.get_x() ) + abs( int_scaled_from_to_vec_b.get_y() ) + abs( int_scaled_from_to_vec_b.get_z() )  != 0.0 ) ) {

		const ptrdiff_t x_diff         = numeric_cast<ptrdiff_t>(int_scaled_from_to_vec_a.get_x()) - numeric_cast<ptrdiff_t>(int_scaled_from_to_vec_b.get_x());
		const size_t    x_diff_squared = numeric_cast<size_t>( x_diff * x_diff );
		squared_distance               = x_diff_squared;

		if (squared_distance < sec_struc_querier::SEC_STRUC_MAX_DIST_SQ_CUTOFF) {

			const ptrdiff_t y_diff          = numeric_cast<ptrdiff_t>(int_scaled_from_to_vec_a.get_y()) - numeric_cast<ptrdiff_t>(int_scaled_from_to_vec_b.get_y());
			const size_t    y_diff_squared  = numeric_cast<size_t>( y_diff * y_diff );
			squared_distance               += y_diff_squared;

			if ( squared_distance < sec_struc_querier::SEC_STRUC_MAX_DIST_SQ_CUTOFF ) {

				const ptrdiff_t z_diff          = numeric_cast<ptrdiff_t>( int_scaled_from_to_vec_a.get_z()) - numeric_cast<ptrdiff_t>(int_scaled_from_to_vec_b.get_z() );
				const size_t    z_diff_squared  = numeric_cast<size_t>( z_diff * z_diff );
				squared_distance               += z_diff_squared;

				if ( squared_distance < sec_struc_querier::SEC_STRUC_MAX_DIST_SQ_CUTOFF ) {
					s_vect = sec_struc_querier::SEC_STRUC_A_VALUE / ( squared_distance + sec_struc_querier::SEC_STRUC_B_VALUE );
				}
			}
		}
	}

	// Score comparison of angles between secondary structures
	size_t s_angle = 0;
	if ( SEC_STRUC_PLANAR_W_ANGLE > 0 &&
	     ( planar_angles_a.get_planar_angle_x() != 0.0 ) && ( planar_angles_a.get_planar_angle_minus_y() != 0.0 ) && ( planar_angles_a.get_planar_angle_z() != 0.0 ) &&
	     ( planar_angles_b.get_planar_angle_x() != 0.0 ) && ( planar_angles_b.get_planar_angle_minus_y() != 0.0 ) && ( planar_angles_b.get_planar_angle_z() != 0.0 ) ) {
		s_angle = ( SEC_STRUC_PLANAR_W_ANGLE * SEC_STRUC_PLANAR_A_ANGLE) / (mean_d_angle + SEC_STRUC_PLANAR_B_ANGLE );
		if ( s_angle < SEC_STRUC_PLANAR_C_ANGLE ) {
			s_angle = 0;
		}
		if ( d_angle1 > 90 || d_angle2 > 90 || d_angle3 > 90.0 ) {
			s_angle = 0;
		}
	}

	return numeric_cast<score_type>( s_vect + s_angle );
}



/// \brief Calculate a normalised, log score by scoring comparison of vectors for each aligned pair along final path
///
/// Currently, this calls context_res() and context_sec() but is otherwise local as far as I can see
///
/// \todo Shift the final scoring to use doubles throughout.
///       This may require making context_res() and context_sec into templates so that
///       the alignment algorithm can continue using ints for a while longer.
///       This will also make it easy to test that the changes in final scores are very
///       small and that alignments are unaffected.
///       Once final scoring is done with doubles, it should be relatively easy to check
///       that making the alignment algorithm use doubles (or maybe floats) only improves final scores.
ssap_scores cath::calculate_log_score(const alignment     &arg_alignment,    ///< The alignment to score
                                      const protein       &arg_protein_a,    ///< The first protein
                                      const protein       &arg_protein_b,    ///< The second protein
                                      const entry_querier &arg_entry_querier ///< The entry_querier to query either residues or secondary structures
                                      ) {
	const size_t length_a = arg_entry_querier.get_length(arg_protein_a);
	const size_t length_b = arg_entry_querier.get_length(arg_protein_b);
	size_t       count    = 0;

	// Matrix in which to do the final score calculations, initialised to zeroes
	const size_t max_alignment_length( length_a + length_b + 10 );
	vector<vector<score_type> > final_score_matrix( max_alignment_length, vector<score_type>( max_alignment_length , 0 ) );

	// FOR EACH CELL ALONG FINAL PATH COMPARE VECTORS TO OTHER RESIDUES ON PATH
	for (size_t alignment_ctr_i = 0; alignment_ctr_i < arg_alignment.length(); ++alignment_ctr_i) {
		const bool i_has_both_posns = has_both_positions_of_index( arg_alignment, alignment_ctr_i );

		if (i_has_both_posns) {
			const aln_posn_type i_posn_a = get_a_offset_1_position_of_index( arg_alignment, alignment_ctr_i );
			const aln_posn_type i_posn_b = get_b_offset_1_position_of_index( arg_alignment, alignment_ctr_i );

			for (size_t alignment_ctr_j = 0; alignment_ctr_j < arg_alignment.length(); ++alignment_ctr_j) {
				const bool j_has_both_posns = has_both_positions_of_index( arg_alignment, alignment_ctr_j );

				if (j_has_both_posns) {
					const aln_posn_type j_posn_a = get_a_offset_1_position_of_index( arg_alignment, alignment_ctr_j );
					const aln_posn_type j_posn_b = get_b_offset_1_position_of_index( arg_alignment, alignment_ctr_j );

					// Score of distance from (
					//   j_posn_a seen from residue_i_a view
					// ) to (
					//   j_posn_b seen from residue_i_b view
					// )

					const bool are_comparable = arg_entry_querier.are_comparable__offset_1(
						arg_protein_a, arg_protein_b,
							 i_posn_a,      i_posn_b,
							 j_posn_a,      j_posn_b
					);
					if ( are_comparable ) {
						++count;
						const score_type pair_score = arg_entry_querier.distance_score__offset_1(
							arg_protein_a, arg_protein_b,
								 i_posn_a,      i_posn_b,
								 j_posn_a,      j_posn_b
						);
//						cerr << i_posn_a << "\t" << j_posn_a << "\t" << i_posn_b << "\t" << j_posn_b << "\t" << pair_score << endl;
						final_score_matrix[j_posn_b][j_posn_a] += pair_score;
					}
				}
			}
		}
	}

	// ACCUMULATE NEW SCORES ALONG FINAL PATH
	bool          prev_had_both_posns      = false;
	bool          is_first_with_both_posns = true;
	aln_posn_type prev_a_position          = 0;
	aln_posn_type prev_b_position          = 0;
	size_t        num_aligned_pairs        = 0;
	score_type    maxscore                 = 0;

//	cerr << "maxscore at start is            : " << maxscore << endl;

	for (size_t alignment_ctr = 0; alignment_ctr < arg_alignment.length(); ++alignment_ctr) {
		const bool has_both_posns = has_both_positions_of_index( arg_alignment, alignment_ctr );

		if ( has_both_posns ) {
			const aln_posn_type a_position = get_a_offset_1_position_of_index( arg_alignment, alignment_ctr );
			const aln_posn_type b_position = get_b_offset_1_position_of_index( arg_alignment, alignment_ctr );

			score_type gap = 0;
			if ( ! is_first_with_both_posns ) {
				if ( ! prev_had_both_posns || ( a_position != prev_a_position + 1 ) || ( b_position != prev_b_position + 1 ) ) {
					gap = numeric_cast<score_type>( get_gap_penalty( arg_entry_querier ) );
				}
			}
			maxscore += final_score_matrix[b_position][a_position] - gap;

			++num_aligned_pairs;
			is_first_with_both_posns = false;
		}

		prev_had_both_posns = has_both_posns;
		if ( prev_had_both_posns ) {
			prev_a_position = get_a_offset_1_position_of_index( arg_alignment, alignment_ctr );
			prev_b_position = get_b_offset_1_position_of_index( arg_alignment, alignment_ctr );
		}
	}

//	cerr << "maxscore after handling gaps is : " << maxscore << endl;

	maxscore = max(numeric_cast<score_type>(0), maxscore);

//	cerr << "maxscore after flooring at 0 is : " << maxscore << endl;
//	cerr << "count    is                     : " << count    << endl;

	const double final_score_scaling  = 1000.0;

	// (note that log() is base-e not base-10 (which is implemented through log10()))
	const double optimum_single_score = arg_entry_querier.optimum_single_score();
	const double max_log              = std::log( optimum_single_score * final_score_scaling );

	ssap_scores local_ssap_scores;

	/// A few notes on the final SSAP scores...
	///
	/// It is helpful to define the following:
	///  * x is the mean dynamic programming matrix score (ie the mean closeness score) divided by the maximum possible
	///  * k is the value by which this is multiplied before taking logs
	///    (which, in terms of the a and b from the original SSAP paper,
	///     is \f$ \frac{1000a}{b} \f$, eg 50,000 for residues and 133,333.33... for secondary structures)
	///
	/// Then the final score is \f$ 100 \left( \frac{\log{(xk)}}{\log{k}} \right) \f$, which is equal to
	/// \f$ 100 \left( \frac{\log{x}}{\log{k}} + 1 \right) \f$, which is equal to
	/// \f$ 100 \left( 1 + \log_k{x} \right) \f$
	///
	/// This means that:
	///  - when \f$ x = 1 \f$ the score is 100;
	///  - when \f$ x = \frac{1}{k} \f$ the score is 0 and
	///  - as x tends to 0 from above, the score tends to \f$ - \infty \f$.

	// CALCULATE AVERAGE SCORE (this doesn't appear to ever be used)
	if ( maxscore && count ) {
		const double scaled_avg_single_score   = numeric_cast<double>(maxscore) * final_score_scaling / numeric_cast<double>(count);
		const double final_score_over_compared = 100.0 * std::log(scaled_avg_single_score) / max_log;
		local_ssap_scores.set_ssap_score_over_compared(final_score_over_compared);
	}

	// CALCULATE GLOBAL SCORE 1 (this is normalised over the smallest protein)
	const size_t min_length              = min(length_a, length_b);
	const size_t num_comparable_over_min = num_comparable(arg_entry_querier, min_length);
	if ( maxscore && num_comparable_over_min > 0 ) {
		const double scaled_avg_single_score   = numeric_cast<double>(maxscore) * final_score_scaling / numeric_cast<double>(num_comparable_over_min);
		const double final_score_over_compared = 100.0 * std::log(scaled_avg_single_score) / max_log;
		local_ssap_scores.set_ssap_score_over_smaller(final_score_over_compared);
	}

	// CALCULATE GLOBAL SCORE 2 (this is normalised over the largest protein)
	const size_t max_length              = max(length_a, length_b);
	const size_t num_comparable_over_max = num_comparable(arg_entry_querier, max_length);
	if ( maxscore && num_comparable_over_max > 0 ) {
		const double scaled_avg_single_score   = numeric_cast<double>(maxscore) * final_score_scaling / numeric_cast<double>(num_comparable_over_max);
		const double final_score_over_compared = 100.0 * std::log(scaled_avg_single_score) / max_log;
		local_ssap_scores.set_ssap_score_over_larger(final_score_over_compared);
	}

	if ( arg_entry_querier.temp_hacky_is_residue() ) {
		if ( num_aligned_pairs > 0 && max_length ) {
			local_ssap_scores.set_percentage_aligned_pairs_over_larger(
				( 100.0 * numeric_cast<double>( num_aligned_pairs ) ) / numeric_cast<double>(max_length)
			);
		}
		local_ssap_scores.set_num_aligned_pairs( num_aligned_pairs );
	}

	// If considering residues, calculate the sequence identity and store it the ssap_scores object
	if ( arg_entry_querier.temp_hacky_is_residue() ) {
		local_ssap_scores.set_seq_id(calculate_sequence_identity(arg_alignment, arg_protein_a, arg_protein_b));
	}

//	cerr << "***************: Scores are\t"    << local_ssap_scores.get_ssap_score_over_larger();
//	cerr << "\t"              << local_ssap_scores.get_ssap_score_over_smaller();
//	cerr << "\t"              << local_ssap_scores.get_ssap_score_over_compared();
//	cerr << " for alignment " << arg_alignment;
//	cerr << endl;

	return local_ssap_scores;
}


/// \brief Calculate the sequence identity from an alignment and the two proteins
///
/// \todo Investigate an apparent bug in the existing code:
///       When SSAP-ing 2ij2B00/2odkA00, the original SSAP appears to output
///       an overlap and sequence identity from a different alignment from
///       the final one that's actually output.
///       (Probably also true of aligned residues but the answer's the same in both cases)
///       Important: be very clear about whether the problem is:
///         - the wrong alignment was output for the SSAP score that was chosen or
///         - the alignment is the correct one for the SSAP score and it's
///           just the sequence identity and overlap that are wrong.
double cath::calculate_sequence_identity(const alignment &arg_alignment, ///< The alignment from which the sequence identity should be calculated
                                         const protein   &arg_protein_a, ///< The first  protein involved in the alignment
                                         const protein   &arg_protein_b  ///< The second protein involved in the alignment
                                         ) {
	// For each residue pair in the alignment, increment num_amino_acid_matches if the amino acids match
	size_t num_amino_acid_matches = 0;
	for (size_t alignment_ctr = 0; alignment_ctr < arg_alignment.length(); ++alignment_ctr) {
		if (has_both_positions_of_index(arg_alignment, alignment_ctr)) {
			const residue &residue_a = get_a_residue_cref_of_index(arg_alignment, arg_protein_a, alignment_ctr);
			const residue &residue_b = get_b_residue_cref_of_index(arg_alignment, arg_protein_b, alignment_ctr);
			if (residue_a.get_amino_acid() == residue_b.get_amino_acid()) {
				++num_amino_acid_matches;
			}
		}
	}

	// Return the count of equal amino acid types as a percentage of the shorter length
	const double min_length = numeric_cast<double>( min( arg_protein_a.get_length(), arg_protein_b.get_length() ) );
	return 100.0 * numeric_cast<double>( num_amino_acid_matches ) / min_length;
}


/// \brief TODOCUMENT
bool cath::save_ssap_scores(const alignment               &arg_alignment,    ///< The alignment for which scores should be output
                            const protein                 &arg_protein_a,    ///< The first protein
                            const protein                 &arg_protein_b,    ///< The second protein
                            const ssap_scores             &arg_ssap_scores,  ///< The scores to be output
                            const old_ssap_options_block  &arg_ssap_options, ///< The old_ssap_options_block to specify how things should be done
                            const data_dirs_spec          &arg_data_dirs     ///< The data directories from which data should be read
                            ) {
	BOOST_LOG_TRIVIAL( info ) << "Function: save_ssap_scores";
	
	// Select get_ssap_score_over_smaller if a local score is needed, or get_ssap_score_over_larger otherwise
	const double select_score = arg_ssap_options.get_use_local_ssap_score() ? arg_ssap_scores.get_ssap_score_over_smaller()
	                                                             : arg_ssap_scores.get_ssap_score_over_larger();
	
	// Do not output files if score is not high enough
	const bool score_is_high_enough = (select_score >= arg_ssap_options.get_min_score_for_writing_files());
	
	// Get RMSD value
	const pair<size_t, double> num_superposed_and_rmsd = superpose(
		arg_protein_a,
		arg_protein_b,
		arg_alignment,
		arg_ssap_options,
		arg_data_dirs,
		score_is_high_enough
	);
	const size_t &num_superposed = num_superposed_and_rmsd.first;
	const double &rmsd           = num_superposed_and_rmsd.second;

	// For Fast SSAP
	if (global_run_counter == 1) {
		snprintf(
			global_ssap_line1,
			SSAP_LINE_LENGTH - 1,
			"%6s  %6s %4zu %4zu %6.2f %4zu %4zu %4zu %6.2f",
			arg_protein_a.get_title().c_str(),
			arg_protein_b.get_title().c_str(),
			arg_protein_a.get_length(),
			arg_protein_b.get_length(),
			select_score,
			arg_ssap_scores.get_num_aligned_pairs(),
			numeric_cast<size_t>( arg_ssap_scores.get_percentage_aligned_pairs_over_larger() ),
			numeric_cast<size_t>( arg_ssap_scores.get_seq_id() ),
			rmsd
		);
		// If the cutoff for superposition was above the default, then output the number of aligned residue pairs
		// used in the superposition
		if (arg_ssap_options.get_min_score_for_superposition() > common_residue_select_min_score_policy::MIN_CUTOFF) {
			const string temp_prev_global_ssap_line1(global_ssap_line1);
			snprintf( global_ssap_line1, SSAP_LINE_LENGTH - 1, "%s %4zu", temp_prev_global_ssap_line1.c_str(), num_superposed );
		}
				
		global_ssap_score1 = select_score;
	}
	// For Slow SSAP
	else if (global_run_counter == 2) {
		snprintf(
			global_ssap_line2,
			SSAP_LINE_LENGTH - 1,
			"%6s  %6s %4zu %4zu %6.2f %4zu %4zu %4zu %6.2f",
			arg_protein_a.get_title().c_str(),
			arg_protein_b.get_title().c_str(),
			arg_protein_a.get_length(),
			arg_protein_b.get_length(),
			select_score,
			arg_ssap_scores.get_num_aligned_pairs(),
			numeric_cast<size_t>( arg_ssap_scores.get_percentage_aligned_pairs_over_larger() ),
			numeric_cast<size_t>( arg_ssap_scores.get_seq_id() ),
			rmsd
		);
		// If the cutoff for superposition was above the default, then output the number of aligned residue pairs
		// used in the superposition
		if (arg_ssap_options.get_min_score_for_superposition() > common_residue_select_min_score_policy::MIN_CUTOFF) {
			const string temp_prev_global_ssap_line2(global_ssap_line2);
			snprintf( global_ssap_line2, SSAP_LINE_LENGTH - 1, "%s %4zu", temp_prev_global_ssap_line2.c_str(), num_superposed );
		}

		global_ssap_score2 = select_score;	
	}

	return score_is_high_enough;
}


/// \brief TODOCUMENT
void cath::save_zero_scores(const protein &arg_protein_a, ///< The first protein
                            const protein &arg_protein_b  ///< The second protein
                            ) {
		BOOST_LOG_TRIVIAL( info ) << "Function: save_zero_scores()";
	
	if (global_run_counter == 1) {
		snprintf(
			global_ssap_line1,
			SSAP_LINE_LENGTH - 1,
			"%6s  %6s %4d %4d %6.2f %4d %4d %4d %6.2f",
			arg_protein_a.get_title().c_str(),
			arg_protein_b.get_title().c_str(),
			0,
			0,
			0.0,
			0,
			0,
			0,
			0.0
		);
		
		global_ssap_score1 = 0.0;
	}
	else if (global_run_counter == 2) {
		snprintf(
			global_ssap_line2,
			SSAP_LINE_LENGTH - 1,
			"%6s  %6s %4d %4d %6.2f %4d %4d %4d %6.2f %6.2f",
			arg_protein_a.get_title().c_str(),
			arg_protein_b.get_title().c_str(),
			0,
			0,
			0.0,
			0,
			0,
			0,
			0.0,
			0.0
		);
	
		global_ssap_score2 = 0.0;
	}
}


/// \brief TODOCUMENT
void cath::print_ssap_scores(ostream            &arg_os,              ///< TODOCUMENT
                             const double       &arg_ssap_score_1,    ///< TODOCUMENT
                             const double       &arg_ssap_score_2,    ///< TODOCUMENT
                             const string       &arg_ssap_line1,      ///< TODOCUMENT
                             const string       &arg_ssap_line2,      ///< TODOCUMENT
                             const ptrdiff_t    &arg_run_counter,     ///< TODOCUMENT
                             const bool         &arg_write_all_scores ///< TODOCUMENT
                             ) {
	if (arg_run_counter == 1) {
		arg_os << arg_ssap_line1 << "\n";
	}
	else if (arg_run_counter == 2) {
		if (arg_write_all_scores) {
			arg_os << arg_ssap_line1 << "\n";
			arg_os << arg_ssap_line2 << "\n";
		}
		else {
			const string best_line = (arg_ssap_score_2 >= arg_ssap_score_1) ? arg_ssap_line2 : arg_ssap_line1;
			arg_os << best_line << "\n";
		}
	}
	else {
		BOOST_LOG_TRIVIAL( warning ) << "There's something strange in your neighbourhood";
	}
}


/// \brief Superpose two structures based on an alignment between them
///
/// AFAIK, this currently reads from the following globals:
///  - global_upper_score_matrix
///  - global_window
/// but is otherwise local.
size_doub_pair cath::superpose(const protein                 &arg_protein_a,           ///< Coordinates for first structure
                               const protein                 &arg_protein_b,           ///< Coordinates for second structure
                               const alignment               &arg_alignment,           ///< The alignment to determine which residues should be as close as possible to which
                               const old_ssap_options_block  &arg_ssap_options,        ///< The old_ssap_options_block to specify how things should be done
                               const data_dirs_spec          &arg_data_dirs,           ///< The data directories from which data should be read
                               const bool                    &arg_score_is_high_enough ///< Whether the score is high enough to justify outputting files
                               ) {
	const auto common_coords = alignment_coord_extractor::get_common_coords(
		arg_alignment,
		arg_protein_a,
		arg_protein_b,
		common_residue_select_min_score_policy( arg_ssap_options.get_min_score_for_superposition() )
	);
	const auto &coord_list_1         = common_coords.first;
	const auto &coord_list_2         = common_coords.second;
	const auto  num_superposed_pairs = coord_list_1.size();
	assert( coord_list_1.size() == coord_list_2.size() );

	// Perform the superposition
	const auto my_superposition = create_pairwise_superposition( coord_list_1, coord_list_2, true, -centre_of_gravity( coord_list_1 ) );
	const auto rmsd             = calc_pairwise_superposition_rmsd( coord_list_1, coord_list_2 );

	// If an XML superposition file has been requested, write one
	if ( arg_score_is_high_enough && arg_ssap_options.get_write_xml_sup() ) {
		const auto xml_outname = arg_protein_a.get_title() + arg_protein_b.get_title() + ".superpose.xml";
		write_xml_sup_filename( my_superposition, xml_outname, { arg_protein_a.get_title(), arg_protein_b.get_title() } );
	}

	// If a sup file has been requested, write one
	if ( arg_score_is_high_enough && has_superposition_dir( arg_ssap_options ) ) {
		const auto pdb1_filename   = find_file( arg_data_dirs, data_file::PDB,  arg_protein_a.get_title() );
		const auto pdb2_filename   = find_file( arg_data_dirs, data_file::PDB,  arg_protein_b.get_title() );
		const auto sup_file_suffix = string( arg_ssap_options.get_write_rasmol_script() ? ".rasc" : ".sup" );
		const auto basename        = ( arg_protein_a.get_title() + arg_protein_b.get_title() + sup_file_suffix );
		const auto outname         = get_superposition_dir( arg_ssap_options ) / basename;

		// Write the superposed PDB file
		write_superposed_pdb_from_files(
			my_superposition,
			outname,
			{ pdb1_filename, pdb2_filename },
			arg_ssap_options.get_write_rasmol_script(),
			true
		);
	}

	// Return the calculated RMSD
	return { num_superposed_pairs, rmsd };
}


/// \brief Prints alignment of structures and score matrices
///
/// A fairly messy subroutine that appears to have quite a lot of interaction with various global variables.
///
/// At some point it decides whether an alignment should be printed and if so does it by calling print_aln().
ssap_scores cath::plot_aln(const protein                 &arg_protein_a,     ///< The first protein
                           const protein                 &arg_protein_b,     ///< The second protein
                           const size_t                  &arg_pass,          ///< The pass of this comparison (where the second typically refines the alignment generated by the first)
                           const entry_querier           &arg_entry_querier, ///< The entry_querier to query either residues or secondary structures
                           const alignment               &arg_alignment,     ///< The alignment to plot
                           const old_ssap_options_block  &arg_ssap_options,  ///< The old_ssap_options_block to specify how things should be done
                           const data_dirs_spec          &arg_data_dirs      ///< The data directories from which data should be read
                           ) {
	const bool res_not_ss__hacky = arg_entry_querier.temp_hacky_is_residue();
	if (res_not_ss__hacky && arg_pass != 2) {
		return ssap_scores();
	}

	const ssap_scores local_ssap_scores = calculate_log_score(
		arg_alignment,
		arg_protein_a,
		arg_protein_b,
		arg_entry_querier
	);

	// Score secondary structure alignment
	if (!res_not_ss__hacky) {
		return ssap_scores();
	}

	// Score and print residue alignment
	global_res_score = true;

	// Select global1 if a local score is required
	const double select_score = arg_ssap_options.get_use_local_ssap_score()
								? local_ssap_scores.get_ssap_score_over_smaller()
								: local_ssap_scores.get_ssap_score_over_larger();

	// Changed print_ssap_scores to save_ssap_scores (v1.14 JEB)
	const bool score_is_high_enough = save_ssap_scores(arg_alignment, arg_protein_a, arg_protein_b, local_ssap_scores, arg_ssap_options, arg_data_dirs);

	// global_score_run1 & global_score_run2 are used to determine whether second alignment should be written out (JEB 12.09.2002 v1.10)
	if (global_doing_fast_ssap) {
		global_score_run1 = select_score;
	}
	else {
		global_score_run2 = select_score;
	}

	BOOST_LOG_TRIVIAL( info ) << "Function: plot_aln:  score_run1 = " << fixed << setprecision(3) << global_score_run1;
	BOOST_LOG_TRIVIAL( info ) << "Function: plot_aln:  score_run2 = " << fixed << setprecision(3) << global_score_run2;
	BOOST_LOG_TRIVIAL( info ) << "Function: plot_aln:  r_fast     = " << boolalpha << global_doing_fast_ssap;;

	// A decision is made here about whether to write an alignment file, based on
	// the score of the alignment. However, this is inconsistent with save_ssap_scores
	// and hence some alignments may not be written when they have a SSAP score. This
	// appears to only affect fairly bad alignments (ssap score < 50)
	if (global_supaln) {
		double out_score = -1.0;

		// Prints SSAP alignment
		// Always for fast run and only for slow run if score is better than for fast run
		if (global_doing_fast_ssap) {
			out_score = global_score_run1;
		}
		if (!global_doing_fast_ssap && global_score_run2 > global_score_run1) {
			out_score = global_score_run2;
		}
		if (out_score > -1.0) {
			BOOST_LOG_TRIVIAL( info ) << "Function: plot_aln: printing alignment (r_fast == 1) || (!r_fast && score_run2 > score_run1)";

			// Prevent writing if score isn't high enough
			if (score_is_high_enough != (out_score >= arg_ssap_options.get_min_score_for_writing_files())) {
				BOOST_THROW_EXCEPTION(out_of_range_exception("Code is inconsistent about what score is high enough"));
			}

			if (score_is_high_enough) {
				BOOST_LOG_TRIVIAL( info ) << "Function: print_aln";
				const path alignment_out_file = path(arg_ssap_options.get_alignment_dir()) / (arg_protein_a.get_title() + arg_protein_b.get_title() + ".list");
				write_alignment_as_cath_ssap_legacy_format(alignment_out_file, arg_alignment, arg_protein_a, arg_protein_b);
			}
		}
	}

	return local_ssap_scores;
}
