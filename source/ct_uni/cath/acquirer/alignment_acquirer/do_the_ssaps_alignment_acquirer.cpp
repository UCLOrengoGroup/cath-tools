/// \file
/// \brief The do_the_ssaps_alignment_acquirer class definitions

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

#include "do_the_ssaps_alignment_acquirer.hpp"

#include <filesystem>

#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/format.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/irange.hpp>

#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>

#include "cath/acquirer/alignment_acquirer/ssap_scores_file_alignment_acquirer.hpp"
#include "cath/alignment/alignment.hpp"
#include "cath/chopping/chopping_format/sillitoe_chopping_format.hpp"
#include "cath/chopping/domain/domain.hpp"
#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/common/exception/not_implemented_exception.hpp"
#include "cath/common/exception/runtime_error_exception.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/common/file/slurp.hpp"
#include "cath/common/file/spew.hpp"
#include "cath/common/matrix/matrix_index.hpp"
#include "cath/common/test_or_exe_run_mode.hpp"
#include "cath/file/strucs_context.hpp"
#include "cath/ssap/options/cath_ssap_options.hpp"
#include "cath/ssap/ssap.hpp"
#include "cath/superposition/options/align_regions_options_block.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::chop;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::opts;

using ::boost::adaptors::transformed;
using ::boost::algorithm::join;
using ::boost::algorithm::trim_right_copy;
using ::boost::format;
using ::boost::irange;
using ::std::filesystem::path;
using ::std::filesystem::temp_directory_path;
using ::std::max;
using ::std::ofstream;
using ::std::pair;
using ::std::string;
using ::std::unique_ptr;

/// \brief A standard do_clone method.
unique_ptr<alignment_acquirer> do_the_ssaps_alignment_acquirer::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Specify that this does require backbone-complete input
bool do_the_ssaps_alignment_acquirer::do_requires_backbone_complete_input() const {
	return true;
}

/// \brief Run the necessary cath-ssaps and then use them to get the alignment and spanning tree
pair<alignment, size_size_pair_vec> do_the_ssaps_alignment_acquirer::do_get_alignment_and_spanning_tree(const strucs_context &prm_strucs_context, ///< The details of the structures for which the alignment and spanning tree is required
                                                                                                        const align_refining &prm_align_refining  ///< How much refining should be done to the alignment
                                                                                                        ) const {
	using ::std::to_string;

	// Get the directory in which the cath-ssaps should be done
	const path ssaps_dir = get_directory_of_joy().value_or(
		make_temp_dir_for_doing_ssaps( prm_strucs_context )
	);

	// Ensure the directory exists
	if ( ! exists( ssaps_dir ) ) {
		::spdlog::info( "About to create directory {}", ssaps_dir );
		if ( ! create_directories( ssaps_dir ) ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception(
				"Unable to create directory "
				+ ssaps_dir.string()
				+ " for temporary cath-ssaps"
			));
		}
	}

	// Perform any necessary cath-ssaps
	path_vec scores_files;
	const size_t num_strucs        = size( prm_strucs_context );
	const size_t num_ssaps         = ( num_strucs * ( max( 1_z, num_strucs ) - 1_z ) ) / 2_z;
	const string num_str           = to_string( num_ssaps );
	const size_t num_str_width     = num_str.length();
	const string num_str_width_str = "%" + to_string( num_str_width ) + "d";
	::spdlog::info( "About to check for and possibly run {} cath-ssaps in directory {}", num_ssaps, ssaps_dir );
	for (const size_t &struc_1_index : indices( num_strucs ) ) {
		for (const size_t &struc_2_index : irange( struc_1_index + 1, num_strucs ) ) {

			const size_t comp_ctr_offset_1 = get_zero_index_of_strict_upper_half_matrix(
				struc_1_index,
				struc_2_index,
				num_strucs
			) + 1;
			const string progress_str = ( format( num_str_width_str ) % comp_ctr_offset_1 ).str()
				+ "/"
				+ num_str;

			// \TODO: Abstract out functions for making these standard file names
			const string id_1   = get_domain_or_specified_or_name_from_acq( prm_strucs_context.get_name_sets()[ struc_1_index ] );
			const string id_2   = get_domain_or_specified_or_name_from_acq( prm_strucs_context.get_name_sets()[ struc_2_index ] );
			const path   scores_file = ssaps_dir / ( id_1 + id_2 + ".scores" );
			const path   alnmnt_file = ssaps_dir / ( id_1 + id_2 + ".list"   );

			scores_files.push_back( scores_file );

			if (   ! exists( scores_file ) || ! exists( alnmnt_file )
				|| is_empty( scores_file ) || is_empty( alnmnt_file ) ) {

				// \TODO Tighten up the interface with the ssap/ssap.cpp code here
				//
				// In particular, the name_set will know the location of the file that each
				// PDB was loaded from (if it was loaded from a PDB) so cath-ssap should have
				// a --pdb-infile option and it should be explicitly used here.
				str_vec cath_ssap_args{
					cath_ssap_options::PROGRAM_NAME,
					"--" + old_ssap_options_block::PO_ALIGN_DIR, ssaps_dir.string()
				};
				for (const size_t &index : { struc_1_index, struc_2_index } ) {
					const auto       &name_set   = prm_strucs_context.get_name_sets()[ index ];
					const domain_opt  opt_domain = get_domain_opt_of_index( prm_strucs_context, index );
					cath_ssap_args.push_back( name_set.get_name_from_acq() );
					if ( opt_domain ) {
						cath_ssap_args.push_back( "--" + align_regions_options_block::PO_ALN_REGIONS );
						cath_ssap_args.push_back( sillitoe_chopping_format{}.write_domain( *opt_domain ) );
					}
				}
				::spdlog::info( "[{}] Running : {} and writing scores to {}",
				                progress_str,
				                join( cath_ssap_args, " " ),
				                scores_file.string() );
				ofstream out_scores_ofstream = open_ofstream( scores_file );
				run_ssap(
					make_and_parse_options<cath_ssap_options>(
						cath_ssap_args,
						( run_mode_flag::value == run_mode::TEST )
							? parse_sources::CMND_LINE_ONLY
							: parse_sources::CMND_ENV_AND_FILE
					),
					std::cout,
					std::cerr,
					ostream_ref( out_scores_ofstream )
				);
				out_scores_ofstream.close();
			}
			else {
				::spdlog::info( "[{}] Skipping {} versus {} - non-empty data files already exist", progress_str, id_1, id_2 );
			}
		}
	}

	// Concatenate the scores files into one big file
	const path scores_file = ssaps_dir / "ssap_scores";
	spew(
		scores_file,
		join(
			scores_files
				| transformed( [] (const path &x) { return trim_right_copy( slurp( x ) ); } ),
			"\n"
		)
	);

	// Use the ssap_scores_file_alignment_acquirer on the directory of data
	return ssap_scores_file_alignment_acquirer{ scores_file }
		.get_alignment_and_spanning_tree( prm_strucs_context, prm_align_refining );
}

/// \brief Ctor for do_the_ssaps_alignment_acquirer
do_the_ssaps_alignment_acquirer::do_the_ssaps_alignment_acquirer(const path_opt &prm_directory_of_joy ///< The directory in which the cath-ssaps should be performed
                                                                 ) : directory_of_joy { prm_directory_of_joy } {
}

/// \brief Getter for the directory in which the cath-ssaps should be performed
const path_opt & do_the_ssaps_alignment_acquirer::get_directory_of_joy() const {
	return directory_of_joy;
}

/// \brief Make a path with a temporary directory in which cath-ssaps for the specified strucs_context should be done
///
/// This aims to create a temporary directory that is likely to persist across multiple identical runs
/// but very unlikely to clash amongst different inputs
///
/// This doesn't actually create the directory
path do_the_ssaps_alignment_acquirer::make_temp_dir_for_doing_ssaps(const strucs_context &prm_strucs_context ///< The strucs_context for which the directory should make a cath-ssaps temporary directory
                                                                    ) {
	return temp_directory_path() / (
		  "cath-tools.tmp."
		+ ( boost::format( R"(%|08X|)" ) % non_crypto_hash_copy( 0, prm_strucs_context ) ).str()
		+ "_"
		+ ( boost::format( R"(%|08X|)" ) % non_crypto_hash_copy( 1, prm_strucs_context ) ).str()
	);
}
