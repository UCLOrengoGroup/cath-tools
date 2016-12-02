/// \file
/// \brief The alignment_input_options_block class definitions

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

#include "alignment_input_options_block.hpp"

#include <boost/optional.hpp>

#include "common/clone/make_uptr_clone.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace std::literals::string_literals;

using boost::filesystem::path;
using boost::none;
using boost::program_options::bool_switch;
using boost::program_options::options_description;
using boost::program_options::value;
using boost::program_options::variables_map;
using std::string;
using std::unique_ptr;

/// \brief The option name for whether to align based on matching residue names
const string alignment_input_options_block::PO_RES_NAME_ALIGN    ( "res-name-align"     );

/// \brief The option name for a file from which to read a FASTA alignment
const string alignment_input_options_block::PO_FASTA_ALIGN_INFILE( "fasta-aln-infile"    );

/// \brief The option name for a file from which to read a legacy-SSAP-format alignment
const string alignment_input_options_block::PO_SSAP_ALIGN_INFILE ( "ssap-aln-infile"    );

/// \brief The option name for a file from which to read a CORA alignment
const string alignment_input_options_block::PO_CORA_ALIGN_INFILE ( "cora-aln-infile"    );

/// \brief The option name for a file from which to read SSAP-scores format data to use to attempt to glue pairwise alignments together
const string alignment_input_options_block::PO_SSAP_SCORE_INFILE ( "ssap-scores-infile" );

/// \brief A standard do_clone method.
unique_ptr<options_block> alignment_input_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string alignment_input_options_block::do_get_block_name() const {
	return "Alignment source";
}

/// \brief Add this block's options to the provided options_description
void alignment_input_options_block::do_add_visible_options_to_description(options_description &arg_desc ///< The options_description to which the options are added
                                                                  ) {
	const string file_varname { "<file>" };

	const auto residue_name_align_notifier   = [&] (const bool &x) { the_alignment_input_spec.set_residue_name_align  ( x ); };
	const auto fasta_alignment_file_notifier = [&] (const path &x) { the_alignment_input_spec.set_fasta_alignment_file( x ); };
	const auto ssap_alignment_file_notifier  = [&] (const path &x) { the_alignment_input_spec.set_ssap_alignment_file ( x ); };
	const auto cora_alignment_file_notifier  = [&] (const path &x) { the_alignment_input_spec.set_cora_alignment_file ( x ); };
	const auto ssap_scores_file_notifier     = [&] (const path &x) { the_alignment_input_spec.set_ssap_scores_file    ( x ); };

	arg_desc.add_options()
		(
			PO_RES_NAME_ALIGN.c_str(),
			bool_switch()
				->notifier     ( residue_name_align_notifier                      )
				->default_value( alignment_input_spec::DEFAULT_RESIDUE_NAME_ALIGN ),
			"Align residues by simply matching their names (numbers+insert)\n(for multiple models of the same structure)"
		)
		(
			PO_FASTA_ALIGN_INFILE.c_str(),
			value<path>()
				->value_name( file_varname                  )
				->notifier  ( fasta_alignment_file_notifier ),
			( "Read FASTA alignment from file " + file_varname ).c_str()
		)
		(
			PO_SSAP_ALIGN_INFILE.c_str(),
			value<path>()
				->value_name( file_varname                 )
				->notifier  ( ssap_alignment_file_notifier ),
			( "Read SSAP alignment from file " + file_varname ).c_str()
		)
		(
			PO_CORA_ALIGN_INFILE.c_str(),
			value<path>()
				->value_name( file_varname                 )
				->notifier  ( cora_alignment_file_notifier ),
			( "Read CORA alignment from file " + file_varname ).c_str()
		)
		(
			PO_SSAP_SCORE_INFILE.c_str(),
			value<path>()
				->value_name( file_varname              )
				->notifier  ( ssap_scores_file_notifier ),
			( "Read SSAP scores from file " + file_varname + "\nAssumes all .list alignment files in same directory" ).c_str()
		);

	static_assert( ! alignment_input_spec::DEFAULT_RESIDUE_NAME_ALIGN,
		"If alignment_input_spec::DEFAULT_RESIDUE_NAME_ALIGN isn't false, it might mess up the bool switch in here" );
}

/// \brief TODOCUMENT
str_opt alignment_input_options_block::do_invalid_string(const variables_map &/*arg_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                         ) const {
	if ( get_num_acquirers( *this ) > 1 ) {
		return "Cannot specify more than one alignment input"s;
	}
	if ( ! the_alignment_input_spec.get_fasta_alignment_file().empty() && ! is_acceptable_input_file( the_alignment_input_spec.get_fasta_alignment_file()    ) ) {
		return "FASTA alignment file " + the_alignment_input_spec.get_ssap_alignment_file().string() + " is not a valid input file";
	}
	if ( ! the_alignment_input_spec.get_ssap_alignment_file().empty()  && ! is_acceptable_input_file( the_alignment_input_spec.get_ssap_alignment_file()    ) ) {
		return "SSAP alignment file " + the_alignment_input_spec.get_ssap_alignment_file().string() + " is not a valid input file";
	}
	if ( ! the_alignment_input_spec.get_cora_alignment_file().empty()  && ! is_acceptable_input_file( the_alignment_input_spec.get_cora_alignment_file()    ) ) {
		return "CORA alignment file " + the_alignment_input_spec.get_cora_alignment_file().string() + " is not a valid input file";
	}
	if ( ! the_alignment_input_spec.get_ssap_scores_file().empty()     && ! is_acceptable_input_file( the_alignment_input_spec.get_ssap_scores_file(), true ) ) {
		return "SSAP scores file "    + the_alignment_input_spec.get_ssap_scores_file().string()    + " is not a valid input file";
	}
	return none;
}

/// \brief TODOCUMENT
const alignment_input_spec & alignment_input_options_block::get_alignment_input_spec() const {
	return the_alignment_input_spec;
}

/// \brief Get the number of alignment_acquirer objects implied by the alignment_input_spec in the specified alignment_input_options_block
///
/// \relates alignment_input_options_block
size_t cath::opts::get_num_acquirers(const alignment_input_options_block &arg_alignment_input_options_block ///< The alignment_input_options_block to query
                                     ) {
	return get_num_acquirers( arg_alignment_input_options_block.get_alignment_input_spec() );
}
