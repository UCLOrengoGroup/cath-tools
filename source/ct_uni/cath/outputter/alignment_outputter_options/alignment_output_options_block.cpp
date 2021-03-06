/// \file
/// \brief The alignment_output_options_block class definitions

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

#include "alignment_output_options_block.hpp"

#include <filesystem>
#include <iostream>
#include <string_view>

#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/display/display_colourer/display_colourer_alignment.hpp"
#include "cath/display_colour/display_colour.hpp"
#include "cath/outputter/alignment_outputter/alignment_outputter_list.hpp"
#include "cath/outputter/alignment_outputter/cath_aln_ostream_alignment_outputter.hpp"
#include "cath/outputter/alignment_outputter/fasta_ostream_alignment_outputter.hpp"
#include "cath/outputter/alignment_outputter/file_alignment_outputter.hpp"
#include "cath/outputter/alignment_outputter/html_ostream_alignment_outputter.hpp"
#include "cath/outputter/alignment_outputter/ssap_ostream_alignment_outputter.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::std;

using ::boost::program_options::bool_switch;
using ::boost::program_options::options_description;
using ::boost::program_options::value;
using ::boost::program_options::variables_map;
using ::std::filesystem::path;
using ::std::nullopt;
using ::std::string_view;

constexpr string_view PO_ALN_TO_CATH_ALN_FILE   ( "aln-to-cath-aln-file"   );
constexpr string_view PO_ALN_TO_CATH_ALN_STDOUT ( "aln-to-cath-aln-stdout" );
constexpr string_view PO_ALN_TO_FASTA_FILE      ( "aln-to-fasta-file"      );
constexpr string_view PO_ALN_TO_FASTA_STDOUT    ( "aln-to-fasta-stdout"    );
constexpr string_view PO_ALN_TO_SSAP_FILE       ( "aln-to-ssap-file"       );
constexpr string_view PO_ALN_TO_SSAP_STDOUT     ( "aln-to-ssap-stdout"     );
constexpr string_view PO_ALN_TO_HTML_FILE       ( "aln-to-html-file"       );
constexpr string_view PO_ALN_TO_HTML_STDOUT     ( "aln-to-html-stdout"     );

/// \brief A standard do_clone method.
unique_ptr<options_block> alignment_output_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string alignment_output_options_block::do_get_block_name() const {
	return "Alignment output";
}

/// \brief Add this block's options to the provided options_description
void alignment_output_options_block::do_add_visible_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                           const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                           ) {
	prm_desc.add_options()
		( string( PO_ALN_TO_CATH_ALN_FILE ).c_str(),   value<path>( &aln_to_cath_aln_file   ),                       "[EXPERIMENTAL] Write the alignment to a CATH alignment file"           )
		( string( PO_ALN_TO_CATH_ALN_STDOUT ).c_str(), bool_switch( &aln_to_cath_aln_stdout )->default_value(false), "[EXPERIMENTAL] Print the alignment to stdout in CATH alignment format" )
		( string( PO_ALN_TO_FASTA_FILE ).c_str(),      value<path>( &aln_to_fasta_file      ),                       "Write the alignment to a FASTA file"                                   )
		( string( PO_ALN_TO_FASTA_STDOUT ).c_str(),    bool_switch( &aln_to_fasta_stdout    )->default_value(false), "Print the alignment to stdout in FASTA format"                         )
		( string( PO_ALN_TO_SSAP_FILE ).c_str(),       value<path>( &aln_to_ssap_file       ),                       "Write the alignment to a SSAP file"                                    )
		( string( PO_ALN_TO_SSAP_STDOUT ).c_str(),     bool_switch( &aln_to_ssap_stdout     )->default_value(false), "Print the alignment to stdout as SSAP"                                 )
		( string( PO_ALN_TO_HTML_FILE ).c_str(),       value<path>( &aln_to_html_file       ),                       "Write the alignment to a HTML file"                                    )
		( string( PO_ALN_TO_HTML_STDOUT ).c_str(),     bool_switch( &aln_to_html_stdout     )->default_value(false), "Print the alignment to stdout as HTML"                                 );
}

/// \brief TODOCUMENT
str_opt alignment_output_options_block::do_invalid_string(const variables_map &/*prm_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                          ) const {
	if ( ! get_aln_to_cath_aln_file().empty() && ! is_acceptable_output_file( get_aln_to_cath_aln_file() ) ) {
		return "Not a valid alignment CATH alignment output file :\"" + get_aln_to_cath_aln_file().string() + "\"";
	}
	if ( ! get_aln_to_fasta_file().empty() && ! is_acceptable_output_file( get_aln_to_fasta_file() ) ) {
		return "Not a valid alignment FASTA output file :\"" + get_aln_to_fasta_file().string() + "\"";
	}
	if ( ! get_aln_to_ssap_file().empty() && ! is_acceptable_output_file( get_aln_to_ssap_file() ) ) {
		return "Not a valid alignment SSAP output file :\"" + get_aln_to_ssap_file().string() + "\"";
	}
	if ( ! get_aln_to_html_file().empty() && ! is_acceptable_output_file( get_aln_to_html_file() ) ) {
		return "Not a valid alignment HTML output file :\"" + get_aln_to_html_file().string() + "\"";
	}
	if ( get_aln_to_html_stdout() && get_aln_to_fasta_stdout() ) {
		return string( "Cannot output alignment to stdout as both HTML and FASTA" );
	}

	return nullopt;
}

/// \brief Return all options names for this block
str_view_vec alignment_output_options_block::do_get_all_options_names() const {
	return {
		PO_ALN_TO_CATH_ALN_FILE,
		PO_ALN_TO_CATH_ALN_STDOUT,
		PO_ALN_TO_FASTA_FILE,
		PO_ALN_TO_FASTA_STDOUT,
		PO_ALN_TO_SSAP_FILE,
		PO_ALN_TO_SSAP_STDOUT,
		PO_ALN_TO_HTML_FILE,
		PO_ALN_TO_HTML_STDOUT,
	};
}

/// \brief TODOCUMENT
path alignment_output_options_block::get_aln_to_cath_aln_file() const {
	return aln_to_cath_aln_file;
}

/// \brief TODOCUMENT
bool alignment_output_options_block::get_aln_to_cath_aln_stdout() const {
	return aln_to_cath_aln_stdout;
}

/// \brief TODOCUMENT
path alignment_output_options_block::get_aln_to_fasta_file() const {
	return aln_to_fasta_file;
}

/// \brief TODOCUMENT
bool alignment_output_options_block::get_aln_to_fasta_stdout() const {
	return aln_to_fasta_stdout;
}

/// \brief TODOCUMENT
path alignment_output_options_block::get_aln_to_ssap_file() const {
	return aln_to_ssap_file;
}

/// \brief TODOCUMENT
bool alignment_output_options_block::get_aln_to_ssap_stdout() const {
	return aln_to_ssap_stdout;
}

/// \brief TODOCUMENT
path alignment_output_options_block::get_aln_to_html_file() const {
	return aln_to_html_file;
}

/// \brief TODOCUMENT
bool alignment_output_options_block::get_aln_to_html_stdout() const {
	return aln_to_html_stdout;
}

/// \brief TODOCUMENT
alignment_outputter_list alignment_output_options_block::get_alignment_outputters(const display_spec &prm_display_spec ///< TODOCUMENT
                                                                                  ) const {
	alignment_outputter_list alignment_outputters;

	if ( ! get_aln_to_cath_aln_file().empty() ) {
		alignment_outputters.push_back( file_alignment_outputter( get_aln_to_cath_aln_file(), cath_aln_ostream_alignment_outputter() ) );
	}
	if ( get_aln_to_cath_aln_stdout() ) {
		alignment_outputters.push_back( cath_aln_ostream_alignment_outputter() );
	}

	if ( ! get_aln_to_fasta_file().empty() ) {
		alignment_outputters.push_back( file_alignment_outputter( get_aln_to_fasta_file(), fasta_ostream_alignment_outputter() ) );
	}
	if ( get_aln_to_fasta_stdout() ) {
		alignment_outputters.push_back( fasta_ostream_alignment_outputter() );
	}

	if ( ! get_aln_to_ssap_file().empty() ) {
		alignment_outputters.push_back( file_alignment_outputter( get_aln_to_ssap_file(), ssap_ostream_alignment_outputter() ) );
	}
	if ( get_aln_to_ssap_stdout() ) {
		alignment_outputters.push_back( ssap_ostream_alignment_outputter() );
	}

	if ( ! get_aln_to_html_file().empty() ) {
		alignment_outputters.push_back( file_alignment_outputter( get_aln_to_html_file(), make_html_ostream_alignment_outputter( prm_display_spec ) ) );
	}
	if ( get_aln_to_html_stdout() ) {
		alignment_outputters.push_back( make_html_ostream_alignment_outputter( prm_display_spec ) );
	}

	return alignment_outputters;
}

/// \brief TODOCUMENT
bool alignment_output_options_block::outputs_to_stdout() const {
	return ( get_aln_to_fasta_stdout() || get_aln_to_ssap_stdout() || get_aln_to_html_stdout() );
}
