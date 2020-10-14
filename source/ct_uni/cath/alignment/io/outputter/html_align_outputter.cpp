/// \file
/// \brief The html_align_outputter class definitions

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

#include "html_align_outputter.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "cath/alignment/alignment.hpp"
#include "cath/alignment/alignment_context.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/common/algorithm/contains.hpp"
#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/display/display_colour_spec/display_colour_spec.hpp"
#include "cath/display/display_colourer/display_colourer.hpp"
#include "cath/display/viewer/viewer.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/file/pdb/pdb_list.hpp"
#include "cath/file/pdb/pdb_residue.hpp"
#include "cath/superposition/superposition_context.hpp"

#include <iomanip>

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::chop;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::sup;
using namespace ::std;

/// \brief Ctor for html_align_outputter
html_align_outputter::html_align_outputter(const alignment          &prm_alignment, ///< The alignment to be output
                                           const strucs_context     &prm_context,   ///< TODOCUMENT
                                           const display_colourer   &prm_colourer   ///< TODOCUMENT
                                           ) : the_alignment { prm_alignment },
                                               context       { prm_context   },
                                               colourer      { prm_colourer  } {
}

/// \brief Getter for the const reference to the alignment
const alignment & html_align_outputter::get_alignment() const {
	return the_alignment;
}

/// \brief Getter for the const reference to the strucs_context
const strucs_context & html_align_outputter::get_strucs_context() const {
	return context;
}

/// \brief TODOCUMENT
const display_colourer & html_align_outputter::get_display_colourer() const {
	return colourer;
}

/// \brief TODOCUMENT
const pdb_list & cath::align::get_pdbs(const html_align_outputter &prm_html_align_outputter ///< TODOCUMENT
                                       ) {
	return prm_html_align_outputter.get_strucs_context().get_pdbs();
}

/// \brief TODOCUMENT
const name_set_list & cath::align::get_name_sets(const html_align_outputter &prm_html_align_outputter ///< TODOCUMENT
                                                 ) {
	return prm_html_align_outputter.get_strucs_context().get_name_sets();
}

/// \brief The specification of the regions of the PDBs to which the alignment refers
const region_vec_opt_vec & cath::align::get_regions(const html_align_outputter &prm_html_align_outputter ///< TODOCUMENT
                                                    ) {
	return prm_html_align_outputter.get_strucs_context().get_regions();
}

/// \brief Make an html_align_outputter from the specified alignment_context and display_colourer
///
/// \relates html_align_outputter
html_align_outputter cath::align::make_html_align_outputter(const alignment_context &prm_alignment_context, ///< The alignment_context from which to build the html_align_outputter
                                                            const display_colourer  &prm_display_colourer   ///< The display_colourer from which to build the html_align_outputter
                                                            ) {
	return {
		prm_alignment_context.get_alignment(),
		prm_alignment_context.get_strucs_context(),
		prm_display_colourer
	};
}

/// \brief Output the alignment to the ostream in horizontal format
///
/// This outputs the alignment itself, ie it outputs position numbers rather the things to which those numbers refer.
/// This means it can be used on a bare alignment without the need for related data.
///
/// The pads all the numbers with spaces so that they all have the same width as the maximum
/// position number.
///
/// \relates html_align_outputter
ostream & cath::align::operator<<(ostream                    &prm_os,                  ///< The ostream to which the alignment should be output
                                  const html_align_outputter &prm_html_align_outputter ///< A html_align_outputter that wraps the alignment to be output
                                  ) {
	// Grab alignment and then some basic information from it
	const alignment           &the_alignment   = prm_html_align_outputter.get_alignment();
	const pdb_list            &pdbs            = get_pdbs     ( prm_html_align_outputter );
	const region_vec_opt_vec  &regions         = get_regions  ( prm_html_align_outputter );
	const name_set_list       &name_sets       = get_name_sets( prm_html_align_outputter );
	const str_vec             &names           = get_alignment_html_names( name_sets );
	const display_colourer    &colourer        = prm_html_align_outputter.get_display_colourer();
	const alignment::size_type length          = the_alignment.length();
	const alignment::size_type num_entries     = the_alignment.num_entries();
	const display_colour_spec  colour_spec     = get_colour_spec( colourer, prm_html_align_outputter.get_strucs_context(), the_alignment );
	const display_colour_vec   colours         = get_all_colours( colour_spec );
	const size_t               num_colours     = colours.size();
	const str_vec              colour_names    = generate_colour_names( num_colours, colour_category::STRUC_OR_RES );

	const size_display_colour_map      &colour_of_pdb_map         = get_clr_of_pdb( colour_spec );
	const size_size_display_colour_map &colour_of_pdb_and_res_map = colour_spec.get_clr_of_pdb_and_res();

	prm_os << "<html>\n";

	if ( false ) {
		prm_os << R"(
<html>

<head>
	<script src="https://code.jquery.com/jquery-3.2.1.min.js"
	        integrity="sha256-hwg4gsxgFZhOsEEamdOYGBf13FyQuiTwlAQgxVSNgt4="
	        crossorigin="anonymous"></script>

	<script>

		$(function() {

			jQuery.fn.order = function(asc, fn) {
				fn = fn || function (el) {
					return $(el).text().replace(/^\s+|\s+$/g, '');
				};
				var T = asc !== false ? 1 : -1,
					F = asc !== false ? -1 : 1;
				this.sort(function (a, b) {
					a = fn(a), b = fn(b);
					if (a == b) return 0;
					return a < b ? F : T;
				});
				this.each(function (i) {
					this.parentNode.appendChild(this);
				});
			};

			$('#aln_container .seq').order(true, function (el) {
				var seq = $(el).find('.seq-name').first().text()
				var seqid = seq.split(/\s+/)[0];
				return seqid;
			});

		});

	</script>

</head>
)";
	}

	prm_os << R"(
<style>
#aln_container {
	font-family : "Liberation Mono", courier, monospace;
	font-size   : 9pt;
	overflow-x  : auto;
	width       : 100%;
}
#aln_container .seq      {
	position    : relative;
	white-space : pre;
	width       : 100%;
}
#aln_container .seq-name {
	float       : left;
	width       : 100px;
}
#aln_container .seq-res  {
	margin-left : 110px;
}
)";

	str_str_map colour_name_of_hex_string;
	for (const size_t &colour_ctr : indices( num_colours ) ) {
		const display_colour &colour      = colours      [ colour_ctr ];
		const string         &colour_name = colour_names [ colour_ctr ];
		const string         &colour_hex  = hex_string_of_colour( colour );
		prm_os << "span." << colour_name << " {background: #" << colour_hex << "}\n";
		colour_name_of_hex_string[ colour_hex ] = colour_name;
	}
	prm_os << R"(span.gap {
	background  : white;
	color       : white;
}
</style>
<body>
<div id="aln_container">
)";

	// Loop over the positions, and output them
	for (const alignment::size_type &entry_ctr : indices( num_entries ) ) {
		const pdb    &the_pdb     = pdbs [ entry_ctr ];
		const string &name        = names[ entry_ctr ];
		const auto   &pdb_regions = regions[ entry_ctr ];
		prm_os << "<div class=\"seq\">";
		prm_os << "<div class=\"seq-name\">&gt;" << name << "</div>";
		prm_os << "<div class=\"seq-res\">";
		for (const alignment::size_type &index_ctr : indices( length ) ) {
			const aln_posn_opt position = the_alignment.position_of_entry_of_index( entry_ctr, index_ctr );
			if ( position ) {
				const char           amino_acid_letter  = get_amino_acid_letter_tolerantly( get_residue_of_region_limited_backbone_complete_index( the_pdb, *position, pdb_regions ) );
				const size_size_pair entry_and_res_pair = make_pair( entry_ctr, *position );
				const bool           has_residue_colour = contains( colour_of_pdb_and_res_map, entry_and_res_pair );
				const display_colour the_colour         = has_residue_colour ? colour_of_pdb_and_res_map.at( entry_and_res_pair )
				                                                             : colour_of_pdb_map.at        ( entry_ctr          );
				const string         the_colour_hex     = hex_string_of_colour( the_colour );
				const string         the_colour_name    = colour_name_of_hex_string.at( the_colour_hex );

				prm_os << "<span class=" << the_colour_name << ">" << amino_acid_letter << "</span>";
			}
			else {
				prm_os << "<span class=gap>-</span>";
			}
		}
		prm_os << "</div></div>\n";
	}
	prm_os << "</div>";
	prm_os << "</body>";
	prm_os << "</html>";

	// Return the specified ostream
	return prm_os;
}
