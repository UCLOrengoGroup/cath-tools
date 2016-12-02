/// \file
/// \brief The horiz_align_outputter class definitions

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

#include "horiz_align_outputter.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "alignment/alignment.hpp"

#include <iomanip>

using namespace cath::align;
using namespace std;

using boost::lexical_cast;
using boost::numeric_cast;

/// \brief Ctor for horiz_align_outputter
horiz_align_outputter::horiz_align_outputter(const alignment &arg_alignment ///< The alignment to be output
                                             ) : the_alignment( arg_alignment ) {
}

/// \brief Getter for the const reference to the alignment
const alignment & horiz_align_outputter::get_alignment() const {
	return the_alignment;
}

/// \brief Output the alignment to the ostream in horizontal format
///
/// This outputs the alignment itself, ie it outputs position numbers rather the things to which those numbers refer.
/// This means it can be used on a bare alignment without the need for related data.
///
/// The pads all the numbers with spaces so that they all have the same width as the maximum
/// position number.
///
/// \relates horiz_align_outputter
ostream & cath::align::operator<<(ostream                     &arg_os,                   ///< The ostream to which the alignment should be output
                                  const horiz_align_outputter &arg_horiz_align_outputter ///< A horiz_align_outputter that wraps the alignment to be output
                                  ) {
	// Grab alignment and then some basic information from it
	const alignment            the_alignment   = arg_horiz_align_outputter.get_alignment();
	const alignment::size_type length          = the_alignment.length();
	const alignment::size_type num_entries     = the_alignment.num_entries();
	const aln_posn_opt         max_index_opt   = get_max_last_present_position( the_alignment );
	const aln_posn_type        max_index       = max_index_opt ? max_index_opt.get() : 0;
	const size_t               max_index_width = lexical_cast<string>( max_index ).length();

	// Loop over the positions, and output them
	arg_os << "alignment[";
	for (alignment::size_type entry_ctr = 0; entry_ctr < num_entries; ++entry_ctr) {
		arg_os << "\n\t";
		for (alignment::size_type index_ctr = 0; index_ctr < length; ++index_ctr) {
			arg_os << " ";
			const aln_posn_opt position = the_alignment.position_of_entry_of_index( entry_ctr, index_ctr );
			if ( position ) {
				ostringstream posn_ss;
				posn_ss << setw( numeric_cast<int>( max_index_width ) ) << setfill( ' ' );
				posn_ss << *position;
				arg_os  << posn_ss.str();
			}
			else {
				arg_os << string( max_index_width - 1, ' ') << "-";
			}
		}
	}
	arg_os << ( ( num_entries > 0 ) ? "\n" : "" );
	arg_os << "]";

	// Return the specified ostream
	return arg_os;
}
