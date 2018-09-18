/// \file
/// \brief The alignment_io class definitions

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

#include "align_scaffold.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/optional.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "alignment/alignment.hpp"
#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/range/indices.hpp"
#include "common/boost_addenda/string_algorithm/split_build.hpp"
#include "common/type_aliases.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace std;

using boost::algorithm::is_any_of;
using boost::algorithm::is_space;
using boost::join;
using boost::make_optional;
using boost::none;

/// \brief Generate data for an alignment entry from a scaffold string
///
/// The position counts up from 0 for all non-gap characters
///
/// Characters are treated is gaps iff they are:
///  * space characters (according to boost::algorithm::is_space()
///  * '-'
///  * '.'
aln_posn_opt_vec cath::align::detail::alignment_entry_of_scaffold_string(const string &prm_scaffold_string ///< The scaffold string defining the entry to be built
                                                                         ) {
	aln_posn_type ctr = 0;
	return transform_build<aln_posn_opt_vec>(
		prm_scaffold_string,
		[&] (const char &x) {
			return ( ! is_space()( x ) && x != '-' && x != '.' ) ? make_optional( ctr++ )
			                                                     : none;
		}
	);
}

/// \brief Generate a scaffold line for the specified entry of the specified alignment
string cath::align::detail::scaffold_line_of_alignment_entry(const alignment &prm_alignment, ///< The alignment containing the entry for which the scaffold line should be made
                                                             const size_t    &prm_entry      ///< The index of the alignment entry for which the scaffold line should be made
                                                             ) {
	return transform_build<string>(
		indices( prm_alignment.length() ),
		[&] (const size_t &x) {
			return has_position_of_entry_of_index( prm_alignment, prm_entry, x ) ? 'X' : ' ';
		}
	);
}

/// \brief Generate an alignment from the specified scaffold_lines
///
/// The position for each entry counts up from 0 for all non-gap characters
///
/// Characters are treated is gaps iff they are:
///  * space characters (according to boost::algorithm::is_space()
///  * '-'
///  * '.'
///
/// \relates alignment
alignment cath::align::alignment_of_scaffold_lines(const str_vec &prm_scaffold_lines ///< The scaffold lines from which to generate the alignment
                                                   ) {
	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return alignment{
		transform_build<aln_posn_opt_vec_vec>(
			prm_scaffold_lines,
			[] (const string &x) {
				return detail::alignment_entry_of_scaffold_string( x );
			}
		)
	};
}

/// \brief Generate an alignment from the specified scaffold string of newline-separated scaffold lines
///
/// The position for each entry counts up from 0 for all non-gap characters
///
/// Characters are treated is gaps iff they are:
///  * space characters (according to boost::algorithm::is_space()
///  * '-'
///  * '.'
///
/// \relates alignment
alignment cath::align::alignment_of_scaffold(const string &prm_scaffold ///< The string of newline-separated scaffold lines from which to generate the alignment
                                             ) {
	return alignment_of_scaffold_lines(
		split_build<str_vec>( prm_scaffold, is_any_of( "\n" ) )
	);
}

/// \brief Generate scaffold lines from the specified alignment
///
/// \returns A vector of strings, one for each entry in the alignment,
///          with each character representing the presence ('X') or absence (' ')
///          of the corresponding position.
///
/// \relates alignment
str_vec cath::align::scaffold_lines_of_alignment(const alignment &prm_alignment ///< The alignment from which the scaffold lines should be generated
                                                 ) {
	return transform_build<str_vec>(
		indices( prm_alignment.num_entries() ),
		[&] (const size_t &entry) {
			return detail::scaffold_line_of_alignment_entry( prm_alignment, entry );
		}
	);
}

/// \brief Generate a newline-separated scaffold of the alignment
///
/// \returns A newline-separated bunch of strings, one for each entry in the alignment,
///          with each character representing the presence ('X') or absence (' ')
///          of the corresponding position.
///
/// \relates alignment
string cath::align::scaffold_of_alignment(const alignment &prm_alignment ///< The alignment from which the scaffold lines should be generated
                                          ) {
	return join( scaffold_lines_of_alignment( prm_alignment ), "\n" );
}
