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

#include "alignment_breaks.h"

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/filtered.hpp>

#include <boost/range/irange.hpp>

#include "alignment/alignment.h"
#include "common/algorithm/copy_build.h"
#include "common/algorithm/sets_are_disjoint.h"
#include "common/boost_addenda/range/adaptor/adjacented.h"
#include "common/size_t_literal.h"

using namespace cath;
using namespace cath::align;
using namespace cath::common;

using boost::adaptors::filtered;
using boost::irange;

/// \brief TODOCUMENT
size_vec cath::align::get_alignment_breaks(const alignment &arg_alignment ///< TODOCUMENT
                                           ) {
	// For each alignment index *after* zero, check whether there's any overlap
	// between the entries that are present at the previous index and the entries that are
	// present at this index. If not, add this index to the list.
	//
	// Note: Would prefer to use Boost Range's adjacent_filtered adaptor here but
	// it doesn't do what's required here: it always lets the first pair through,
	// even if they fail the predicate.
	return copy_build<size_vec>(
		irange( 1_z, arg_alignment.length() )
			| filtered (
				[&] (const size_t &curr_idx) {
					// Return whether the set of present entries for the previous index
					// is disjoint with the set of present entries for this index
					return sets_are_disjoint(
						present_positions_of_index( arg_alignment, curr_idx - 1 ),
						present_positions_of_index( arg_alignment, curr_idx     )
					);
				}
			)
	);
}
