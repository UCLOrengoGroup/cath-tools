/// \file
/// \brief The gen_dyn_prog_string_aligner class definitions

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

#include "gen_dyn_prog_string_aligner.hpp"

#include <boost/lexical_cast.hpp>

#include "cath/alignment/align_type_aliases.hpp"
#include "cath/alignment/alignment.hpp"
#include "cath/alignment/dyn_prog_align/dyn_prog_aligner.hpp"
#include "cath/alignment/dyn_prog_align/dyn_prog_score_source/sequence_string_dyn_prog_score_source.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"

#include <string>

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::align::detail;
using namespace ::cath::align::gap;
using namespace ::cath::common;
using namespace ::std;

using ::boost::lexical_cast;

/// \brief Check that there is dyn_prog_aligner stored in dyn_prog_aligner_ptr
///
/// \pre There must be a dyn_prog_aligner stored in dyn_prog_aligner_ptr, else this will throw an out_of_range_exception
void gen_dyn_prog_string_aligner::check_dyn_prog_aligner_ptr() const {
	if ( dyn_prog_aligner_ptr == nullptr ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception("Attempt to use gen_dyn_prog_string_aligner with unpopulated dyn_prog_aligner_ptr"));
	}
}

/// \brief TODOCUMENT
dyn_prog_aligner & gen_dyn_prog_string_aligner::get_dyn_prog_aligner() {
	return *dyn_prog_aligner_ptr;
}

/// \brief TODOCUMENT
const dyn_prog_aligner & gen_dyn_prog_string_aligner::get_dyn_prog_aligner() const {
	return *dyn_prog_aligner_ptr;
}

/// \brief Concrete definition of do_align() that uses the held dyn_prog_aligner to perform a generic alignment on the two sequences
str_str_pair gen_dyn_prog_string_aligner::do_align(const string      &prm_string_a,   ///< The first  string to be aligned
                                                   const string      &prm_string_b,   ///< The second string to be aligned
                                                   const gap_penalty &prm_gap_penalty ///< The gap penalty to be applied for each gap step (ie for opening OR extending a gap)
                                                   ) const {
	// Check that there is dyn_prog_aligner stored in dyn_prog_aligner_ptr
	check_dyn_prog_aligner_ptr();

	// Align the pair of strings (with the specified gap penalty) using the dyn_prog_aligner_ptr
	const string::size_type max_length = max( prm_string_a.length(), prm_string_b.length() );
	const score_alignment_pair score_and_alignment = get_dyn_prog_aligner().align(
		sequence_string_dyn_prog_score_source(
			prm_string_a,
			prm_string_b
		),
		prm_gap_penalty,
		max_length
	);

	// Grab a reference to the alignment and use it to form two alignment strings and return the result
	const alignment    &the_alignment   = score_and_alignment.second;
	const str_str_pair  aligned_strings = format_alignment_strings(
		the_alignment,
		prm_string_a,
		prm_string_b
	);

	// If running in debug mode, check that the score is the same as a freshly calculated one
#ifndef NDEBUG
	const score_type &score              = score_and_alignment.first;
	const score_type  recalculated_score = get_score_of_aligned_sequence_strings(
		aligned_strings.first,
		aligned_strings.second,
		prm_gap_penalty
	);
	if ( score != recalculated_score ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception(
			"Score retrieved from alignment \""
			+ aligned_strings.first
			+ "\" <-> \""
			+ aligned_strings.second
			+ "\" was "
			+ lexical_cast<string>(score)
			+ " but, based on a recalculation, it should be "
			+ lexical_cast<string>(recalculated_score)
		));
	}
#endif /* #ifndef NDEBUG */

	return aligned_strings;
}

/// \brief Ctor for gen_dyn_prog_string_aligner
gen_dyn_prog_string_aligner::gen_dyn_prog_string_aligner(const dyn_prog_aligner &prm_dyn_prog_aligner ///< The dyn_prog_aligner that should be used to perform the alignments
                                                         ) : dyn_prog_aligner_ptr(prm_dyn_prog_aligner.clone()) {
}

