/// \file
/// \brief The aligned_pair_score class definitions

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

#include "aligned_pair_score.hpp"

#include <boost/logic/tribool.hpp>
#include <boost/range/algorithm/find.hpp>

#include "cath/alignment/alignment.hpp"
#include "cath/alignment/pair_alignment.hpp"
#include "cath/common/algorithm/contains.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/boost_addenda/ptr_container/unique_ptr_functions.hpp"
#include "cath/common/clone/check_uptr_clone_against_this.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"
#include "cath/score/aligned_pair_score_list/aligned_pair_score_list.hpp"
#include "cath/score/aligned_pair_score_list/aligned_pair_score_list_factory.hpp"
#include "cath/score/detail/score_name_helper.hpp"
#include "cath/structure/protein/protein.hpp"
#include "cath/structure/protein/residue.hpp"

using namespace cath;
using namespace cath::align;
using namespace cath::common;
using namespace cath::score;
using namespace cath::score::detail;
using namespace std;

using boost::logic::tribool;

/// \brief Return a map of the different concrete aligned_pair_score types by the human_friendly_short_name string they return
///
/// \relates aligned_pair_score
str_aligned_pair_score_pmap cath::score::detail::get_aligned_pair_score_of_id_name() {
	str_aligned_pair_score_pmap aligned_pair_score_of_id_name;
	for (const aligned_pair_score &the_score : make_one_of_each_type_aligned_pair_score_list() ) {
		insert( aligned_pair_score_of_id_name, the_score.id_name(), the_score.clone() );
	}
	return aligned_pair_score_of_id_name;
}

/// \brief A default implementation for providing an empty string for the long name, which is
///        interpreted as returning the full_short_name() as the long_name() too.
string aligned_pair_score::do_long_name() const {
	return "";
}

/// \brief A default implementation for providing an empty string where it isn't appropriate to
///        provide a reference to a publication for the score.
string aligned_pair_score::do_reference() const {
	return "";
}


/// \brief Standard approach to achieving a virtual copy-ctor
unique_ptr<aligned_pair_score> aligned_pair_score::clone() const {
	return check_uptr_clone_against_this( do_clone(), *this );
}

/// \brief An NVI pass-through to the concrete class's do_higher_is_better() which defines whether a higher score
///        (rather than a lower score) generally reflects a better alignment and/or more similar structures.
tribool aligned_pair_score::higher_is_better() const {
	return do_higher_is_better();
}

/// \brief An NVI pass-through to the concrete class's do_calculate() which defines the method of calculating a score for
///        a pair alignment and two associated proteins.
///
/// \pre Alignment must be a pair alignment (ie prm_alignment.num_entries is alignment::NUM_ENTRIES_IN_PAIR_ALIGNMENT (2) )
///      else an invalid_argument_exception will be thrown.
///
/// \pre The alignment must contain at least one position for each protein
///      else an invalid_argument_exception will be thrown.
///
/// \pre The alignment must not overrun the end of the two associated proteins
///      else an invalid_argument_exception will be thrown.
score_value aligned_pair_score::calculate(const alignment &prm_alignment, ///< The pair alignment to be scored
                                          const protein   &prm_protein_a, ///< The protein associated with the first  half of the alignment
                                          const protein   &prm_protein_b  ///< The protein associated with the second half of the alignment
                                          ) const {
	// Sanity check that the input alignment has two entries
	if ( alignment::NUM_ENTRIES_IN_PAIR_ALIGNMENT != prm_alignment.num_entries() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Alignment passed to aligned_pair_score but be a pair alignment"));
	}

	// Sanity check that the first half of the alignment doesn't overrun prm_protein_a
	// (implicitly checks that the alignment contains at least one position for the first half)
	if ( get_last_present_a_position(prm_alignment) >= prm_protein_a.get_length() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("First half of alignment passed to aligned_pair_score overruns prm_protein_a"));
	}
	// Sanity check that the second half of the alignment doesn't overrun prm_protein_b
	// (implicitly checks that the alignment contains at least one position for the second half)
	if ( get_last_present_b_position(prm_alignment) >= prm_protein_b.get_length() ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Second half of alignment passed to aligned_pair_score overruns prm_protein_b"));
	}

	// Return the result of the concrete score's calculation
	return do_calculate( prm_alignment, prm_protein_a, prm_protein_b );
}

/// \brief An NVI pass-through to the concrete class's do_long_name() which defines a free text description, describing the score.
string aligned_pair_score::description() const {
	return do_description();
}

/// \brief TODOCUMENT
string aligned_pair_score::id_name() const {
	return do_id_name();
}

/// \brief TODOCUMENT
str_bool_pair_vec aligned_pair_score::short_name_suffixes() const {
	return do_short_name_suffixes();
}

/// \brief TODOCUMENT
///
/// \post The returned string will be non-empty and will contain no spaces
///       (else an out_of_range_exception will be thrown)
string aligned_pair_score::human_friendly_short_name() const {
	return score_name_helper::human_friendly_short_name( id_name(), short_name_suffixes() );
}

/// \brief TODOCUMENT
///
/// \post The returned string will be non-empty and will contain no spaces
///       (else an out_of_range_exception will be thrown)
string aligned_pair_score::full_short_name() const {
	return score_name_helper::full_short_name( id_name(), short_name_suffixes() );
}

/// \brief An NVI pass-through to the concrete class's do_long_name() which defines a long name that may expand acronyms and contain spaces.
///        do_long_name() may return an empty string (default), in which case the human_friendly_short_name() will also be returned as the long_name().
string aligned_pair_score::long_name() const {
	const string the_long_name = do_long_name();
	return the_long_name.empty() ? full_short_name() : the_long_name;
}

/// \brief An NVI pass-through to the concrete class's do_reference() which defines a free text description, providing a string
///        describing a publication reference for the score, if appropriate, or returning an empty string otherwise (default).
string aligned_pair_score::reference() const {
	return do_reference();
}

///// \brief An NVI pass-through to the concrete class's do_build_from_short_name_spec() which defines how to build an aligned_pair_score
/////        of the concrete type from a short_name_spec string
//unique_ptr<aligned_pair_score> aligned_pair_score::build_from_short_name_spec(const string &prm_short_name_spec ///< TODOCUMENT
//                                                                              ) const {
//	return do_build_from_short_name_spec( prm_short_name_spec );
//}

/// \brief An NVI pass-through to the concrete class's do_less_than_with_same_dynamic_type(),
///        which defines the less-than operator when the argument's known to have the same dynamic type
bool aligned_pair_score::less_than_with_same_dynamic_type(const aligned_pair_score &prm_aligned_pair_score ///< TODOCUMENT
                                                          ) const {
	assert( typeid( *this ) == typeid( prm_aligned_pair_score ) );
	return do_less_than_with_same_dynamic_type( prm_aligned_pair_score );
}

///// \brief Build an aligned_pair_score from the short name it produces
/////
///// \relates aligned_pair_score
//unique_ptr<aligned_pair_score> cath::score::make_aligned_pair_score_from_full_short_name(const string &prm_short_name ///< The human_friendly_short_name that specifies which aligned_pair_score to build
//                                                                                         ) {
//	// Create a string containing the start of prm_short_name up to the first full-stop
//	const auto first_full_stop_itr = find( prm_short_name, '.' );
//	const auto the_id_name         = string( common::cbegin( prm_short_name ), first_full_stop_itr );
//
//	// Grab all concrete types of aligned_pair_score by their id_name and throw if the_id_name isn't in it
//	const auto aligned_pair_score_of_id_name = detail::get_aligned_pair_score_of_id_name();
//	if ( ! cath::common::contains( aligned_pair_score_of_id_name, the_id_name ) ) {
//		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot find any aligned_pair_score with the specified human_friendly_short_name beginning " + the_id_name.substr( 0, 100 )));
//	}
//
//	// Grab the aligned_pair_score associated with the_id_name
//	const auto &the_score = aligned_pair_score_of_id_name.at( the_id_name );
//
//	// Grab any further characters if a full-stop was found or an empty string otherwise
//	const bool has_more_chars  = ( first_full_stop_itr != common::cend( prm_short_name ) );
//	const auto short_name_spec = has_more_chars ? string( next( first_full_stop_itr ), common::cend( prm_short_name ) )
//	                                            : string();
//
//	// Return the result of getting the concrete type to build from the remainder of the string
//	return the_score.build_from_short_name_spec( short_name_spec );
//}

/// \brief Simple insertion operator for aligned_pair_score
///
/// \relates aligned_pair_score
ostream & cath::score::operator<<(ostream                  &prm_os,                ///< The ostream to which the aligned_pair_score should be output
                                  const aligned_pair_score &prm_aligned_pair_score ///< The aligned_pair_score to output
                                  ) {
	prm_os << ( "aligned_pair_score[" + prm_aligned_pair_score.human_friendly_short_name() + "]" );
	return prm_os;
}
