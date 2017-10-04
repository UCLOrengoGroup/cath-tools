/// \file
/// \brief The true_false_pos_neg class definitions

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

#include "true_false_pos_neg.hpp"

#include "exception/invalid_argument_exception.hpp"
#include "score/true_pos_false_neg/true_false_pos_neg.hpp"

#include <iomanip>
#include <sstream>

using namespace cath;
using namespace cath::common;
using namespace cath::score;
using namespace std;

/// \brief Ctor for true_false_pos_neg from the numbers of true/false positives/negatives
true_false_pos_neg::true_false_pos_neg(const size_t &arg_num_true_positives,  ///< The number of true  positives
                                       const size_t &arg_num_true_negatives,  ///< The number of true  negatives
                                       const size_t &arg_num_false_positives, ///< The number of false positives
                                       const size_t &arg_num_false_negatives  ///< The number of false negatives
                                       ) : num_true_positives ( arg_num_true_positives  ),
                                           num_true_negatives ( arg_num_true_negatives  ),
                                           num_false_positives( arg_num_false_positives ),
                                           num_false_negatives( arg_num_false_negatives ) {

}

/// \brief Incrementer for the specified classn_outcome
true_false_pos_neg & true_false_pos_neg::operator+=(const classn_outcome &arg_outcome ///< The outcome for which the number should be incremented
                                                    ) {
	switch ( arg_outcome ) {
		case ( classn_outcome::TRUE_POSITIVE  ) : { increment_num< classn_outcome::TRUE_POSITIVE  >( *this ); break; }
		case ( classn_outcome::TRUE_NEGATIVE  ) : { increment_num< classn_outcome::TRUE_NEGATIVE  >( *this ); break; }
		case ( classn_outcome::FALSE_POSITIVE ) : { increment_num< classn_outcome::FALSE_POSITIVE >( *this ); break; }
		case ( classn_outcome::FALSE_NEGATIVE ) : { increment_num< classn_outcome::FALSE_NEGATIVE >( *this ); break; }
	}
	return *this;
}

/// \brief Getter for the number of true positives in a true_false_pos_neg
///
/// \relates true_false_pos_neg
const size_t & cath::score::get_num_true_positives(const true_false_pos_neg &arg_true_false_pos_neg ///< The true_false_pos_neg to query
                                                   ) {
	return arg_true_false_pos_neg.get_num<classn_outcome::TRUE_POSITIVE>();
}

/// \brief Getter for the number of true negatives in a true_false_pos_neg
///
/// \relates true_false_pos_neg
const size_t & cath::score::get_num_true_negatives(const true_false_pos_neg &arg_true_false_pos_neg ///< The true_false_pos_neg to query
                                                   ) {
	return arg_true_false_pos_neg.get_num<classn_outcome::TRUE_NEGATIVE>();
}

/// \brief Getter for the number of false positives in a true_false_pos_neg
///
/// \relates true_false_pos_neg
const size_t & cath::score::get_num_false_positives(const true_false_pos_neg &arg_true_false_pos_neg ///< The true_false_pos_neg to query
                                                    ) {
	return arg_true_false_pos_neg.get_num<classn_outcome::FALSE_POSITIVE>();
}

/// \brief Getter for the number of false negatives in a true_false_pos_neg
///
/// \relates true_false_pos_neg
const size_t & cath::score::get_num_false_negatives(const true_false_pos_neg &arg_true_false_pos_neg ///< The true_false_pos_neg to query
                                                    ) {
	return arg_true_false_pos_neg.get_num<classn_outcome::FALSE_NEGATIVE>();
}

/// \brief Setter for the number of true positives in the specified true_false_pos_neg
///
/// \relates true_false_pos_neg
void cath::score::set_num_true_positives(true_false_pos_neg &arg_true_false_pos_neg, ///< The true_false_pos_neg to update
                                         const size_t       &arg_number              ///< The value to which to the number of true positives should be set
                                         ) {
	arg_true_false_pos_neg.set_num<classn_outcome::TRUE_POSITIVE>( arg_number );
}

/// \brief Setter for the number of true negatives in the specified true_false_pos_neg
///
/// \relates true_false_pos_neg
void cath::score::set_num_true_negatives(true_false_pos_neg &arg_true_false_pos_neg, ///< The true_false_pos_neg to update
                                         const size_t       &arg_number              ///< The value to which to the number of true negatives should be set
                                         ) {
	arg_true_false_pos_neg.set_num<classn_outcome::TRUE_NEGATIVE>( arg_number );
}

/// \brief Setter for the number of false positives in the specified true_false_pos_neg
///
/// \relates true_false_pos_neg
void cath::score::set_num_false_positives(true_false_pos_neg &arg_true_false_pos_neg, ///< The true_false_pos_neg to update
                                          const size_t       &arg_number              ///< The value to which to the number of false positives should be set
                                          ) {
	arg_true_false_pos_neg.set_num<classn_outcome::FALSE_POSITIVE>( arg_number );
}

/// \brief Setter for the number of false negatives in the specified true_false_pos_neg
///
/// \relates true_false_pos_neg
void cath::score::set_num_false_negatives(true_false_pos_neg &arg_true_false_pos_neg, ///< The true_false_pos_neg to update
                                          const size_t       &arg_number              ///< The value to which to the number of false negatives should be set
                                          ) {
	arg_true_false_pos_neg.set_num<classn_outcome::FALSE_NEGATIVE>( arg_number );
}



/// \brief Incrementer for the number of true positives in the specified true_false_pos_neg
///
/// \relates true_false_pos_neg
void cath::score::increment_num_true_positives(true_false_pos_neg &arg_true_false_pos_neg ///< The true_false_pos_neg to update
                                               ) {
	return increment_num<classn_outcome::TRUE_POSITIVE>( arg_true_false_pos_neg );
}

/// \brief Incrementer for the number of true negatives in the specified true_false_pos_neg
///
/// \relates true_false_pos_neg
void cath::score::increment_num_true_negatives(true_false_pos_neg &arg_true_false_pos_neg ///< The true_false_pos_neg to update
                                               ) {
	return increment_num<classn_outcome::TRUE_NEGATIVE>( arg_true_false_pos_neg );
}

/// \brief Incrementer for the number of false positives in the specified true_false_pos_neg
///
/// \relates true_false_pos_neg
void cath::score::increment_num_false_positives(true_false_pos_neg &arg_true_false_pos_neg ///< The true_false_pos_neg to update
                                                ) {
	return increment_num<classn_outcome::FALSE_POSITIVE>( arg_true_false_pos_neg );
}

/// \brief Incrementer for the number of false positives in the specified true_false_pos_neg
///
/// \relates true_false_pos_neg
void cath::score::increment_num_false_negatives(true_false_pos_neg &arg_true_false_pos_neg ///< The true_false_pos_neg to update
                                                ) {
	return increment_num<classn_outcome::FALSE_NEGATIVE>( arg_true_false_pos_neg );
}

/// \brief Update a true_false_pos_neg to reflect a switch in some predictions from negative to positive
///
/// In other words, this updates the numbers to reflect moving:
///  * arg_num_positives from false_negatives to true_positives
///  * arg_num_negatives from true_negatives  to false_positives
///
/// \relates true_false_pos_neg
void cath::score::update_with_predicted_positives(true_false_pos_neg &arg_true_false_pos_neg, ///< The true_false_pos_neg to update
                                                  const size_t       &arg_num_positives,      ///< The number of positive cases that were previously predicted as negatives but are now predicted as positives
                                                  const size_t       &arg_num_negatives       ///< The number of negative cases that were previously predicted as negatives but are now predicted as positives
                                                  ) {
	const size_t &num_true_positives  = arg_true_false_pos_neg.get_num< classn_outcome::TRUE_POSITIVE  >();
	const size_t &num_true_negatives  = arg_true_false_pos_neg.get_num< classn_outcome::TRUE_NEGATIVE  >();
	const size_t &num_false_positives = arg_true_false_pos_neg.get_num< classn_outcome::FALSE_POSITIVE >();
	const size_t &num_false_negatives = arg_true_false_pos_neg.get_num< classn_outcome::FALSE_NEGATIVE >();

	if ( arg_num_positives > num_false_negatives ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot move more positive cases out of being predicted negative than there are false_negatives"));
	}
	if ( arg_num_negatives > num_true_negatives  ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot move more negative cases out of being predicted negative than there are true_negatives"));
	}

	arg_true_false_pos_neg.set_num< classn_outcome::TRUE_POSITIVE  >( num_true_positives  + arg_num_positives );
	arg_true_false_pos_neg.set_num< classn_outcome::TRUE_NEGATIVE  >( num_true_negatives  - arg_num_negatives );
	arg_true_false_pos_neg.set_num< classn_outcome::FALSE_POSITIVE >( num_false_positives + arg_num_negatives );
	arg_true_false_pos_neg.set_num< classn_outcome::FALSE_NEGATIVE >( num_false_negatives - arg_num_positives );
}

/// \brief Simple insertion operator for true_false_pos_neg
///
/// \relates true_false_pos_neg
ostream & cath::score::operator<<(ostream                  &arg_os,                ///< The ostream to which the true_false_pos_neg should be output
                                  const true_false_pos_neg &arg_true_false_pos_neg ///< The true_false_pos_neg to output
                                  ) {
	ostringstream temp_ss;
	temp_ss << "true_false_pos_neg[TP:";
	temp_ss << right << setw( 4 ) << arg_true_false_pos_neg.get_num<classn_outcome::TRUE_POSITIVE>();
	temp_ss << "; TN: ";
	temp_ss << right << setw( 4 ) << arg_true_false_pos_neg.get_num<classn_outcome::TRUE_NEGATIVE>();
	temp_ss << "; FP: ";
	temp_ss << right << setw( 4 ) << arg_true_false_pos_neg.get_num<classn_outcome::FALSE_POSITIVE>();
	temp_ss << "; FN: ";
	temp_ss << right << setw( 4 ) << arg_true_false_pos_neg.get_num<classn_outcome::FALSE_NEGATIVE>();
	temp_ss << "]";
	arg_os << temp_ss.str();
	return arg_os;
}

