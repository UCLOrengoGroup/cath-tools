/// \file
/// \brief The ssap_scores class definitions

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

#include "ssap_scores.h"

#include <ostream>

using namespace cath;
using namespace std;

/// \brief Setter for percentage_aligned_pairs_over_larger
void ssap_scores::set_percentage_aligned_pairs_over_larger(const double &arg_percentage_aligned_pairs_over_larger ///< The value to which percentage_aligned_pairs_over_larger should be set
                                                           ) {
	percentage_aligned_pairs_over_larger = arg_percentage_aligned_pairs_over_larger;
}

/// \brief Setter for seq_id
void ssap_scores::set_seq_id(const double &arg_seq_id ///< The value to which seq_id should be set
                             ) {
	seq_id = arg_seq_id;
}

/// \brief Setter for num_aligned_pairs
void ssap_scores::set_num_aligned_pairs(const size_t &arg_num_aligned_pairs ///< The value to which num_aligned_pairs should be set
                                        ) {
	num_aligned_pairs = arg_num_aligned_pairs;
}

/// \brief Setter for ssap_score_over_compared
void ssap_scores::set_ssap_score_over_compared(const double &arg_ssap_score_over_compared ///< The value to which ssap_score_over_compared should be set
                                               ) {
	ssap_score_over_compared = arg_ssap_score_over_compared;
}

/// \brief Setter for ssap_score_over_smaller
void ssap_scores::set_ssap_score_over_smaller(const double &arg_ssap_score_over_smaller ///< The value to which ssap_score_over_smaller should be set
                                              ) {
	ssap_score_over_smaller = arg_ssap_score_over_smaller;
}

/// \brief Setter for ssap_score_over_larger
void ssap_scores::set_ssap_score_over_larger(const double &arg_ssap_score_over_larger ///< The value to which ssap_score_over_larger should be set
                                             ) {
	ssap_score_over_larger = arg_ssap_score_over_larger;
}

/// \brief Getter for percentage_aligned_pairs_over_larger
double ssap_scores::get_percentage_aligned_pairs_over_larger() const {
	return percentage_aligned_pairs_over_larger;
}

/// \brief Getter for seq_id
double ssap_scores::get_seq_id() const {
	return seq_id;
}

/// \brief Getter for num_aligned_pairs
size_t ssap_scores::get_num_aligned_pairs() const {
	return num_aligned_pairs;
}

/// \brief Getter for ssap_score_over_compared
double ssap_scores::get_ssap_score_over_compared() const {
	return ssap_score_over_compared;
}

/// \brief Getter for ssap_score_over_smaller
double ssap_scores::get_ssap_score_over_smaller() const {
	return ssap_score_over_smaller;
}

/// \brief Getter for ssap_score_over_larger
double ssap_scores::get_ssap_score_over_larger() const {
	return ssap_score_over_larger;
}

/// \brief Basic insertion operator to output a rough summary of an ssap_scores object to an ostream
///
/// \relates ssap_scores
ostream & cath::operator<<(ostream           &arg_os,         ///< The ostream to which the ssap_scores should be output
                           const ssap_scores &arg_ssap_scores ///< The ssap_scores to output
                           ) {
	arg_os << "ssap_scores[";
	arg_os << "\tget_num_aligned_pairs                    : " << arg_ssap_scores.get_num_aligned_pairs()                    << "\n";
	arg_os << "\tget_percentage_aligned_pairs_over_larger : " << arg_ssap_scores.get_percentage_aligned_pairs_over_larger() << "\n";
	arg_os << "\tget_seq_id                               : " << arg_ssap_scores.get_seq_id()                               << "\n";
	arg_os << "\tget_ssap_score_over_compared             : " << arg_ssap_scores.get_ssap_score_over_compared()             << "\n";
	arg_os << "\tget_ssap_score_over_smaller              : " << arg_ssap_scores.get_ssap_score_over_smaller()              << "\n";
	arg_os << "\tget_ssap_score_over_larger               : " << arg_ssap_scores.get_ssap_score_over_larger()               << "\n";
	arg_os << "]";
	return arg_os;
}

