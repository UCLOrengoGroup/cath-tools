/// \file
/// \brief The filter_vs_full_score_less class definitions

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

#include "filter_vs_full_score_less.h"

#include "structure/view_cache/filter/filter_vs_full_score.h"

using namespace cath::index::filter::detail;

/// \brief A less-than operator for filter_vs_full_scores that uses the filter_score
bool filter_score_less::operator()(const filter_vs_full_score &arg_filter_vs_full_score_lhs, ///< The first filter_vs_full_score to compare
                                   const filter_vs_full_score &arg_filter_vs_full_score_rhs  ///< The second filter_vs_full_score to compare
                                   ) const {
	return ( arg_filter_vs_full_score_lhs.get_filter_score() < arg_filter_vs_full_score_rhs.get_filter_score() );
}

/// \brief A less-than operator for filter_vs_full_score that compares the filter_vs_full_score's
///        filter score against an explicit filter score value
///
/// This can be useful for, say, using lower_bound to find the position for a filter value in a
/// filter-value-sorted range of filter_vs_full_scores
bool filter_score_less::operator()(const filter_vs_full_score &arg_filter_vs_full_score_lhs, ///< The first filter_vs_full_score to compare
                                   const double               &arg_filter_score_rhs          ///< The explicit filter value against which to compare
                                   ) const {
	return ( arg_filter_vs_full_score_lhs.get_filter_score() < arg_filter_score_rhs );
}

/// \brief A less-than operator for filter_vs_full_scores that uses the full score
bool full_score_less::operator()(const filter_vs_full_score &arg_filter_vs_full_score_lhs, ///< The first filter_vs_full_score to compare
                                 const filter_vs_full_score &arg_filter_vs_full_score_rhs  ///< The second filter_vs_full_score to compare
                                 ) const {
return ( arg_filter_vs_full_score_lhs.get_full_score() < arg_filter_vs_full_score_rhs.get_full_score() );
}

/// \brief A less-than operator for filter_vs_full_score that compares the filter_vs_full_score's
///        full score against an explicit full score value
///
/// This can be useful for, say, using lower_bound to find the position for a full value in a
/// full-value-sorted range of filter_vs_full_scores
bool full_score_less::operator()(const filter_vs_full_score &arg_filter_vs_full_score_lhs, ///< The first filter_vs_full_score to compare
                                 const double               &arg_full_score_rhs            ///< The explicit full value against which to compare
                                 ) const {
return ( arg_filter_vs_full_score_lhs.get_full_score() < arg_full_score_rhs );
}
