/// \file
/// \brief The filter_vs_full_score class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_FILTER_FILTER_VS_FULL_SCORE_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_FILTER_FILTER_VS_FULL_SCORE_HPP

#include <boost/operators.hpp>

#include "cath/score/true_pos_false_neg/classn_outcome.hpp"

#include <iosfwd>

namespace cath::index::filter {

	/// \brief A score produced by a filter along with the full score produced by complete analysis
	///
	/// Given a score that's slow to calculate for entries (eg SSAP score for previously unaligned
	/// pairs of structures), it can be helpful to use a fast filter to generate a score that
	/// identifies entries likely to get good scores. That are worth doing the slow calculations on.
	///
	/// This class is useful in analysing how effectively a given filter's score
	/// find good entries whilst rejecting bad entries.
	///
	/// In addition to representing an actual filter/full score pair for some entry, it can also
	/// represent an attempted classification: attempting to find all entries with full scores >= some
	/// value by choosing all entries with filter scores >= some value
	class filter_vs_full_score final : private boost::equality_comparable<filter_vs_full_score> {
	private:
		/// \brief The score generated by the filter for a particular entry
		double filter_score = 0.0;

		/// \brief The full score that would be generated
		double full_score   = 0.0;

	public:
		filter_vs_full_score() = default;
		filter_vs_full_score(const double &,
		                     const double &);

		[[nodiscard]] const double &get_filter_score() const;
		[[nodiscard]] const double &get_full_score() const;
	};

	bool operator==(const filter_vs_full_score &,
	                const filter_vs_full_score &);

	score::classn_outcome assess_real_scores_on_filter_attempt(const filter_vs_full_score &,
	                                                           const filter_vs_full_score &);

	std::istream & operator>>(std::istream &,
	                          filter_vs_full_score &);

	std::ostream & operator<<(std::ostream &,
	                          const filter_vs_full_score &);

} // namespace cath::index::filter

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_VIEW_CACHE_FILTER_FILTER_VS_FULL_SCORE_HPP
