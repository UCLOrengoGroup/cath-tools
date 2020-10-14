/// \file
/// \brief The filter_vs_full_score_less class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_STRUCTURE_VIEW_CACHE_FILTER_DETAIL_FILTER_VS_FULL_SCORE_LESS_HPP
#define _CATH_TOOLS_SOURCE_UNI_STRUCTURE_VIEW_CACHE_FILTER_DETAIL_FILTER_VS_FULL_SCORE_LESS_HPP

namespace cath { namespace index { namespace filter { class filter_vs_full_score; } } }

namespace cath {
	namespace index {
		namespace filter {
			namespace detail {

				/// \brief A functor for sorting filter_vs_full_scores by filter_score
				class filter_score_less {
				public:
					bool operator()(const filter_vs_full_score &,
					                const filter_vs_full_score &) const;

					bool operator()(const filter_vs_full_score &,
					                const double &) const;
				};

				/// \brief A functor for sorting filter_vs_full_scores by full_score
				class full_score_less {
				public:
					bool operator()(const filter_vs_full_score &,
					                const filter_vs_full_score &) const;

					bool operator()(const filter_vs_full_score &,
					                const double &) const;
				};

			} // namespace detail
		} // namespace filter
	} // namespace index
} // namespace cath

#endif
