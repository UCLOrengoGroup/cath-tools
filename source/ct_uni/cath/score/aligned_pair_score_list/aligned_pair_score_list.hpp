/// \file
/// \brief The aligned_pair_score_list class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST_ALIGNED_PAIR_SCORE_LIST_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST_ALIGNED_PAIR_SCORE_LIST_HPP

#include <boost/ptr_container/ptr_vector.hpp>

#include "cath/common/type_aliases.hpp"

namespace cath {
	namespace score {
		class aligned_pair_score;

		/// \brief Store a list of aligned_pair_scores
		class aligned_pair_score_list final {
		private:
			/// \brief TODOCUMENT
			using aligned_pair_score_ptr_vec = boost::ptr_vector<aligned_pair_score>;

			/// \brief The list of the scores
			aligned_pair_score_ptr_vec scores;

		public:
			/// \brief TODOCUMENT
			using const_iterator = aligned_pair_score_ptr_vec::const_iterator;

			/// @{
			aligned_pair_score_list() = default;
			explicit aligned_pair_score_list(const boost::ptr_vector<aligned_pair_score> &);
			/// @}

			void add_aligned_pair_score(const aligned_pair_score &);

			/// @{
			[[nodiscard]] size_t       size() const;
			aligned_pair_score & operator[](const size_t &);
			const aligned_pair_score & operator[](const size_t &) const;
			/// @}

			/// @{
			[[nodiscard]] const_iterator begin() const;
			[[nodiscard]] const_iterator end() const;
			/// @}
		};

		void warn_on_duplicate_human_friendly_names(const aligned_pair_score_list &);

		str_vec get_short_names(const aligned_pair_score_list &);
		str_vec get_long_names(const aligned_pair_score_list &);
		str_vec get_descriptions(const aligned_pair_score_list &);
		str_vec get_references(const aligned_pair_score_list &);
	} // namespace score
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LIST_ALIGNED_PAIR_SCORE_LIST_HPP
