/// \file
/// \brief The alignment_residue_scores class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_RESIDUE_SCORE_ALIGNMENT_RESIDUE_SCORES_H
#define _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_RESIDUE_SCORE_ALIGNMENT_RESIDUE_SCORES_H

#include "common/type_aliases.hpp"

namespace cath {
	namespace align {
		class alignment;

		/// \brief TODOCUMENT
		///
		/// TODOCUMENT other_present_entries / all_other_entries
		///
		/// TODOCUMENT Doesn't manage presence
		class alignment_residue_scores final {
		private:
			/// \brief TODOCUMENT
			size_t            num_entries;

			/// \brief TODOCUMENT
			size_vec          num_present_entries_by_index;

			/// \brief TODOCUMENT
			score_opt_vec_vec scores_to_other_present_entries;

			void sanity_check() const;

		public:
			alignment_residue_scores(const size_t &,
			                         size_vec,
			                         score_opt_vec_vec);

			size_t get_num_entries() const;
			size_t get_length() const;
			size_t get_num_present_entries_of_index(const size_t &) const;

			score_opt get_opt_score_to_other_present_entries(const size_t &,
			                                                 const size_t &) const;
		};

		bool has_score(const alignment_residue_scores &,
		               const size_t &,
		               const size_t &);

		float_score_type get_unnormalised_score(const alignment_residue_scores &,
		                                        const size_t &,
		                                        const size_t &,
		                                        const bool &);

		float_score_type get_max_score(const alignment_residue_scores &,
		                               const bool &);

		float_score_type get_score(const alignment_residue_scores &,
		                           const size_t &,
		                           const size_t &,
		                           const bool &,
		                           const bool &);

		float_score_type get_score_to_present_entries_of_index(const alignment_residue_scores &,
		                                                       const size_t &,
		                                                       const size_t &);

		float_score_type get_score_to_all_entries(const alignment_residue_scores &,
		                                          const size_t &,
		                                          const size_t &);

		float_score_type get_normalised_score_to_present_entries_of_index(const alignment_residue_scores &,
		                                                                  const size_t &,
		                                                                  const size_t &);

		float_score_type get_normalised_score_to_all_entries(const alignment_residue_scores &,
		                                                     const size_t &,
		                                                     const size_t &);

		alignment_residue_scores make_alignment_residue_scores(const alignment &,
		                                                       const score_opt_vec_vec &);

	} // namespace align
} // namespace cath

#endif
