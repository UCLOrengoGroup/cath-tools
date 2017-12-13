/// \file
/// \brief The num_aligned_length_getter class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_SCORE_LENGTH_GETTER_NUM_ALIGNED_LENGTH_GETTER_HPP
#define _CATH_TOOLS_SOURCE_UNI_SCORE_LENGTH_GETTER_NUM_ALIGNED_LENGTH_GETTER_HPP

#include "score/aligned_pair_score/detail/score_common_coord_handler.hpp"
#include "score/length_getter/length_getter.hpp"

namespace cath {
	namespace score {
		class num_aligned_length_getter;

		bool operator<(const num_aligned_length_getter &,
		               const num_aligned_length_getter &);


		/// \brief TODOCUMENT
		class num_aligned_length_getter final : public length_getter {
		private:
			friend bool operator<(const num_aligned_length_getter &,
			                      const num_aligned_length_getter &);

			/// \brief TODOCUMENT
			detail::score_common_coord_handler common_coord_handler;

			std::unique_ptr<length_getter> do_clone() const final;

			boost::logic::tribool do_higher_is_better() const final;

			size_t do_get_length(const align::alignment &,
			                     const protein &,
			                     const protein &) const final;

			length_getter_category do_get_length_getter_category() const final;

			std::string do_id_name() const final;

			str_bool_pair_vec do_short_name_suffixes() const final;

			std::string do_long_name() const final;

			std::string do_description() const final;

//			virtual std::string do_short_suffix_string() const;

//			virtual std::string do_long_suffix_string() const;

			const std::string do_description_brackets_string() const final;

			bool do_less_than_with_same_dynamic_type(const length_getter &) const final;

			const detail::score_common_coord_handler & get_common_coord_handler() const;

		public:
			num_aligned_length_getter() = default;
			explicit num_aligned_length_getter(const align::common_residue_selection_policy &);
		};

		bool operator<(const num_aligned_length_getter &,
		               const num_aligned_length_getter &);
	} // namespace score
} // namespace cath

#endif
