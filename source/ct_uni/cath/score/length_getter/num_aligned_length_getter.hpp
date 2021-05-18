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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_LENGTH_GETTER_NUM_ALIGNED_LENGTH_GETTER_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_LENGTH_GETTER_NUM_ALIGNED_LENGTH_GETTER_HPP

#include "cath/score/aligned_pair_score/detail/score_common_coord_handler.hpp"
#include "cath/score/length_getter/length_getter.hpp"

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

			[[nodiscard]] std::unique_ptr<length_getter> do_clone() const final;

			[[nodiscard]] boost::logic::tribool do_higher_is_better() const final;

			[[nodiscard]] size_t do_get_length( const align::alignment &, const protein &, const protein & ) const final;

			[[nodiscard]] length_getter_category do_get_length_getter_category() const final;

			[[nodiscard]] std::string do_id_name() const final;

			[[nodiscard]] str_bool_pair_vec do_short_name_suffixes() const final;

			[[nodiscard]] std::string do_long_name() const final;

			[[nodiscard]] std::string do_description() const final;

			//			virtual std::string do_short_suffix_string() const;

			//			virtual std::string do_long_suffix_string() const;

			[[nodiscard]] std::string do_description_brackets_string() const final;

			[[nodiscard]] bool do_less_than_with_same_dynamic_type( const length_getter & ) const final;

			[[nodiscard]] const detail::score_common_coord_handler &get_common_coord_handler() const;

		  public:
			num_aligned_length_getter() = default;
			explicit num_aligned_length_getter(const align::common_residue_selection_policy &);
		};

		bool operator<(const num_aligned_length_getter &,
		               const num_aligned_length_getter &);
	} // namespace score
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_LENGTH_GETTER_NUM_ALIGNED_LENGTH_GETTER_HPP
