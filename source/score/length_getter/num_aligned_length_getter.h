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

#ifndef NUM_ALIGNED_LENGTH_GETTER_H_INCLUDED
#define NUM_ALIGNED_LENGTH_GETTER_H_INCLUDED

#include "score/aligned_pair_score/detail/score_common_coord_handler.h"
#include "score/length_getter/length_getter.h"

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

			virtual std::unique_ptr<length_getter> do_clone() const override final;

			virtual boost::logic::tribool do_higher_is_better() const override final;

			virtual size_t do_get_length(const align::alignment &,
			                             const protein &,
			                             const protein &) const override final;

			virtual length_getter_category do_get_length_getter_category() const override final;

			virtual std::string do_id_name() const override final;

			virtual str_bool_pair_vec do_short_name_suffixes() const override final;

			virtual std::string do_long_name() const override final;

			virtual std::string do_description() const override final;

//			virtual std::string do_short_suffix_string() const;

//			virtual std::string do_long_suffix_string() const;

			virtual const std::string do_description_brackets_string() const override final;

			virtual bool do_less_than_with_same_dynamic_type(const length_getter &) const override final;

			const detail::score_common_coord_handler & get_common_coord_handler() const;

		public:
			num_aligned_length_getter() = default;
			num_aligned_length_getter(const align::common_residue_selection_policy &);
			virtual ~num_aligned_length_getter() noexcept = default;
		};

		bool operator<(const num_aligned_length_getter &,
		               const num_aligned_length_getter &);
	}
}

#endif
