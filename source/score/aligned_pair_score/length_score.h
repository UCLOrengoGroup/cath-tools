/// \file
/// \brief The length_score class header

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

#ifndef LENGTH_SCORE_H_INCLUDED
#define LENGTH_SCORE_H_INCLUDED

#include "common/clone/clone_ptr.h"
#include "score/aligned_pair_score/aligned_pair_score.h"
#include "score/length_getter/num_aligned_length_getter.h"

namespace cath {
	namespace score {

		/// \brief Calculate (and represent) the number of aligned residues in an alignment.
		///
		/// \todo Update the names and description based on the details of the score_common_coord_handler.
		class length_score : public aligned_pair_score {
		private:
			friend class boost::serialization::access;

			template<class archive> void serialize(archive &ar,
			                                       const size_t /*version*/
			                                       ) {
				ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( aligned_pair_score );
				ar & BOOST_SERIALIZATION_NVP( length_getter_ptr );
			}

			/// \brief TODOCUMENT
			common::clone_ptr<length_getter> length_getter_ptr;

			virtual std::unique_ptr<aligned_pair_score> do_clone() const override final;

			virtual boost::logic::tribool do_higher_is_better() const override final;
			virtual score_value do_calculate(const align::alignment &,
			                                 const protein &,
			                                 const protein &) const override final;
			virtual std::string do_description() const override final;
			virtual std::string do_id_name() const override final;
			virtual str_bool_pair_vec do_short_name_suffixes() const override final;
			virtual std::string do_long_name() const override final;

//			virtual std::unique_ptr<aligned_pair_score> do_build_from_short_name_spec(const std::string &) const override final;

			virtual bool do_less_than_with_same_dynamic_type(const aligned_pair_score &) const override final;

		public:
			explicit length_score(const length_getter &);
			virtual ~length_score() noexcept = default;

			const length_getter & get_length_getter() const;

			std::string description_brackets_string() const;
		};

		bool operator<(const length_score &,
		               const length_score &);

	}
}

#endif
