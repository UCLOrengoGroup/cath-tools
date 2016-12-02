/// \file
/// \brief The sequence_similarity_score class header

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

#ifndef _CATH_TOOLS_SOURCE_SCORE_ALIGNED_PAIR_SCORE_SEQUENCE_SIMILARITY_SCORE_H
#define _CATH_TOOLS_SOURCE_SCORE_ALIGNED_PAIR_SCORE_SEQUENCE_SIMILARITY_SCORE_H

#include "common/clone/clone_ptr.hpp"
#include "common/cpp14/make_unique.hpp"
#include "score/aligned_pair_score/aligned_pair_score.hpp"
#include "score/aligned_pair_score/substitution_matrix/identity_substitution_matrix.hpp"
#include "score/aligned_pair_score/substitution_matrix/substitution_matrix.hpp"
#include "score/length_getter/length_getter_make_clone.hpp"
#include "score/length_getter/length_of_longer_getter.hpp"
#include "structure/protein/amino_acid.hpp"
#include "structure/structure_type_aliases.hpp"

#include <memory>

namespace cath { namespace score { class length_getter; } }

namespace cath {
	namespace score {

		/// \brief Concrete aligned_pair_score to calculate the sequence similarity between an aligned pair
		///
		/// This can be configured on:
		///  * the length over which the score is normalised (eg first, second, longest, shortest, mean, geometric_mean, aligned etc)
		///  * the substitution matrix to use to calculate the score (eg identity, dayhoff, pam70, blosum62)
		class sequence_similarity_score final : public aligned_pair_score {
		private:
//			friend class boost::serialization::access;

//			template<class archive> void serialize(archive &ar,
//			                                       const size_t /*version*/
//			                                       ) {
//				ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( aligned_pair_score );
//				ar & BOOST_SERIALIZATION_NVP( scores );
//				ar & BOOST_SERIALIZATION_NVP( length_getter_ptr );
//			}

			/// \brief The substitution matrix with which the sequence similarity should be scored
			substitution_matrix scores;

			/// \brief The length getter with which the normalisation length should be acquired
			common::clone_ptr<const length_getter> length_getter_ptr = { common::make_unique<length_of_longer_getter>() };

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
			explicit sequence_similarity_score(const substitution_matrix & = make_subs_matrix_identity());
			sequence_similarity_score(const length_getter &,
			                          const substitution_matrix & = make_subs_matrix_identity());
			
			const substitution_matrix & get_substitution_matrix() const;
			const length_getter & get_length_getter() const;
		};

		bool operator<(const sequence_similarity_score &,
		               const sequence_similarity_score &);

	} // namespace score
} // namespace cath
#endif
