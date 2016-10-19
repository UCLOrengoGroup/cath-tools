/// \file
/// \brief The ssap_score class header

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

#ifndef _CATH_TOOLS_SOURCE_SCORE_ALIGNED_PAIR_SCORE_SSAP_SCORE_H
#define _CATH_TOOLS_SOURCE_SCORE_ALIGNED_PAIR_SCORE_SSAP_SCORE_H

#include "common/clone/clone_ptr.h"
#include "common/cpp14/make_unique.h"
#include "score/aligned_pair_score/aligned_pair_score.h"
#include "score/aligned_pair_score/detail/score_common_coord_handler.h"
#include "score/aligned_pair_score/ssap_score/ssap_score_accuracy.h"
#include "score/aligned_pair_score/ssap_score/ssap_score_post_processing.h"
#include "score/length_getter/length_getter_make_clone.h"
#include "score/length_getter/length_of_longer_getter.h"
#include "ssap/distance_score_formula.h"

namespace cath {
	namespace score {

		/// \brief Calculate (and represent) SSAP score, a measure that scores an aligned pair
		///        by rewarding all pairs of alignment positions that correspond to two residues in
		///        one structure with a similar "view" to the equivalent two residues in the other.
		///        (where a view is the position of the second residue relative to the position/orientation
		///         of the first's backbone atoms)
		///
		/// Suitable length getters include:
		///  * geometric_mean_length_getter
		///  * length_of_longer_getter
		///  * length_of_shorter_getter
		///  * mean_length_getter
		///  * num_aligned_length_getter
		class ssap_score : public aligned_pair_score {
		public:
			/// \brief TODOCUMENT
			static constexpr ssap_score_post_processing default_post_processing       = ssap_score_post_processing::COMPLX_NORMLS_THEN_LOG;

			/// \brief TODOCUMENT
			static constexpr ssap_score_accuracy        default_accuracy              = ssap_score_accuracy::LOW;

			/// \brief TODOCUMENT
			static constexpr size_t                     default_num_excluded_on_sides = 5;

			/// \brief TODOCUMENT
			static constexpr distance_score_formula     default_distance_formula      = distance_score_formula::USED_IN_PREVIOUS_CODE;

		private:
			friend class boost::serialization::access;

			/// \brief A minimum for the score of adding up the quads and subtracting the gap penalty
			///
			/// The score is floored at this extremely low value because a value of 0 or lower will
			/// break the log that (frequently) comes shortly afterwards.
			static constexpr score_value min_total_score = 1e-30;

			template<class archive> void serialize(archive &ar,
			                                       const size_t /*version*/
			                                       ) {
				ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( aligned_pair_score );
				ar & BOOST_SERIALIZATION_NVP( the_coord_handler     );
				ar & BOOST_SERIALIZATION_NVP( length_getter_ptr     );
				ar & BOOST_SERIALIZATION_NVP( post_processing       );
				ar & BOOST_SERIALIZATION_NVP( accuracy              );
				ar & BOOST_SERIALIZATION_NVP( num_excluded_on_sides );
			}

			/// \brief TODOCUMENT
			///
			/// \todo It doesn't look like this is used as yet
			detail::score_common_coord_handler the_coord_handler;

			/// \brief The method to use for generating the normalisation length
			common::clone_ptr<const length_getter> length_getter_ptr = { common::make_unique<length_of_longer_getter>() };

			/// \brief How to post process the basic score once it has been calculated
			ssap_score_post_processing post_processing = default_post_processing;

			/// \brief Whether to use high accuracy (floating point numbers; no explicit rounding) or low accuracy (ints)
			ssap_score_accuracy accuracy = default_accuracy;

			/// \brief The number of residues on each side that are excluded from view calculations
			size_t num_excluded_on_sides = default_num_excluded_on_sides;

			/// \brief TODOCUMENT
			distance_score_formula distance_formula = default_distance_formula;

			static score_size_pair calculate_total_score_and_num_quads(const align::alignment &,
			                                                           const protein &,
			                                                           const protein &,
			                                                           const ssap_score_accuracy &,
			                                                           const size_t &,
																	   const distance_score_formula &);
			static score_value log_score_copy(const score_value &);
			static score_value simple_normalise(const score_value &,
			                                    const size_t &,
			                                    const size_t &);
			static score_value complex_normalise(const score_value &,
			                                     const size_t &,
			                                     const size_t &,
			                                     const size_t &);

			virtual std::unique_ptr<aligned_pair_score> do_clone() const override final;

			virtual boost::logic::tribool do_higher_is_better() const override final;
			virtual score_value do_calculate(const align::alignment &,
			                                 const protein &,
			                                 const protein &) const override final;
			virtual std::string do_description() const override final;
			virtual std::string do_id_name() const override final;
			virtual str_bool_pair_vec do_short_name_suffixes() const override final;
			virtual std::string do_long_name() const override final;
			virtual std::string do_reference() const override final;

//			virtual std::unique_ptr<aligned_pair_score> do_build_from_short_name_spec(const std::string &) const override final;

			virtual bool do_less_than_with_same_dynamic_type(const aligned_pair_score &) const override final;

		public:
			explicit ssap_score(const ssap_score_post_processing &arg_post_processing       = default_post_processing,
			                    const ssap_score_accuracy        &arg_accuracy              = default_accuracy,
			                    const size_t                     &arg_num_excluded_on_sides = default_num_excluded_on_sides);

			explicit ssap_score(const length_getter &,
			                    const ssap_score_post_processing &arg_post_processing       = default_post_processing,
			                    const ssap_score_accuracy        &arg_accuracy              = default_accuracy,
			                    const size_t                     &arg_num_excluded_on_sides = default_num_excluded_on_sides);

			ssap_score(const length_getter &,
			           const align::common_residue_selection_policy &,
			           const align::common_atom_selection_policy &,
			           const ssap_score_post_processing &arg_post_processing       = default_post_processing,
			           const ssap_score_accuracy        &arg_accuracy              = default_accuracy,
			           const size_t                     &arg_num_excluded_on_sides = default_num_excluded_on_sides,
					   const distance_score_formula     &arg_distance_formula      = default_distance_formula);

			virtual ~ssap_score() noexcept = default;

			const length_getter & get_length_getter() const;
			const ssap_score_post_processing & get_post_processing() const;
			const ssap_score_accuracy & get_accuracy() const;
			const size_t & get_num_excluded_on_sides() const;

			const detail::score_common_coord_handler & get_score_common_coord_handler() const;
		};

		bool operator<(const ssap_score &,
		               const ssap_score &);

	} // namespace score
} // namespace cath
#endif
