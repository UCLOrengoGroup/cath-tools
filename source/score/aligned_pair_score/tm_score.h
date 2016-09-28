/// \file
/// \brief The tm_score class header

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

#ifndef TM_SCORE_H_INCLUDED
#define TM_SCORE_H_INCLUDED

#include "common/clone/clone_ptr.h"
#include "common/cpp14/make_unique.h"
#include "score/aligned_pair_score/aligned_pair_score.h"
#include "score/aligned_pair_score/rmsd_score.h"
#include "score/length_getter/length_of_shorter_getter.h"
#include "score/length_getter/num_aligned_length_getter.h"
#include "structure/structure_type_aliases.h"

namespace cath {
	namespace score {

		/// \brief Calculate (and represent) match index (MI), a measure that attempts to
		///        balance the RMSD according to fraction of residues that have been aligned
		///
		/// \todo The tm_score and si_score classes should probably be merged into one
		class tm_score : public aligned_pair_score {
		private:
			/// \brief TODOCUMENT
			detail::score_common_coord_handler the_coord_handler;

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

			static score_value score_for_target_length(const std::pair<geom::coord_list_vec, geom::coord_list_vec> &,
			                                           const score_value &);

		public:
			tm_score() = default;
			tm_score(const align::common_residue_selection_policy &,
			         const align::common_atom_selection_policy &);
			virtual ~tm_score() noexcept = default;

			const detail::score_common_coord_handler & get_score_common_coord_handler() const;
		};

		bool operator<(const tm_score &,
		               const tm_score &);

	}
}
#endif
