/// \file
/// \brief The structal_score class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_STRUCTAL_SCORE_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_STRUCTAL_SCORE_HPP

#include "cath/common/clone/clone_ptr.hpp"
#include "cath/common/cpp14/make_unique.hpp"
#include "cath/score/aligned_pair_score/aligned_pair_score.hpp"
#include "cath/score/aligned_pair_score/rmsd_score.hpp"
#include "cath/score/length_getter/length_of_shorter_getter.hpp"
#include "cath/score/length_getter/num_aligned_length_getter.hpp"
#include "cath/structure/structure_type_aliases.hpp"

namespace cath::score {

	/// \brief Calculate (and represent) match index (MI), a measure that attempts to
	///        balance the RMSD according to fraction of residues that have been aligned
	///
	/// \todo The structal_score and si_score classes should probably be merged into one
	class structal_score : public aligned_pair_score {
	private:
		/// \brief TODOCUMENT
		detail::score_common_coord_handler the_coord_handler;

		[[nodiscard]] std::unique_ptr<aligned_pair_score> do_clone() const final;

		[[nodiscard]] boost::logic::tribool do_higher_is_better() const final;
		[[nodiscard]] score_value do_calculate( const align::alignment &, const protein &, const protein & ) const final;
		[[nodiscard]] std::string       do_description() const final;
		[[nodiscard]] std::string       do_id_name() const final;
		[[nodiscard]] str_bool_pair_vec do_short_name_suffixes() const final;
		[[nodiscard]] std::string       do_long_name() const final;
		[[nodiscard]] std::string       do_reference() const final;

		// std::unique_ptr<aligned_pair_score> do_build_from_short_name_spec(const std::string &) const final;

		[[nodiscard]] bool do_less_than_with_same_dynamic_type( const aligned_pair_score & ) const final;

		static score_value score_for_target_length(const std::pair<geom::coord_list_vec, geom::coord_list_vec> &,
		                                           const score_value &);

	public:
		structal_score() = default;
		structal_score(const align::common_residue_selection_policy &,
		               const align::common_atom_selection_policy &);

		[[nodiscard]] const detail::score_common_coord_handler &get_score_common_coord_handler() const;
	};

	bool operator<(const structal_score &,
	               const structal_score &);

} // namespace cath::score

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_STRUCTAL_SCORE_HPP
