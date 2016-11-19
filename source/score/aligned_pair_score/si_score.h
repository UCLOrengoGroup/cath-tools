/// \file
/// \brief The si_score class header

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

#ifndef _CATH_TOOLS_SOURCE_SCORE_ALIGNED_PAIR_SCORE_SI_SCORE_H
#define _CATH_TOOLS_SOURCE_SCORE_ALIGNED_PAIR_SCORE_SI_SCORE_H

#include "common/clone/clone_ptr.h"
#include "common/cpp14/make_unique.h"
#include "score/aligned_pair_score/aligned_pair_score.h"
#include "score/aligned_pair_score/rmsd_score.h"
#include "score/length_getter/length_getter_make_clone.h"
#include "score/length_getter/length_of_shorter_getter.h"
#include "score/length_getter/num_aligned_length_getter.h"

namespace cath {
	namespace score {

		/// \brief Calculate (and represent) similarity index (SI), a measure that attempts to
		///        balance the RMSD according to fraction of residues that have been aligned
		///
		/// \todo The mi_score and si_score classes should probably be merged into one
		class si_score final : public  aligned_pair_score {
		private:

//			friend class boost::serialization::access;

//			template<class archive> void serialize(archive &ar,
//			                                       const size_t /*version*/
//			                                       ) {
//				ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( aligned_pair_score );
//				ar & BOOST_SERIALIZATION_NVP( rmsd );
//				ar & BOOST_SERIALIZATION_NVP( num_aligned_getter );
//			}

			/// \brief TODOCUMENT
			rmsd_score rmsd;

			/// \brief TODOCUMENT
			num_aligned_length_getter num_aligned_getter;

			/// \brief TODOCUMENT
			common::clone_ptr<const sym_protein_only_length_getter> full_length_getter_ptr = { common::make_unique<length_of_shorter_getter>() };

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
			si_score() = default;
			explicit si_score(const sym_protein_only_length_getter &);
			si_score(const sym_protein_only_length_getter &,
			         const align::common_residue_selection_policy &,
			         const align::common_atom_selection_policy &);

			const rmsd_score & get_rmsd() const;
			const num_aligned_length_getter & get_num_aligned_getter() const;
			const sym_protein_only_length_getter & get_full_length_getter() const;
		};

		si_score make_simax_score(const align::common_residue_selection_policy & = align::common_residue_select_all_policy(),
		                          const align::common_atom_selection_policy & = align::common_atom_select_ca_policy() );

		bool operator<(const si_score &,
		               const si_score &);

	} // namespace score
} // namespace cath
#endif
