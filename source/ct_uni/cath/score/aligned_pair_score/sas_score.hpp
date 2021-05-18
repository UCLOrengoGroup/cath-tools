/// \file
/// \brief The sas_score class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_SAS_SCORE_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_SAS_SCORE_HPP

#include "cath/score/aligned_pair_score/aligned_pair_score.hpp"
#include "cath/score/aligned_pair_score/length_score.hpp"
#include "cath/score/aligned_pair_score/rmsd_score.hpp"
#include "cath/score/length_getter/num_aligned_length_getter.hpp"

namespace cath {
	namespace score {

		/// \brief Calculate (and represent) structural alignment score (SAS), a measure that attempts
		///        to balance the RMSD based on the number of aligned residues.
		class sas_score : public aligned_pair_score {
		private:
			friend class boost::serialization::access;

			template<class archive> void serialize(archive &ar,
			                                       const size_t /*version*/
			                                       ) {
				ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( aligned_pair_score );
				ar & BOOST_SERIALIZATION_NVP( rmsd );
				ar & BOOST_SERIALIZATION_NVP( num_aligned_residues );
			}

			/// \brief TODOCUMENT
			rmsd_score   rmsd;

			/// \brief TODOCUMENT
			length_score num_aligned_residues { num_aligned_length_getter() };

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

		  public:
			/// \brief Default ctor for sas_score that uses the defaults for rmsd_score and length_score (selecting CA atoms for all aligned residues)
			sas_score() = default;
			sas_score(const align::common_residue_selection_policy &,
			          const align::common_atom_selection_policy &);

			[[nodiscard]] const rmsd_score &  get_rmsd_score() const;
			[[nodiscard]] const length_score &get_num_aligned_residues() const;
		};

		bool operator<(const sas_score &,
		               const sas_score &);

	} // namespace score
} // namespace cath
#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_SAS_SCORE_HPP
