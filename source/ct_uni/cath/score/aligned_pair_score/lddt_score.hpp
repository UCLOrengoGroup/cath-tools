/// \file
/// \brief The lddt_score class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LDDT_SCORE_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LDDT_SCORE_HPP

#include "cath/score/aligned_pair_score/aligned_pair_score.hpp"
#include "cath/score/aligned_pair_score/detail/score_common_coord_handler.hpp"

namespace cath::score {

	/// \brief Represent the standard distance thresholds that are used in lDDT scores
	enum class lddt_distance_threshold : char {
		DEFAULT_AVG, ///< The average of the results for 0.5, 1.0, 2.0 and 4.0 angstroms
		HALF_A,      ///< 0.5 angstroms
		ONE_A,       ///< 1.0 angstroms
		TWO_A,       ///< 1.0 angstroms
		FOUR_A       ///< 4.0 angstroms
	};

	/// \brief Calculate (and represent) structural alignment score (SAS), a measure that attempts
	///        to balance the RMSD based on the number of aligned residues.
	///
	/// \ingroup cath_score_aligned_pair_score_group
	class lddt_score : public aligned_pair_score {
	private:
		friend class boost::serialization::access;

		template<class archive> void serialize(archive &ar,
		                                       const size_t /*version*/
		                                       ) {
			ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( aligned_pair_score );
			ar & BOOST_SERIALIZATION_NVP( the_coord_handler );
			ar & BOOST_SERIALIZATION_NVP( threshold_values );
		}

		/// \brief TODOCUMENT
		detail::score_common_coord_handler the_coord_handler;

		/// \brief The threshold values over which the results should be averaged
		doub_vec threshold_values;

		// A threshold value that will ensure everything gets considered
		static constexpr double THRESHOLD_FOR_ALL = std::numeric_limits<double>::max();

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

		void init_distance_threshold_values(const lddt_distance_threshold &);

	public:
		lddt_score();
		explicit lddt_score(const lddt_distance_threshold &);
		lddt_score(const lddt_distance_threshold &,
		           const align::common_residue_selection_policy &,
		           const align::common_atom_selection_policy &);

		[[nodiscard]] const doub_vec &get_threshold_values() const;

		[[nodiscard]] const detail::score_common_coord_handler &get_score_common_coord_handler() const;
	};

	bool operator<(const lddt_score &,
	               const lddt_score &);

	std::string get_thresholds_summary_string(const lddt_score &);

	std::vector<lddt_distance_threshold> get_all_lddt_distance_thresholds();

} // namespace cath::score

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_ALIGNED_PAIR_SCORE_LDDT_SCORE_HPP
