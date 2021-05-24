/// \file
/// \brief The cath_id_score_category header

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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FILE_CATH_ID_SCORE_CATEGORY_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FILE_CATH_ID_SCORE_CATEGORY_HPP

#include <boost/utility/string_ref_fwd.hpp>

namespace cath::rslv {

	/// \brief Represent the type of match ID, wrt the CATH-specific scoring that Jon uses
	enum class cath_id_score_category : char {
		NORMAL, ///< A normal match (use this as the default if unsure)
		DC_TYPE ///< An ID like dc_c869189e57e572c71376c2f3dfe7dc9c that is handled differently
	};

	cath_id_score_category cath_score_category_of_id(const boost::string_ref &,
	                                                 const bool &);

	std::string to_string(const cath_id_score_category &);

	std::ostream & operator<<(std::ostream &,
	                          const cath_id_score_category &);

	/// \brief Return whether the specified HMMER evalues are suspicious under the CATH-Gene3D rules
	inline constexpr bool hmmer_evalues_are_suspicious(const double &prm_cond_evalue, ///< The HMMER conditional evalue
	                                                   const double &prm_indp_evalue  ///< The HMMER independent evalue
	                                                   ) {
		// constexpr double EVALUE_CUTOFF = 0.001; //< Not permitted by GCC 4.9.2
		return (
			prm_cond_evalue <= 0.001
			&&
			prm_indp_evalue >  0.001
		);
	}

	/// \brief Return the value by which the bitscore should be divided under the CATH-Gene3D rules
	inline constexpr double bitscore_divisor(const bool                   &prm_apply_cath_policies, ///< Whether the CATH-Gene3D rules are to be applied (if not, then this will return 1.0 )
	                                         const bool                   &prm_evalues_are_susp     ///< Whether the HMMER evalues were deemed suspicious
	                                         ) {
		return ( prm_apply_cath_policies && prm_evalues_are_susp ) ? 4.0 : 1.0;
	}

} // namespace cath::rslv

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FILE_CATH_ID_SCORE_CATEGORY_HPP
