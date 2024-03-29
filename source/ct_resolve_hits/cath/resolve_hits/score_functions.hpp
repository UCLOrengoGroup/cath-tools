/// \file
/// \brief The score functions header

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

#ifndef CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_SCORE_FUNCTIONS_HPP
#define CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_SCORE_FUNCTIONS_HPP

#include "cath/common/debug_numeric_cast.hpp"
#include "cath/resolve_hits/options/spec/crh_score_spec.hpp"
#include "cath/resolve_hits/resolve_hits_type_aliases.hpp"

#include <cmath>

namespace cath::rslv {

	/// \brief A pseudo bitscore 
	inline double pseudo_bitscore_of_evalue(const double &prm_evalue ///< The evalue from which to convert
	                                        ) {
		return std::log( prm_evalue ) / -0.7;
	}

	/// \brief The value by which scores for hits of the specified length should be multiplied to apply the specified length-preference
	inline double length_multiplier(const size_t &prm_length,           ///< The total length of the hit
	                                const double &prm_length_preference ///< The degree to which long hits are preferred (where 0 is neutral and negative values are allowed)
	                                ) {
		constexpr double longish_length = 400;
		return std::pow( debug_numeric_cast<double>( prm_length ) / longish_length, prm_length_preference );
	}

	/// \brief The value by which scores for hits of the specified length should be multiplied to apply the specified score_spec
	inline double length_multiplier(const size_t         &prm_length,    ///< The total length of the hit
	                                const crh_score_spec &prm_score_spec ///< The crh_score_spec to apply
	                                ) {
		return length_multiplier(
			prm_length,
			prm_score_spec.get_long_domains_preference()
		);
	}

	/// \brief The cath-resolve-hits score for the specified (pseudo)bitscore after applying the specified highscore-preference
	///        but not any length-preference
	inline double length_agnostic_crh_score_of_pseudo_bitscore(const double &prm_pseudo_bitscore,     ///< The pseudo_bitscore from which to convert
	                                                           const double &prm_highscore_preference ///< The degree to which higher-scoring hits are preferred (where 0 uses unaltered scores and negative scores are allowed)
	                                                           ) {
		constexpr double power_base   =    1.41421356237; // sqrt( 2 ) \todo Come a constexpr sqrt in std, use it here : ` = constexpr_sqrt( 2.0 );`
		constexpr double max_bitscore = 2000.0;
		constexpr double max_score    =  100.0;
		return max_score * std::pow( prm_pseudo_bitscore / max_bitscore, std::pow( power_base, prm_highscore_preference ) );
		// return ( prm_pseudo_bitscore > max_bitscore )
			// ? max_score
			// : max_score * std::pow( prm_pseudo_bitscore / max_bitscore, std::pow( power_base, prm_highscore_preference ) );
	}

	/// \brief The cath-resolve-hits score for the specified (pseudo)bitscore after applying the specified highscore-preference
	///        but not any length-preference
	inline double length_agnostic_crh_score_of_pseudo_bitscore(const double         &prm_pseudo_bitscore, ///< The pseudo_bitscore from which to convert
	                                                           const crh_score_spec &prm_score_spec       ///< The crh_score_spec to apply
	                                                           ) {
		return length_agnostic_crh_score_of_pseudo_bitscore(
			prm_pseudo_bitscore,
			prm_score_spec.get_high_scores_preference()
		);
	}

	/// \brief Calculate the cath-resolve-hits score for the specified (pseudo)bitscore given the specified length,
	///        highscore_preference and length_preference
	inline resscr_t crh_score_of_pseudo_bitscore(const double &prm_pseudo_bitscore,      ///< The pseudo_bitscore from which to convert
	                                             const size_t &prm_length,               ///< The total length of the hit
	                                             const double &prm_highscore_preference, ///< The degree to which higher-scoring hits are preferred (where 0 uses unaltered scores and negative scores are allowed)
	                                             const double &prm_length_preference     ///< The degree to which long hits are preferred (where 0 is neutral and negative values are allowed)
	                                             ) {
		return static_cast<resscr_t>(
			length_multiplier( prm_length, prm_length_preference )
			*
			length_agnostic_crh_score_of_pseudo_bitscore( prm_pseudo_bitscore, prm_highscore_preference )
		);
	}

	/// \brief Calculate the cath-resolve-hits score for the specified (pseudo)bitscore given the specified length,
	///        highscore_preference and length_preference
	inline resscr_t crh_score_of_pseudo_bitscore(const double         &prm_pseudo_bitscore, ///< The pseudo_bitscore from which to convert
	                                             const size_t         &prm_length,          ///< The total length of the hit
	                                             const crh_score_spec &prm_score_spec       ///< The crh_score_spec to apply
	                                             ) {
		return crh_score_of_pseudo_bitscore(
			prm_pseudo_bitscore,
			prm_length,
			prm_score_spec.get_high_scores_preference(),
			prm_score_spec.get_long_domains_preference()
		);
	}

	/// \brief Calculate the cath-resolve-hits score for the specified evalue given the specified length,
	///        highscore_preference and length_preference
	inline resscr_t crh_score_of_evalue(const double &prm_evalue,               ///< The evalue from which to convert
	                                    const size_t &prm_length,               ///< The total length of the hit
	                                    const double &prm_highscore_preference, ///< The degree to which higher-scoring hits are preferred (where 0 uses unaltered scores and negative scores are allowed)
	                                    const double &prm_length_preference     ///< The degree to which long hits are preferred (where 0 is neutral and negative values are allowed)
	                                    ) {
		return crh_score_of_pseudo_bitscore(
			pseudo_bitscore_of_evalue( prm_evalue ),
			prm_length,
			prm_highscore_preference,
			prm_length_preference
		);
	}

	/// \brief Calculate the cath-resolve-hits score for the specified evalue given the specified length,
	///        highscore_preference and length_preference
	inline resscr_t crh_score_of_evalue(const double         &prm_evalue,    ///< The evalue from which to convert
	                                    const size_t         &prm_length,    ///< The total length of the hit
	                                    const crh_score_spec &prm_score_spec ///< The crh_score_spec to apply
	                                    ) {
		return crh_score_of_evalue(
			prm_evalue,
			prm_length,
			prm_score_spec.get_high_scores_preference(),
			prm_score_spec.get_long_domains_preference()
		);
	}

} // namespace cath::rslv

#endif // CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_SCORE_FUNCTIONS_HPP
