/// \file
/// \brief The full_hit functions header

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

#ifndef CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FULL_HIT_FNS_HPP
#define CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FULL_HIT_FNS_HPP

#include <optional>

#include "cath/resolve_hits/full_hit.hpp"
#include "cath/resolve_hits/full_hit_list.hpp"
#include "cath/resolve_hits/hit_output_format.hpp"
#include "cath/resolve_hits/options/spec/crh_score_spec.hpp"
#include "cath/resolve_hits/options/spec/crh_segment_spec.hpp"
#include "cath/resolve_hits/options/spec/hit_boundary_output.hpp"
#include "cath/resolve_hits/score_functions.hpp"
#include "cath/resolve_hits/trim/trim_spec.hpp"

namespace cath::rslv {

	seq::seq_seg_vec get_segments(const full_hit &,
	                              const trim_spec_opt & = ::std::nullopt);

	std::string get_segments_string(const full_hit &,
	                                const trim_spec_opt & = ::std::nullopt);

	cath::str_vec get_field_headers(const full_hit &,
	                                const bool &,
	                                const bool &);

	std::string to_string(const full_hit &,
	                      const hit_output_format & = hit_output_format::CLASS,
	                      const std::string & = std::string{},
	                      const crh_segment_spec_opt & = ::std::nullopt,
	                      const full_hit_list_opt & = ::std::nullopt,
	                      const hit_boundary_output & = hit_boundary_output::ORIG);

	/// \brief Get the crh-score associated with the specified hit, calculated according to the specified crh_score_spec
	inline resscr_t get_crh_score(const full_hit       &prm_full_hit,  ///< The full_hit to query
	                              const crh_score_spec &prm_score_spec ///< The crh_score_spec to use in any score calculations that are required
	                              ) {
		switch ( prm_full_hit.get_score_type() ) {
			case( hit_score_type::FULL_EVALUE ) : {
				return crh_score_of_evalue(
					prm_full_hit.get_score(),
					get_total_length( prm_full_hit ),
					prm_score_spec
				);
			}
			case( hit_score_type::BITSCORE    ) : {
				return crh_score_of_pseudo_bitscore(
					prm_full_hit.get_score(),
					get_total_length( prm_full_hit ),
					prm_score_spec
				);
			}
			case( hit_score_type::CRH_SCORE   ) : {
				return debug_numeric_cast<resscr_t>( prm_full_hit.get_score() );
			}
		}
		BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of hit_score_type not recognised whilst getting crh_score from full_hit"));
	}

} // namespace cath::rslv

#endif // CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FULL_HIT_FNS_HPP
