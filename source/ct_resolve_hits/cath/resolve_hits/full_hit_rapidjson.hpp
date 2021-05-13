/// \file
/// \brief The full_hit class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FULL_HIT_RAPIDJSON_HPP
#define _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FULL_HIT_RAPIDJSON_HPP

#include <optional>

#include "cath/common/rapidjson_addenda/to_rapidjson_string.hpp"
#include "cath/resolve_hits/full_hit.hpp"
#include "cath/resolve_hits/full_hit_fns.hpp"
#include "cath/resolve_hits/full_hit_list.hpp"
#include "cath/resolve_hits/full_hit_list_fns.hpp"
#include "cath/resolve_hits/options/spec/crh_segment_spec.hpp"

namespace cath {
	namespace rslv {

		/// \brief Write the specified full_hit to the specified rapidjson_writer
		template <common::json_style Style>
		void write_to_rapidjson(common::rapidjson_writer<Style> &prm_writer,                        ///< The rapidjson_writer to which the full_hit should be written
		                        const full_hit                  &prm_full_hit,                      ///< The full_hit to write
		                        const crh_segment_spec_opt      &prm_segment_spec = ::std::nullopt, ///< An optional crh_segment_spec which can be used for including each full_hit's trimmed boundaries and resolved boundaries
		                        const full_hit_list_opt         &prm_hits         = ::std::nullopt  ///< An optional full_hit_list (from which the specified full_hit is drawn), which can be used for including the full_hit's resolved boundaries
		                        ) {
			prm_writer.start_object();
			prm_writer.write_key_value( full_hit::get_label_name(),       prm_full_hit.get_label()                   );
			prm_writer.write_key_value( full_hit::get_score_name(),       prm_full_hit.get_score()                   );
			prm_writer.write_key_value( full_hit::get_score_type_name(),  to_string( prm_full_hit.get_score_type() ) );
			prm_writer.write_key( full_hit::get_segments_name()   );
			write_to_rapidjson( prm_writer, prm_full_hit.get_segments() );
			if ( prm_segment_spec ) {
				prm_writer.write_key( full_hit::get_trimmed_name() );
				write_to_rapidjson( prm_writer, get_segments( prm_full_hit, prm_segment_spec->get_overlap_trim_spec() ) );
				if ( prm_hits ) {
					prm_writer.write_key( full_hit::get_resolved_name() );
					write_to_rapidjson( prm_writer, get_present_segments( resolve_all_boundaries( prm_full_hit, *prm_hits, *prm_segment_spec ) ) );
				}
			}
			const auto extras_store = prm_full_hit.get_extras_store();
			for (const auto &extra_info_pair : extras_store) {
				prm_writer.write_key( to_string( extra_info_pair.first ) );
				invoke_for_hit_extra_info( [&] (const auto &x) { prm_writer.write_value( x ); }, extra_info_pair );
			}
			prm_writer.end_object();
		}


		/// \brief Write the specified full_hit_list to the specified rapidjson_writer
		template <common::json_style Style>
		void write_to_rapidjson_with_compact_fullhits(common::rapidjson_writer<Style> &prm_writer,        ///< The rapidjson_writer to which the full_hit_list should be written
		                                              const full_hit_list             &prm_full_hit_list, ///< The full_hit_list to write
		                                              const crh_segment_spec_opt      &prm_segment_spec   ///< An optional crh_segment_spec which can be used for including each full_hit's trimmed boundaries and resolved boundaries
		                                              ) {
			prm_writer.start_array();
			for (const auto &the_full_hit : prm_full_hit_list) {
				prm_writer.write_raw_string(
					common::to_rapidjson_string<common::json_style::COMPACT>(
						the_full_hit,
						0,
						prm_segment_spec,
						::std::make_optional( prm_full_hit_list )
					)
				);
			}
			prm_writer.end_array();
		}

	} // namespace rslv
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_RESOLVE_HITS_CATH_RESOLVE_HITS_FULL_HIT_RAPIDJSON_HPP
