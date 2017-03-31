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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_FULL_HIT_RAPIDJSON_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_FULL_HIT_RAPIDJSON_H

#include "common/rapidjson_addenda/to_rapidjson_string.hpp"
#include "resolve_hits/full_hit.hpp"
#include "resolve_hits/full_hit_fns.hpp"
#include "resolve_hits/full_hit_list.hpp"
#include "resolve_hits/full_hit_list_fns.hpp"
#include "resolve_hits/options/spec/crh_segment_spec.hpp"

namespace cath {
	namespace rslv {

		/// \brief Write the specified full_hit to the specified rapidjson_writer
		template <common::json_style Style>
		void write_to_rapidjson(common::rapidjson_writer<Style> &arg_writer,                     ///< The rapidjson_writer to which the full_hit should be written
		                        const full_hit                  &arg_full_hit,                   ///< The full_hit to write
		                        const crh_segment_spec_opt      &arg_segment_spec = boost::none, ///< An optional crh_segment_spec which can be used for including each full_hit's trimmed boundaries and resolved boundaries
		                        const full_hit_list_opt         &arg_hits         = boost::none  ///< An optional full_hit_list (from which the specified full_hit is drawn), which can be used for including the full_hit's resolved boundaries
		                        ) {
			arg_writer.start_object();
			arg_writer.write_key_value( full_hit::get_label_name(),       arg_full_hit.get_label()                   );
			arg_writer.write_key_value( full_hit::get_score_name(),       arg_full_hit.get_score()                   );
			arg_writer.write_key_value( full_hit::get_score_type_name(),  to_string( arg_full_hit.get_score_type() ) );
			arg_writer.write_key( full_hit::get_segments_name()   );
			write_to_rapidjson( arg_writer, arg_full_hit.get_segments() );
			if ( arg_segment_spec ) {
				arg_writer.write_key( full_hit::get_trimmed_name() );
				write_to_rapidjson( arg_writer, get_segments( arg_full_hit, arg_segment_spec->get_overlap_trim_spec() ) );
				if ( arg_hits ) {
					arg_writer.write_key( full_hit::get_resolved_name() );
					write_to_rapidjson( arg_writer, get_present_segments( resolve_all_boundaries( arg_full_hit, *arg_hits, *arg_segment_spec ) ) );
				}
			}
			const auto extras_store = arg_full_hit.get_extras_store();
			for (const auto &extra_info_pair : extras_store) {
				arg_writer.write_key( to_string( extra_info_pair.first ) );
				invoke_for_hit_extra_info( [&] (const auto &x) { arg_writer.write_value( x ); }, extra_info_pair );
			}
			arg_writer.end_object();
		}


		/// \brief Write the specified full_hit_list to the specified rapidjson_writer
		template <common::json_style Style>
		void write_to_rapidjson_with_compact_fullhits(common::rapidjson_writer<Style> &arg_writer,        ///< The rapidjson_writer to which the full_hit_list should be written
		                                              const full_hit_list             &arg_full_hit_list, ///< The full_hit_list to write
		                                              const crh_segment_spec_opt      &arg_segment_spec   ///< An optional crh_segment_spec which can be used for including each full_hit's trimmed boundaries and resolved boundaries
		                                              ) {
			arg_writer.start_array();
			for (const auto &the_full_hit : arg_full_hit_list) {
				arg_writer.write_raw_string(
					common::to_rapidjson_string<common::json_style::COMPACT>(
						the_full_hit,
						0,
						arg_segment_spec,
						boost::make_optional( arg_full_hit_list )
					)
				);
			}
			arg_writer.end_array();
		}

	} // namespace rslv
} // namespace cath

#endif
