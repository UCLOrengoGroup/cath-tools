/// \file
/// \brief The crh_spec class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_OPTIONS_SPEC_CRH_SPEC_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_OPTIONS_SPEC_CRH_SPEC_H

#include "resolve_hits/options/spec/crh_filter_spec.hpp"
#include "resolve_hits/options/spec/crh_html_spec.hpp"
#include "resolve_hits/options/spec/crh_input_spec.hpp"
#include "resolve_hits/options/spec/crh_output_spec.hpp"
#include "resolve_hits/options/spec/crh_score_spec.hpp"
#include "resolve_hits/options/spec/crh_segment_spec.hpp"

namespace cath {
	namespace rslv {

		/// \brief Wrap the different specs involved in cath-resolve-hits
		class crh_spec final {
		private:
			/// \brief The input spec
			crh_input_spec   the_input_spec;

			/// \brief The segment spec
			crh_segment_spec the_segment_spec;

			/// \brief The score spec
			crh_score_spec   the_score_spec;

			/// \brief The filter spec
			crh_filter_spec  the_filter_spec;

			/// \brief The output spec
			crh_output_spec  the_output_spec;

			/// \brief The html spec
			crh_html_spec    the_html_spec;

		public:
			crh_spec() = default;
			explicit crh_spec(const crh_input_spec &,
			                  const crh_segment_spec & = crh_segment_spec{},
			                  const crh_score_spec & = crh_score_spec{},
			                  const crh_filter_spec & = crh_filter_spec{},
			                  const crh_output_spec & = crh_output_spec{},
			                  const crh_html_spec & = crh_html_spec{});

			crh_input_spec & get_input_spec();
			const crh_input_spec & get_input_spec() const;
			crh_segment_spec & get_segment_spec();
			const crh_segment_spec & get_segment_spec() const;
			crh_score_spec & get_score_spec();
			const crh_score_spec & get_score_spec() const;
			crh_filter_spec & get_filter_spec();
			const crh_filter_spec & get_filter_spec() const;
			crh_output_spec & get_output_spec();
			const crh_output_spec & get_output_spec() const;
			crh_html_spec & get_html_spec();
			const crh_html_spec & get_html_spec() const;

			crh_spec & set_input_spec(const crh_input_spec &);
			crh_spec & set_segment_spec(const crh_segment_spec &);
			crh_spec & set_score_spec(const crh_score_spec &);
			crh_spec & set_filter_spec(const crh_filter_spec &);
			crh_spec & set_output_spec(const crh_output_spec &);
			crh_spec & set_html_spec(const crh_html_spec &);
		};


	} // namespace rslv
} // namespace cath

#endif
