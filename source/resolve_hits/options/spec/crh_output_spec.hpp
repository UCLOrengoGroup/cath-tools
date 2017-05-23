/// \file
/// \brief The crh_output_spec class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_OPTIONS_SPEC_CRH_OUTPUT_SPEC_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_OPTIONS_SPEC_CRH_OUTPUT_SPEC_H

#include <boost/optional.hpp>

#include "common/path_type_aliases.hpp"
#include "common/type_aliases.hpp"
#include "resolve_hits/options/spec/hit_boundary_output.hpp"

namespace cath {
	namespace rslv {

		/// \brief Specify the output for cath-resolve-hits
		class crh_output_spec final {
		private:
			/// \brief Any files to which hits text should be output
			path_vec            hits_text_files;

			/// \brief Whether to suppress the default output of hits text to stdout
			bool                quiet           = DEFAULT_QUIET;

			/// \brief Whether to output the hits starts/stops *after* trimming
			hit_boundary_output boundary_output = DEFAULT_BOUNDARY_OUTPUT;

			/// \brief Any files to which a summary of the input should be output
			path_vec            summarise_files;

			/// \brief Any files to which HTML should be output
			path_vec            html_output_files;

			/// \brief Any files to which JSON should be output
			path_vec            json_output_files;

			/// \brief Any files to which the HTML's CSS should be output
			path_opt            export_css_file;

			/// \brief Whether to output a summary of the hmmsearch output alignment
			bool                output_hmmsearch_aln = DEFAULT_OUTPUT_HMMSEARCH_ALN;

		public:
			/// \brief The default value for whether to suppress the default output of hits text to stdout
			static constexpr bool                DEFAULT_QUIET           = false;

			/// \brief The default value for whether to output the hits starts/stops *after* trimming
			static constexpr hit_boundary_output DEFAULT_BOUNDARY_OUTPUT = hit_boundary_output::ORIG;

			/// \brief The default value for whether to output a summary of the hmmsearch output alignment
			static constexpr bool                DEFAULT_OUTPUT_HMMSEARCH_ALN = false;

			const path_vec & get_hits_text_files() const;
			const bool & get_quiet() const;
			const hit_boundary_output & get_boundary_output() const;
			const path_vec & get_summarise_files() const;
			const path_vec & get_html_output_files() const;
			const path_vec & get_json_output_files() const;
			const path_opt & get_export_css_file() const;
			const bool & get_output_hmmsearch_aln() const;

			crh_output_spec & set_hits_text_files(const path_vec &);
			crh_output_spec & set_quiet(const bool &);
			crh_output_spec & set_boundary_output(const hit_boundary_output &);
			crh_output_spec & set_summarise_files(const path_vec &);
			crh_output_spec & set_html_output_files(const path_vec &);
			crh_output_spec & set_json_output_files(const path_vec &);
			crh_output_spec & set_export_css_file(const path_opt &);
			crh_output_spec & set_output_hmmsearch_aln(const bool &);
		};

		path_vec get_all_output_paths(const crh_output_spec &);

		str_opt get_invalid_description(const crh_output_spec &);

		crh_output_spec & set_output_trimmed_hits(crh_output_spec &,
		                                          const bool &);

	} // namespace rslv
} // namespace cath

#endif
