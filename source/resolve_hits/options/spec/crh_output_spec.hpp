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

		/// \brief Represent the different cath-resolve-hits output formats
		enum class crh_out_format {
			STANDARD,
			SUMMARY,
			HTML,
			JSON
		};

		/// \brief Type alias for an optional crh_out_format
		using crh_out_format_vec = std::vector<crh_out_format>;

		/// \brief Type alias for an optional crh_out_format
		using crh_out_format_opt = boost::optional<crh_out_format>;

		/// \brief Type alias for a vector of crh_out_format_opt
		using crh_out_format_opt_vec = std::vector<crh_out_format_opt>;

		std::string to_string(const crh_out_format &);

		/// \brief Specify the output for cath-resolve-hits
		class crh_output_spec final {
		private:
			/// \brief The output file to which data should be written
			path_opt            output_file;

			/// \brief Whether to output the hits starts/stops *after* trimming
			hit_boundary_output boundary_output           = DEFAULT_BOUNDARY_OUTPUT;

			/// \brief Whether to output a summary of the input data
			bool                summarise                 = DEFAULT_SUMMARISE;

			/// \brief Whether to output HTML describing the hits and the results
			bool                generate_html_output      = DEFAULT_GENERATE_HTML_OUTPUT;

			/// \brief Whether to output the results in JSON format
			bool                json_output               = DEFAULT_JSON_OUTPUT;

			/// \brief An optional file to which the cath-resolve-hits CSS should be dumped
			path_opt            export_css_file;

			/// \brief Whether to output a summary of the hmmsearch output alignment
			bool                output_hmmsearch_aln      = DEFAULT_OUTPUT_HMMSEARCH_ALN;

		public:
			/// \brief The default value for whether to output the hits starts/stops *after* trimming
			static constexpr hit_boundary_output DEFAULT_BOUNDARY_OUTPUT           = hit_boundary_output::ORIG;

			/// \brief The default value for whether to output a summary of the input data
			static constexpr bool                DEFAULT_SUMMARISE                 = false;

			/// \brief The default value for whether to output HTML describing the hits and the results
			static constexpr bool                DEFAULT_GENERATE_HTML_OUTPUT      = false;

			/// \brief The default value for whether to output the results in JSON format
			static constexpr bool                DEFAULT_JSON_OUTPUT      = false;

			/// \brief The default value for whether to output a summary of the hmmsearch output alignment
			static constexpr bool                DEFAULT_OUTPUT_HMMSEARCH_ALN      = false;

			const path_opt & get_output_file() const;
			const hit_boundary_output & get_boundary_output() const;
			const bool & get_summarise() const;
			const bool & get_generate_html_output() const;
			const bool & get_json_output() const;
			const path_opt & get_export_css_file() const;
			const bool & get_output_hmmsearch_aln() const;

			crh_output_spec & set_output_file(const boost::filesystem::path &);
			crh_output_spec & set_boundary_output(const hit_boundary_output &);
			crh_output_spec & set_summarise(const bool &);
			crh_output_spec & set_generate_html_output(const bool &);
			crh_output_spec & set_json_output(const bool &);
			crh_output_spec & set_export_css_file(const path_opt &);
			crh_output_spec & set_output_hmmsearch_aln(const bool &);
		};

		crh_out_format get_out_format(const crh_output_spec &);

		str_opt get_invalid_description(const crh_output_spec &);

		crh_output_spec & set_output_trimmed_hits(crh_output_spec &,
		                                          const bool &);

	} // namespace rslv
} // namespace cath

#endif
