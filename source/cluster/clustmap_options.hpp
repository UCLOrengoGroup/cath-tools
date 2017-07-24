/// \file
/// \brief The clustmap_options class header

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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_OPTIONS_CLUSTMAP_OPTIONS_H
#define _CATH_TOOLS_SOURCE_CLUSTER_OPTIONS_CLUSTMAP_OPTIONS_H

#include "cluster/options/options_block/clust_mapping_options_block.hpp"
#include "cluster/options/options_block/clustmap_input_options_block.hpp"
#include "cluster/options/options_block/clustmap_output_options_block.hpp"
#include "options/executable/executable_options.hpp"

#include <iosfwd>

namespace cath { namespace clust {class clustmap_spec; } }

namespace cath {
	namespace clust {

		/// \brief Implement the executable_options for cath-resolve-hits
		class clustmap_options final : public opts::executable_options {
		private:
			using super = opts::executable_options;

			static const std::string STANDARD_USAGE_ERROR_STRING;

			/// \brief The cath-resolve-hits input options_block
			clustmap_input_options_block         the_input_ob;

			/// \brief The cath-resolve-hits segment options_block
			clust_mapping_options_block          the_mapping_ob;

			/// \brief The cath-resolve-hits output options_block
			clustmap_output_options_block        the_output_ob;

			std::string do_get_program_name() const final;
			boost::program_options::positional_options_description get_positional_options() final;
			str_opt do_get_error_or_help_string() const final;

			std::string do_get_help_prefix_string() const final;
			std::string do_get_help_suffix_string() const final;
			std::string do_get_overview_string() const final;

		public:
			clustmap_options();

			clustmap_spec get_clustmap_spec() const;

			const clustmap_input_spec & get_clustmap_input_spec() const;
			const clust_mapping_spec & get_clust_mapping_spec() const;
			const clustmap_output_spec & get_clustmap_output_spec() const;

			static const std::string PROGRAM_NAME;
		};

		// std::string get_clustmap_raw_format_help_string();
		// std::string get_clustmap_cath_rules_help_string();

	} // namespace clust
} // namespace cath

#endif
