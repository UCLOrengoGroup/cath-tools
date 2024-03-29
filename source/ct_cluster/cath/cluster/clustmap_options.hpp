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

#ifndef CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_CLUSTMAP_OPTIONS_HPP
#define CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_CLUSTMAP_OPTIONS_HPP

#include <iosfwd>
#include <string_view>

#include "cath/cluster/options/options_block/clust_mapping_options_block.hpp"
#include "cath/cluster/options/options_block/clustmap_input_options_block.hpp"
#include "cath/cluster/options/options_block/clustmap_output_options_block.hpp"
#include "cath/options/executable/executable_options.hpp"
#include "cath/options/options_block/detail_help_options_block.hpp"

// clang-format off
namespace cath::clust {class clustmap_spec; }
// clang-format on

namespace cath::clust {

	/// \brief Implement the executable_options for cath-resolve-hits
	class clustmap_options final : public opts::executable_options {
	private:
		using super = opts::executable_options;

		static std::map<std::string, str_str_pair> detail_help_spec();

		/// \brief The cath-resolve-hits input options_block
		clustmap_input_options_block         the_input_ob;

		/// \brief The cath-resolve-hits segment options_block
		clust_mapping_options_block          the_mapping_ob;

		/// \brief The cath-resolve-hits output options_block
		clustmap_output_options_block        the_output_ob;

		/// \brief The detailed help options_block
		opts::detail_help_options_block      the_detail_help_ob;

		[[nodiscard]] std::string_view                         do_get_program_name() const final;
		boost::program_options::positional_options_description get_positional_options() final;
		[[nodiscard]] str_opt                                  do_get_error_or_help_string() const final;

		[[nodiscard]] std::string do_get_help_prefix_string() const final;
		[[nodiscard]] std::string do_get_help_suffix_string() const final;
		[[nodiscard]] std::string do_get_overview_string() const final;

	  public:
		clustmap_options();

		[[nodiscard]] clustmap_spec get_clustmap_spec() const;

		[[nodiscard]] const clustmap_input_spec & get_clustmap_input_spec() const;
		[[nodiscard]] const clust_mapping_spec &  get_clust_mapping_spec() const;
		[[nodiscard]] const clustmap_output_spec &get_clustmap_output_spec() const;

		/// The name of the program that uses this executable_options
		static constexpr ::std::string_view PROGRAM_NAME{ "cath-map-clusters" };
	};

	std::string get_cmc_sorting_criteria_help_string();

	// std::string get_clustmap_raw_format_help_string();
	// std::string get_clustmap_cath_rules_help_string();

} // namespace cath::clust

#endif // CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_CLUSTMAP_OPTIONS_HPP
