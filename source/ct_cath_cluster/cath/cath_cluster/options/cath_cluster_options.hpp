/// \file
/// \brief The cath_cluster_options class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_CATH_CLUSTER_OPTIONS_HPP
#define _CATH_TOOLS_SOURCE_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_CATH_CLUSTER_OPTIONS_HPP

#include "cath/cath_cluster/options/options_block/cath_cluster_clustering_options_block.hpp"
#include "cath/cath_cluster/options/options_block/cath_cluster_input_options_block.hpp"
#include "cath/cath_cluster/options/options_block/cath_cluster_output_options_block.hpp"
#include "cath/options/executable/executable_options.hpp"

// clang-format off
namespace cath::clust {class cath_cluster_spec; }
// clang-format on

namespace cath::clust {

	/// \brief Implement the executable_options for cath-resolve-hits
	class cath_cluster_options final : public opts::executable_options {
	private:
		using super = opts::executable_options;

		/// \brief The input options_block
		cath_cluster_input_options_block  the_input_ob;

		/// \brief The clustering options_block
		cath_cluster_clustering_options_block clustering_ob;

		/// \brief The output options_block
		cath_cluster_output_options_block the_output_ob;

		[[nodiscard]] std::string_view                         do_get_program_name() const final;
		boost::program_options::positional_options_description get_positional_options() final;
		[[nodiscard]] str_opt                                  do_get_error_or_help_string() const final;

		[[nodiscard]] std::string do_get_help_prefix_string() const final;
		[[nodiscard]] std::string do_get_help_suffix_string() const final;
		[[nodiscard]] std::string do_get_overview_string() const final;

	  public:
		cath_cluster_options();

		[[nodiscard]] const cath_cluster_input_spec &     get_cath_cluster_input_spec() const;
		[[nodiscard]] const cath_cluster_clustering_spec &get_cath_cluster_clustering_spec() const;
		[[nodiscard]] const cath_cluster_output_spec &    get_cath_cluster_output_spec() const;

		/// The name of the program that uses this executable_options
		static constexpr ::std::string_view PROGRAM_NAME{ "cath-cluster" };
	};

} // namespace cath::clust

#endif // _CATH_TOOLS_SOURCE_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_CATH_CLUSTER_OPTIONS_HPP
