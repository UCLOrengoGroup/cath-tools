/// \file
/// \brief The clustmap_output_options_block class header

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

#ifndef CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_OPTIONS_OPTIONS_BLOCK_CLUSTMAP_OUTPUT_OPTIONS_BLOCK_HPP
#define CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_OPTIONS_OPTIONS_BLOCK_CLUSTMAP_OUTPUT_OPTIONS_BLOCK_HPP

#include <string_view>

#include "cath/cluster/options/spec/clustmap_output_spec.hpp"
#include "cath/options/options_block/options_block.hpp"

namespace cath::clust {

	/// \brief Define an options_block for options specifying how cath-map-clusters should write the output
	class clustmap_output_options_block final : public opts::options_block {
	private:
		using super = opts::options_block;

		// /// \brief The spec this options_block configures
		clustmap_output_spec the_spec;

		[[nodiscard]] std::unique_ptr<opts::options_block> do_clone() const final;
		[[nodiscard]] std::string                          do_get_block_name() const final;
		void do_add_visible_options_to_description(boost::program_options::options_description &,
		                                           const size_t &) final;
		[[nodiscard]] str_opt do_invalid_string( const boost::program_options::variables_map & ) const final;
		[[nodiscard]] str_view_vec do_get_all_options_names() const final;

	  public:
		[[nodiscard]] const clustmap_output_spec &get_clustmap_output_spec() const;

		/// \brief The option name for a batch ID to optionally append as a last column to the output
		static constexpr ::std::string_view PO_APPEND_BATCH_ID{ "append-batch-id" };

		/// \brief The option name for an optional file to which output should be redirected
		static constexpr ::std::string_view PO_OUTPUT_TO_FILE{ "output-to-file" };

		/// \brief The option name for an optional file to which a Markdown summary should be written
		static constexpr ::std::string_view PO_SUMMARISE_TO_FILE{ "summarise-to-file" };

		/// \brief The option name for whether to print out per-domain mapping info
		static constexpr ::std::string_view PO_PRINT_DOMAIN_MAPPING{ "print-entry-results" };
	};

} // namespace cath::clust

#endif // CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_OPTIONS_OPTIONS_BLOCK_CLUSTMAP_OUTPUT_OPTIONS_BLOCK_HPP
