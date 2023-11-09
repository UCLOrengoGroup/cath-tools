/// \file
/// \brief The cath_cluster_clustering_options_block class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_OPTIONS_BLOCK_CATH_CLUSTER_CLUSTERING_OPTIONS_BLOCK_HPP
#define _CATH_TOOLS_SOURCE_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_OPTIONS_BLOCK_CATH_CLUSTER_CLUSTERING_OPTIONS_BLOCK_HPP

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>

#include "cath/cath_cluster/options/spec/cath_cluster_clustering_spec.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/options/options_block/options_block.hpp"

namespace cath::clust {

	/// \brief Define an options_block for options specifying how cath-cluster should write the output
	class cath_cluster_clustering_options_block final : public opts::options_block {
	private:
		using super = opts::options_block;

		/// \brief The spec this options_block configures
		cath_cluster_clustering_spec the_spec;

		[[nodiscard]] std::unique_ptr<opts::options_block> do_clone() const final;
		[[nodiscard]] std::string                          do_get_block_name() const final;
		void do_add_visible_options_to_description(boost::program_options::options_description &,
		                                           const size_t &) final;
		[[nodiscard]] str_opt do_invalid_string( const boost::program_options::variables_map & ) const final;
		[[nodiscard]] str_view_vec do_get_all_options_names() const final;

	  public:
		[[nodiscard]] const cath_cluster_clustering_spec &get_cath_cluster_clustering_spec() const;

		/// \brief The option name for the levels at which the clustering should be performed
		static constexpr ::std::string_view PO_LEVELS{ "levels" };
	};

} // namespace cath::clust

#endif // _CATH_TOOLS_SOURCE_CT_CATH_CLUSTER_CATH_CATH_CLUSTER_OPTIONS_OPTIONS_BLOCK_CATH_CLUSTER_CLUSTERING_OPTIONS_BLOCK_HPP
