/// \file
/// \brief The cath_assign_domains_options_block class header

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

#ifndef CATH_TOOLS_SOURCE_CT_CATH_ASSIGN_DOMAINS_CATH_CATH_ASSIGN_DOMAINS_OPTIONS_CATH_ASSIGN_DOMAINS_OPTIONS_BLOCK_HPP
#define CATH_TOOLS_SOURCE_CT_CATH_ASSIGN_DOMAINS_CATH_CATH_ASSIGN_DOMAINS_OPTIONS_CATH_ASSIGN_DOMAINS_OPTIONS_BLOCK_HPP

#include <cstddef>
#include <filesystem>
#include <memory>
#include <string>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>

#include "cath/common/type_aliases.hpp"
#include "cath/options/options_block/options_block.hpp"

// clang-format off
namespace cath::opts { class pdbs_acquirer; }
// clang-format on

namespace cath::opts {

	/// \brief Options block for the main cath-assign-domains options
	class cath_assign_domains_options_block final : public cath::opts::options_block {
	private:
		/// \brief The SVM-light RBF model file
		::std::filesystem::path rbf_svm_file;

		/// \brief The file containing the list of PRC/SSAP data files
		::std::filesystem::path data_data_file;

		/// \brief The file containing superfamily of domain
		::std::filesystem::path sf_of_dom_file;

		/// \brief The list of CATH nodes forbidden for assignment
		str_vec forbidden_nodes;

		[[nodiscard]] std::unique_ptr<options_block> do_clone() const final;
		[[nodiscard]] std::string                    do_get_block_name() const final;
		void do_add_visible_options_to_description(boost::program_options::options_description &,
		                                           const size_t &) final;
		[[nodiscard]] str_opt do_invalid_string( const boost::program_options::variables_map & ) const final;
		[[nodiscard]] str_view_vec do_get_all_options_names() const final;

	  public:
		[[nodiscard]] const ::std::filesystem::path &get_rbf_svm_file() const;
		[[nodiscard]] const ::std::filesystem::path &get_data_data_file() const;
		[[nodiscard]] const ::std::filesystem::path &get_sf_of_dom_file() const;
		[[nodiscard]] const str_vec &                get_forbidden_nodes() const;
	};

} // namespace cath::opts

#endif // CATH_TOOLS_SOURCE_CT_CATH_ASSIGN_DOMAINS_CATH_CATH_ASSIGN_DOMAINS_OPTIONS_CATH_ASSIGN_DOMAINS_OPTIONS_BLOCK_HPP
