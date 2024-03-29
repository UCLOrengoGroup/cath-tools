/// \file
/// \brief The align_regions_options_block class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_SUPERPOSITION_OPTIONS_ALIGN_REGIONS_OPTIONS_BLOCK_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_SUPERPOSITION_OPTIONS_ALIGN_REGIONS_OPTIONS_BLOCK_HPP

#include <string_view>

#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/options/options_block/options_block.hpp"

namespace cath::opts {

	/// \brief Define an options_block for options specifying how cath-resolve-hits should handle the segments
	class align_regions_options_block final : public opts::options_block {
	private:
		using super = opts::options_block;

		/// \brief The regions this options_block configures
		chop::domain_vec align_domains;

		[[nodiscard]] std::unique_ptr<opts::options_block> do_clone() const final;
		[[nodiscard]] std::string                          do_get_block_name() const final;
		void do_add_visible_options_to_description(boost::program_options::options_description &,
		                                           const size_t &) final;
		[[nodiscard]] str_opt do_invalid_string( const boost::program_options::variables_map & ) const final;
		[[nodiscard]] str_view_vec do_get_all_options_names() const final;

	  public:
		[[nodiscard]] const chop::domain_vec &get_align_domains() const;

		/// \brief The option name for the regions of the alignment (either to align or that have been aligned)
		static constexpr ::std::string_view PO_ALN_REGIONS{ "align-regions" };
	};

} // namespace cath::opts

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_SUPERPOSITION_OPTIONS_ALIGN_REGIONS_OPTIONS_BLOCK_HPP
