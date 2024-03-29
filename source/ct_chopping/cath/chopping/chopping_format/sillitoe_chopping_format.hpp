/// \file
/// \brief The sillitoe_chopping_format class header

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

#ifndef CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_CHOPPING_FORMAT_SILLITOE_CHOPPING_FORMAT_HPP
#define CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_CHOPPING_FORMAT_SILLITOE_CHOPPING_FORMAT_HPP

#include <memory>
#include <string>
#include <string_view>
#include <utility>

#include "cath/biocore/residue_name.hpp"
#include "cath/chopping/chopping_format/chopping_format.hpp"
#include "cath/chopping/region/region.hpp"
#include "cath/common/type_aliases.hpp"

namespace cath::chop {

	/// \brief TODOCUMENT
	class sillitoe_chopping_format final : public chopping_format {
	  private:
		[[nodiscard]] std::unique_ptr<chopping_format> do_clone() const final;

		[[nodiscard]] bool do_represents_fragments() const final;

		static std::pair<::std::string_view, str_citr> parse_to_start_of_regions( const std::string & );

		[[nodiscard]] domain do_parse_domain( const std::string & ) const final;

		[[nodiscard]] std::string do_write_region( const region & ) const final;

		[[nodiscard]] std::string do_write_domain( const domain & ) const final;

	  public:
		[[nodiscard]] region parse_segment( const ::std::string_view & ) const;

		[[nodiscard]] residue_name parse_residue( const ::std::string_view & ) const;
	};

} // namespace cath::chop

#endif // CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_CHOPPING_FORMAT_SILLITOE_CHOPPING_FORMAT_HPP
