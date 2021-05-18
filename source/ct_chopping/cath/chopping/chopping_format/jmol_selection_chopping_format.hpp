/// \file
/// \brief The jmol_selection_chopping_format class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_CHOPPING_FORMAT_JMOL_SELECTION_CHOPPING_FORMAT_HPP
#define _CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_CHOPPING_FORMAT_JMOL_SELECTION_CHOPPING_FORMAT_HPP

#include "cath/chopping/chopping_format/chopping_format.hpp"

namespace cath {
	namespace chop {

		/// \brief TODOCUMENT
		class jmol_selection_chopping_format final : public chopping_format {
		  private:
			[[nodiscard]] std::unique_ptr<chopping_format> do_clone() const final;

			[[nodiscard]] bool do_represents_fragments() const final;

			[[nodiscard]] domain do_parse_domain( const std::string & ) const final;

			[[nodiscard]] std::string do_write_region( const region & ) const final;

			[[nodiscard]] std::string do_write_domain( const domain & ) const final;
		};

	} // namespace chop
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_CHOPPING_FORMAT_JMOL_SELECTION_CHOPPING_FORMAT_HPP
