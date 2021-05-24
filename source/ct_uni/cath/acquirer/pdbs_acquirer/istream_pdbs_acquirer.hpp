/// \file
/// \brief The istream_pdbs_acquirer class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_PDBS_ACQUIRER_ISTREAM_PDBS_ACQUIRER_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_PDBS_ACQUIRER_ISTREAM_PDBS_ACQUIRER_HPP

#include "cath/acquirer/pdbs_acquirer/pdbs_acquirer.hpp"

#include <string>

namespace cath::opts {

	/// \brief TODOCUMENT
	class istream_pdbs_acquirer final : public pdbs_acquirer {
	  private:
		[[nodiscard]] std::unique_ptr<pdbs_acquirer> do_clone() const final;
		file::pdb_list_name_set_list_pair            do_get_pdbs_and_names( std::istream & ) const final;
	};

} // namespace cath::opts

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_PDBS_ACQUIRER_ISTREAM_PDBS_ACQUIRER_HPP
