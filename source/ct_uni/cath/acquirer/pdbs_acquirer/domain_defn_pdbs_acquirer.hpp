/// \file
/// \brief The domain_defn_pdbs_acquirer class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_PDBS_ACQUIRER_DOMAIN_DEFN_PDBS_ACQUIRER_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_PDBS_ACQUIRER_DOMAIN_DEFN_PDBS_ACQUIRER_HPP

#include <filesystem>

#include "cath/acquirer/pdbs_acquirer/pdbs_acquirer.hpp"

namespace cath::opts {

	/// \brief TODOCUMENT
	class domain_defn_pdbs_acquirer final : public pdbs_acquirer {
	private:
		/// \brief TODOCUMENT
		::std::filesystem::path domain_defn_file;

		[[nodiscard]] std::unique_ptr<pdbs_acquirer> do_clone() const final;
		file::pdb_list_name_set_list_pair do_get_pdbs_and_names(std::istream &) const final;

	  public:
		explicit domain_defn_pdbs_acquirer( ::std::filesystem::path );
	};

} // namespace cath::opts

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_ACQUIRER_PDBS_ACQUIRER_DOMAIN_DEFN_PDBS_ACQUIRER_HPP
