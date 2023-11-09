/// \file
/// \brief The domain_definition class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_DOMAIN_DOMAIN_DEFINITION_HPP
#define _CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_DOMAIN_DOMAIN_DEFINITION_HPP

#include <string>

#include "cath/chopping/domain/domain.hpp"

// clang-format off
namespace cath::file { class pdb; }
namespace cath::opts { class data_dirs_spec; }
// clang-format on

namespace cath::chop {

	/// \brief TODOCUMENT
	class domain_definition final {
	private:
		/// \brief TODOCUMENT
		chop::domain the_domain;

		/// \brief STODOCUMENT
		std::string pdb_name;

	public:
		domain_definition(chop::domain,
		                  std::string);

		[[nodiscard]] const chop::domain &get_domain() const;
		[[nodiscard]] const std::string & get_pdb_name() const;
	};

} // namespace cath::chop

#endif // _CATH_TOOLS_SOURCE_CT_CHOPPING_CATH_CHOPPING_DOMAIN_DOMAIN_DEFINITION_HPP
