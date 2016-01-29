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

#ifndef DOMAIN_DEFINITION_H_INCLUDED
#define DOMAIN_DEFINITION_H_INCLUDED

#include "chopping/domain/domain.h"

namespace cath { namespace file { class pdb; } }
namespace cath { namespace opts { class data_dirs_spec; } }

namespace cath {
	namespace chop {

		/// \brief TODOCUMENT
		class domain_definition final {
		private:
			/// \brief TODOCUMENT
			chop::domain the_domain;

			/// \brief STODOCUMENT
			std::string pdb_name;

		public:
			domain_definition(const chop::domain &,
			                  const std::string &);

			const chop::domain & get_domain() const;
			const std::string & get_pdb_name() const;
		};

		file::pdb read_domain_from_pdb(const domain_definition &,
		                               const opts::data_dirs_spec &);

	}
}

#endif
