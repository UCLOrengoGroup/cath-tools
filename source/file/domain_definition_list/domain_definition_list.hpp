/// \file
/// \brief The domain_definition_list class header

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

#ifndef _CATH_TOOLS_SOURCE_FILE_DOMAIN_DEFINITION_LIST_DOMAIN_DEFINITION_LIST_H
#define _CATH_TOOLS_SOURCE_FILE_DOMAIN_DEFINITION_LIST_DOMAIN_DEFINITION_LIST_H

#include <boost/filesystem/path.hpp>

#include "chopping/chopping_type_aliases.hpp"
#include "file/file_type_aliases.hpp"

namespace cath { namespace opts { class data_dirs_spec; } }

namespace cath {
	namespace file {

		/// \brief TODOCUMENT
		class domain_definition_list final {
		private:
			/// \brief TODOCUMENT
			chop::domain_definition_vec domain_definitions;

		public:
			using const_iterator = chop::domain_definition_vec::const_iterator;

			explicit domain_definition_list(const chop::domain_definition_vec &);

			size_t size() const;

			const_iterator begin() const;
			const_iterator end() const;
		};

		domain_definition_list parse_domain_definition_file(const boost::filesystem::path &);
		domain_definition_list parse_domain_definition_file(std::istream &);

		pdb_list_str_vec_pair read_domains_from_pdbs(const domain_definition_list &,
		                                             const opts::data_dirs_spec &);

	} // namespace file
} // namespace cath

#endif
