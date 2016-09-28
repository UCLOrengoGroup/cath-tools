/// \file
/// \brief The pdb_list class header

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

#ifndef PDB_LIST_H_INCLUDED
#define PDB_LIST_H_INCLUDED

#include <boost/operators.hpp>
#include <boost/range.hpp>

#include "common/path_type_aliases.h"
#include "file/file_type_aliases.h"
#include "structure/structure_type_aliases.h"

#include <iostream>
#include <vector>

namespace cath { namespace file {class pdb; } }
namespace cath { class protein_list; }

namespace cath {
	namespace file {

		/// \brief TODOCUMENT
		class pdb_list final {
		private:
			/// \brief TODOCUMENT
			pdb_vec pdbs;

		public:
			pdb_list() = default;
			pdb_list(const pdb_vec &);

			void push_back(const pdb &);
			void reserve(const size_t &);

			size_t size() const;
			bool empty() const;

			pdb & operator[](const size_t &);
			const pdb & operator[](const size_t &) const;

			// Provide iterators to make this into a range
//			using iterator = std::vector<pdb>::iterator;
			using const_iterator = std::vector<pdb>::const_iterator;
//			iterator begin();
//			iterator end();
			const_iterator begin() const;
			const_iterator end() const;
		};

		pdb_list read_pdb_files(const path_vec &);

		pdb_list make_pdb_list(const pdb_vec &);

		pdb_list pdb_list_of_backbone_complete_subset_pdbs(const pdb_list &,
		                                                   std::ostream &arg_ostream = std::cerr);

		protein_list build_protein_list_of_pdb_list(const pdb_list &);

		protein_list build_protein_list_of_pdb_list_and_names(const pdb_list &,
		                                                      const str_vec &);

		amino_acid_vec_vec get_amino_acid_lists(const pdb_list &);

		residue_name_vec_vec get_residue_names_of_first_chains__backbone_unchecked(const pdb_list &);

		residue_name_vec_vec get_backbone_complete_residue_names_of_first_chains(const pdb_list &);
	}
}

#endif
