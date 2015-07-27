/// \file
/// \brief The pdbs_acquirer class header

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

#ifndef PDBS_ACQUIRER_H_INCLUDED
#define PDBS_ACQUIRER_H_INCLUDED

#include "common/type_aliases.h"
#include "file/file_type_aliases.h"

#include <memory>
#include <utility>

namespace cath { namespace opts { class cath_superpose_options; } }
namespace cath { namespace file { class pdb_list; } }

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		class pdbs_acquirer {
		private:
			/// \brief Pure virtual method with which each concrete pdbs_acquirer must define how to create a clone of itself
			virtual std::unique_ptr<pdbs_acquirer> do_clone() const = 0;
			
			/// \brief TODOCUMENT
			virtual std::pair<cath::file::pdb_list, str_vec> do_get_pdbs_and_names(std::istream &) const = 0;

		public:
			std::unique_ptr<pdbs_acquirer> clone() const;
			virtual ~pdbs_acquirer() noexcept = default;

			file::pdb_list_str_vec_pair get_pdbs_and_names(std::istream &,
			                                               const bool &) const;
		};

		/// \brief Function to make pdbs_acquirer meet the Clonable concept (used in ptr_container)
		///
		/// NOTE: Don't call this yourself. Call the object's clone() method instead because that returns a
		///       smart pointer rather than the raw pointer this has to return to meet the Clonable concept.
		///
		/// This gets the smart pointer from the clone() method and then calls release on it.
		///
		/// \returns A raw pointer to a new copy of the pdbs_acquirer argument, with the same dynamic type.
		///          The caller is responsible for deleting this new object.
		inline pdbs_acquirer * new_clone(const pdbs_acquirer &arg_pdbs_acquirer ///< The pdbs_acquirer to clone
		                                 ) {
			return arg_pdbs_acquirer.clone().release();
		}

	}
}

#endif
