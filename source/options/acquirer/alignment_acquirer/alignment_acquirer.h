/// \file
/// \brief The alignment_acquirer class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef ALIGNMENT_ACQUIRER_H_INCLUDED
#define ALIGNMENT_ACQUIRER_H_INCLUDED

#include <boost/tuple/tuple.hpp>

#include "common/type_aliases.h"
#include "file/pdb/pdb_list.h"

#include <memory>
#include <utility>

namespace cath { namespace align { class alignment; } }
namespace cath { namespace opts {class pdbs_acquirer; } }
namespace cath { namespace sup { class superpose_orderer; } }

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		class alignment_acquirer {
		private:
			/// \brief Pure virtual method with which each concrete alignment_acquirer must define how to create a clone of itself
			virtual std::unique_ptr<alignment_acquirer> do_clone() const = 0;
			
			/// \brief TODOCUMENT
			virtual std::pair<align::alignment, sup::superpose_orderer> do_get_alignment_and_orderer(const file::pdb_list &) const = 0;

		protected:
			/// \brief The minimum number of residues that are required in "residue name" aligning
			///
			/// \todo Decide: Is it correct that this is currently stored here rather than in the residue_name_alignment_acquirer
			///       so that it can be used in the error message when trying to construct the tree?
			///       Would it be better for the residue_name_alignment_acquirer to identify if it has failed to make a tree
			///       for this reason and then let the alignment_acquirer give a more generic message?
			static const size_t MIN_NUM_COMMON_RESIDUES_TO_SUPERPOSE_PAIR;

		public:
			std::unique_ptr<alignment_acquirer> clone() const;
			virtual ~alignment_acquirer() noexcept = default;

			std::pair<align::alignment, size_size_pair_vec> get_alignment_and_spanning_tree(const file::pdb_list &) const;
		};

		/// \brief Function to make alignment_acquirer meet the Clonable concept (used in ptr_container)
		///
		/// NOTE: Don't call this yourself. Call the object's clone() method instead because that returns a
		///       smart pointer rather than the raw pointer this has to return to meet the Clonable concept.
		///
		/// This gets the smart pointer from the clone() method and then calls release on it.
		///
		/// \returns A raw pointer to a new copy of the alignment_acquirer argument, with the same dynamic type.
		///          The caller is responsible for deleting this new object.
		inline alignment_acquirer * new_clone(const alignment_acquirer &arg_alignment_acquirer ///< The alignment_acquirer to clone
		                                      ) {
			return arg_alignment_acquirer.clone().release();
		}

	}
}

#endif
