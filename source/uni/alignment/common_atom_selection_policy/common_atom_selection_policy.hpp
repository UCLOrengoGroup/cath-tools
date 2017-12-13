/// \file
/// \brief The common_atom_selection_policy class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_COMMON_ATOM_SELECTION_POLICY_COMMON_ATOM_SELECTION_POLICY_HPP
#define _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_COMMON_ATOM_SELECTION_POLICY_COMMON_ATOM_SELECTION_POLICY_HPP

#include <boost/operators.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "common/polymorphic_comparison/polymorphic_less_than_comparable.hpp"
#include "structure/structure_type_aliases.hpp"

#include <memory>
#include <string>
#include <utility>

namespace cath { namespace geom { class coord_list; } }
namespace cath { class residue; }

namespace cath {
	namespace align {

		/// \brief Provide ABC interface for policy for selecting common coordinates between residues
		///
		/// When comparing residues, it is common to want to compare equivalent coordinates. This is often
		/// done by taking the coordinates of the alpha carbon atom but this interface allows that choice
		/// to be made more abstract.
		///
		/// It's quite likely that all policies may be included in one class (or at least in one template class)
		/// but it's worth having this ABC interface for the sake of separating the abstract interface that clients
		/// use for selecting common atoms from the implementation details
		///
		/// \todo ATM, this isn't really much use because the residue class is only capable of storing
		///       the CA and CB atoms. (Though it's worth having this in place so that client code can begin to
		///       depend on this abstraction (Dependency inversion principle).
		///       Once the residue class is expanded (see todo note in the notes for residue)
		///       the range of options here should be expanded.
		///
		/// \todo Consider whether this same code could also be used for varying the atom that's used for
		///       distance measures in the SSAP score.
		class common_atom_selection_policy : private cath::common::polymorphic_less_than_comparable<common_atom_selection_policy>,
		                                     private boost::equivalent<common_atom_selection_policy,
		                                             boost::totally_ordered<common_atom_selection_policy> > {
		private:
			/// \brief Pure virtual method with which each concrete policy must define how it extracts common atoms from a pair of residues
			virtual void do_append_common_atoms_to_coord_lists(geom::coord_list_coord_list_pair &,
			                                                   const residue &,
			                                                   const residue &) const = 0;

			/// \brief Pure virtual method with which each concrete policy must descriptively name itself
			virtual std::string do_get_descriptive_name() const = 0;

			/// \brief Pure virtual method with which each concrete policy must define how to create a clone of itself
			virtual std::unique_ptr<common_atom_selection_policy> do_clone() const = 0;

			/// \brief TODOCUMENT
			virtual bool do_less_than_with_same_dynamic_type(const common_atom_selection_policy &) const = 0;

		public:
			common_atom_selection_policy() = default;
			virtual ~common_atom_selection_policy() noexcept = default;

			common_atom_selection_policy(const common_atom_selection_policy &) = default;
			common_atom_selection_policy(common_atom_selection_policy &&) noexcept = default;
			common_atom_selection_policy & operator=(const common_atom_selection_policy &) = default;
			common_atom_selection_policy & operator=(common_atom_selection_policy &&) noexcept = default;

			void append_common_atoms_to_coord_lists(geom::coord_list_coord_list_pair &,
			                                        const residue &,
			                                        const residue &) const;
			std::string get_descriptive_name() const;

			std::unique_ptr<common_atom_selection_policy> clone() const;

			bool less_than_with_same_dynamic_type(const common_atom_selection_policy &) const;
		};

		geom::coord_list_coord_list_pair select_common_atoms(const common_atom_selection_policy &,
		                                                     const residue &,
		                                                     const residue &);

		boost::ptr_vector<common_atom_selection_policy> get_all_common_atom_selection_policies();

		std::unique_ptr<common_atom_selection_policy> make_default_common_atom_selection_policy();

		bool is_default_policy(const common_atom_selection_policy &);

		/// \brief Function to make common_atom_selection_policy meet the Clonable concept (used in ptr_container)
        ///
        /// NOTE: Don't call this yourself. Call the object's clone() method instead because that returns a
        ///       smart pointer rather than the raw pointer this has to return to meet the Clonable concept.
        ///
        /// This gets the smart pointer from the clone() method and then calls release on it.
        ///
        /// \returns A raw pointer to a new copy of the common_atom_selection_policy argument, with the same dynamic type.
        ///          The caller is responsible for deleting this new object.
        inline common_atom_selection_policy * new_clone(const common_atom_selection_policy &arg_common_atom_selection_policy ///< The common_atom_selection_policy to clone
                                                        ) {
        	return arg_common_atom_selection_policy.clone().release();
        }
	} // namespace align
} // namespace cath

#endif
