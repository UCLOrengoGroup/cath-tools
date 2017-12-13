/// \file
/// \brief The common_residue_selection_policy class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_COMMON_RESIDUE_SELECTION_POLICY_COMMON_RESIDUE_SELECTION_POLICY_HPP
#define _CATH_TOOLS_SOURCE_UNI_ALIGNMENT_COMMON_RESIDUE_SELECTION_POLICY_COMMON_RESIDUE_SELECTION_POLICY_HPP

#include <boost/ptr_container/ptr_vector.hpp>

#include "alignment/alignment.hpp"
#include "common/polymorphic_comparison/polymorphic_less_than_comparable.hpp"

#include <memory>
#include <string>
#include <vector>

namespace cath {
	namespace align {

		/// \brief TODOCUMENT
		class common_residue_selection_policy : private cath::common::polymorphic_less_than_comparable<common_residue_selection_policy>,
		                                        private boost::equivalent<common_residue_selection_policy,
		                                                boost::totally_ordered<common_residue_selection_policy> > {
		private:
			/// \brief TODOCUMENT
			virtual size_vec do_select_common_residues(const alignment &,
			                                           const std::vector<alignment::size_type> &,
			                                           const alignment::size_type &,
			                                           const alignment::size_type &) const = 0;
			
			/// \brief TODOCUMENT
			virtual std::string do_get_descriptive_name() const = 0;
			
			/// \brief Pure virtual method with which each concrete common_residue_selection_policy must define how to create a clone of itself
			virtual std::unique_ptr<common_residue_selection_policy> do_clone() const = 0;

			/// \brief TODOCUMENT
			virtual bool do_less_than_with_same_dynamic_type(const common_residue_selection_policy &) const = 0;

		public:
			common_residue_selection_policy() = default;
			virtual ~common_residue_selection_policy() noexcept = default;

			common_residue_selection_policy(const common_residue_selection_policy &) = default;
			common_residue_selection_policy(common_residue_selection_policy &&) noexcept = default;
			common_residue_selection_policy & operator=(const common_residue_selection_policy &) = default;
			common_residue_selection_policy & operator=(common_residue_selection_policy &&) noexcept = default;

			std::vector<alignment::size_type> select_common_residues(const alignment &,
			                                                         const alignment::size_type &,
			                                                         const alignment::size_type &) const;
			std::string get_descriptive_name() const;

			std::unique_ptr<common_residue_selection_policy> clone() const;

			bool less_than_with_same_dynamic_type(const common_residue_selection_policy &) const;
		};

		boost::ptr_vector<common_residue_selection_policy> get_all_common_residue_selection_policies();

		std::vector<alignment::size_type> select_common_residues_of_pair_alignment(const common_residue_selection_policy &,
		                                                                           const alignment &);

		std::unique_ptr<common_residue_selection_policy> make_default_common_residue_selection_policy();

		bool is_default_policy(const common_residue_selection_policy &);

		std::ostream & operator<<(std::ostream &,
		                          const common_residue_selection_policy &);

		/// \brief Function to make common_residue_selection_policy meet the Clonable concept (used in ptr_container)
		///
		/// NOTE: Don't call this yourself. Call the object's clone() method instead because that returns a
		///       smart pointer rather than the raw pointer this has to return to meet the Clonable concept.
		///
		/// This gets the smart pointer from the clone() method and then calls release on it.
		///
		/// \returns A raw pointer to a new copy of the common_residue_selection_policy argument, with the same dynamic type.
		///          The caller is responsible for deleting this new object.
		inline common_residue_selection_policy * new_clone(const common_residue_selection_policy &arg_common_residue_selection_policy ///< The common_residue_selection_policy to clone
		                                                ) {
			return arg_common_residue_selection_policy.clone().release();
		}
	} // namespace align
} // namespace cath

#endif
