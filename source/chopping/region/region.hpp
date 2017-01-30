/// \file
/// \brief The region class header

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

#ifndef _CATH_TOOLS_SOURCE_CHOPPING_REGION_REGION_H
#define _CATH_TOOLS_SOURCE_CHOPPING_REGION_REGION_H

#include "biocore/chain_label.hpp"
#include "chopping/region/region_comparison.hpp"
#include "chopping/residue_location/residue_location.hpp"

#include <cstddef>

namespace cath {
	namespace chop {

		/// \brief TODOCUMENT
		///
		/// Note that the region includes both the start and stop residues
		/// (ie it's a fully-closed interval unlike the standard C++ ranges)
		///
		/// Invariants:
		///  * If the residues have indices, the start index must not exceed the stop index
		///  * The start_residue and stop_residue must have matching residue_locating
		///  * If residues are specified with names then a valid chain label must also be specified
		///  * If residues are not specified with names then a chain label may not be specified
		///    (though this may be dropped if there's a justifying context)
		///
		/// The first two invariants are checked in sanity_check().
		///
		/// The last three invariants are maintained through limited ctors and no modifiers.
		class region final {
		private:
			/// \brief TODOCUMENT
			chain_label_opt the_chain_label;

			/// \brief TODOCUMENT
			residue_location start_residue;

			/// \brief TODOCUMENT
			residue_location stop_residue;

			void sanity_check() const;

		public:
			region(const chain_label &,
			       const residue_name &,
			       const residue_name &);
			region(const chain_label &,
			       const residue_name &,
			       const size_t &,
			       const residue_name &,
			       const size_t &);
			region(const size_t &,
			       const size_t &);

			const chain_label_opt & get_opt_chain_label() const;

			const residue_location & get_start_residue() const;
			const residue_location & get_stop_residue() const;
		};

		bool has_chain_label(const region &);
		chain_label get_chain_label(const region &);

		const residue_name_opt & get_opt_start_name(const region &);
		const residue_name_opt & get_opt_stop_name(const region &);
		const size_opt & get_opt_start_index(const region &);
		const size_opt & get_opt_stop_index(const region &);

		bool has_names(const region &);
		bool has_indices(const region &);

		void check_has_names(const region &);
		void check_has_indices(const region &);

		residue_name_opt get_start_name(const region &);
		residue_name_opt get_stop_name(const region &);
		size_t get_start_index(const region &);
		size_t get_stop_index(const region &);

		residue_locating get_residue_locating(const region &);

		size_t get_length(const region &);

		region_comparison compare_locations(const region &,
		                                    const region &);

	} // namespace chop
} // namespace cath

#endif
