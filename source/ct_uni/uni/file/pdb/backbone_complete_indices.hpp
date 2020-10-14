/// \file
/// \brief The backbone_complete_indices class header

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

#ifndef _CATH_TOOLS_SOURCE_UNI_FILE_PDB_BACKBONE_COMPLETE_INDICES_HPP
#define _CATH_TOOLS_SOURCE_UNI_FILE_PDB_BACKBONE_COMPLETE_INDICES_HPP

#include <boost/algorithm/cxx11/is_sorted.hpp>
#include <boost/optional.hpp>
#include <boost/range/algorithm/lower_bound.hpp>

#include "common/cpp14/cbegin_cend.hpp"
#include "common/type_aliases.hpp"

#include <cassert>

namespace cath {
	namespace file {

		/// \brief Store the indices of the backbone_complete residues
		////       so that they can be looked up quickly without having to rescan through
		///
		/// If pdb gets moved to a different namespace, then move this too
		class backbone_complete_indices final {
		private:
			/// \brief The indices of the backbone_compete residues
			///
			/// eg if residues at indices 0, 2 and 3 of five residues are backbone complete,
			/// this would contain: { 0, 2, 3 }
			size_vec indices;

		public:
			/// \brief A const_iterator type alias as part of making this a range over the indices
			using const_iterator = decltype(indices)::const_iterator;

			/// \brief Default ctor
			backbone_complete_indices() = default;
			explicit backbone_complete_indices(size_vec);

			bool empty() const;
			size_t size() const;

			const size_t & operator[](const size_t &) const;

			const_iterator begin() const;
			const_iterator end() const;
		};

		size_t get_index_of_backbone_complete_index(const backbone_complete_indices &,
		                                            const size_t &);

		size_opt get_backbone_complete_index_of_index(const backbone_complete_indices &,
		                                              const size_t &);

		/// \brief Ctor from a vector of indices
		inline backbone_complete_indices::backbone_complete_indices(size_vec prm_indices ///< The vector of indices
		                                                            ) : indices{ prm_indices } {
			// Test that the indices are *strictly* increasing
			assert( boost::algorithm::is_strictly_increasing( indices ) );
		}

		/// \brief Whether this contains any backbone_complete indices
		inline bool backbone_complete_indices::empty() const {
			return indices.empty();
		}

		/// \brief The number of backbone_complete indices that are indexed in this
		inline size_t backbone_complete_indices::size() const {
			return indices.size();
		}

		/// \brief Access the absolute index of the i-th backbone-complete residue
		inline const size_t & backbone_complete_indices::operator[](const size_t &prm_index ///< The backbone-complete index of the residue of interest
		                                                            ) const {
			return indices[ prm_index ];
		}

		/// \brief Standard const begin() method, as part of making this a range over indices
		inline auto backbone_complete_indices::begin() const -> const_iterator {
			return common::cbegin( indices );
		}

		/// \brief Standard const end() method, as part of making this a range over indices
		inline auto backbone_complete_indices::end() const -> const_iterator {
			return common::cend  ( indices );
		}

		/// \brief Get the absolute index i-th backbone-complete residue
		///
		/// \relates backbone_complete_indices
		inline size_t get_index_of_backbone_complete_index(const backbone_complete_indices &prm_bb_compl_indices, ///< The backbone_complete_indices
		                                                   const size_t                    &prm_index             ///< The index of the required residue
		                                                   ) {
			return prm_bb_compl_indices[ prm_index ];
		}

		/// \brief Get the backbone-complete index of the i-th residue
		///
		/// \relates backbone_complete_indices
		inline size_opt get_backbone_complete_index_of_index(const backbone_complete_indices &prm_bb_compl_indices, ///< The backbone_complete_indices
		                                                     const size_t                    &prm_index             ///< The index of the required residue
		                                                     ) {
			const auto lower_bound_itr = boost::range::lower_bound(
				prm_bb_compl_indices,
				prm_index
			);
			return
				( lower_bound_itr != common::cend( prm_bb_compl_indices ) && *lower_bound_itr == prm_index )
				? boost::optional<size_t>( std::distance(
					common::cbegin( prm_bb_compl_indices ),
					lower_bound_itr
				) )
				: boost::none;
		}

	} // namespace file
} // namespace cath

#endif
