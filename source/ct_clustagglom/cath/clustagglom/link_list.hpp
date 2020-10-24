/// \file
/// \brief The link_list class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_LINK_LIST_HPP
#define _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_LINK_LIST_HPP

#include "cath/clustagglom/link_dirn.hpp"
#include "cath/clustagglom/link.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"

namespace cath {
	namespace clust {

		/// \brief A list of links from some implied source to a bunch of destination
		class link_list final {
		private:
			/// \brief The vector of links
			link_vec links;

		public:
			/// \brief A const_iterator type alias as part of making this a range over links
			using const_iterator = link_vec::const_iterator;

			/// \brief An iterator type alias as part of making this a range over links
			///
			/// \TODOSOON ******* SHOULD THIS BE HERE? *******
			/// \TODOSOON ******* SHOULD THIS BE HERE? *******
			/// \TODOSOON ******* SHOULD THIS BE HERE? *******
			/// \TODOSOON ******* SHOULD THIS BE HERE? *******
			/// \TODOSOON ******* SHOULD THIS BE HERE? *******
			using iterator       = link_vec::iterator;

			link_list() = default;

			explicit link_list(link_vec);

			bool empty() const;
			size_t size() const;

			const link & operator[](const size_t &) const;

			template <typename... Ts>
			void emplace_back(Ts &&...);

			void clear();
			void shrink_to_fit();

			/// \TODOSOON ****** REMOVE THIS ******
			/// \TODOSOON ****** REMOVE THIS ******
			/// \TODOSOON ****** REMOVE THIS ******
			/// \TODOSOON ****** REMOVE THIS ******
			template <typename... Ts>
			auto erase(Ts &&...vars) -> decltype(auto) {
				return links.erase( std::forward<Ts>( vars )... );
			}

			// int get_best_active() {
			// }

			/// \TODOSOON ******* SHOULD THIS BE HERE? *******
			/// \TODOSOON ******* SHOULD THIS BE HERE? *******
			/// \TODOSOON ******* SHOULD THIS BE HERE? *******
			/// \TODOSOON ******* SHOULD THIS BE HERE? *******
			/// \TODOSOON ******* SHOULD THIS BE HERE? *******
			iterator begin() {
				return std::begin( links );
			}
			/// \TODONOW ******* SHOULD THIS BE HERE? *******
			/// \TODONOW ******* SHOULD THIS BE HERE? *******
			/// \TODONOW ******* SHOULD THIS BE HERE? *******
			/// \TODONOW ******* SHOULD THIS BE HERE? *******
			/// \TODONOW ******* SHOULD THIS BE HERE? *******
			iterator end() {
				return std::end( links );
			}

			const_iterator begin() const;
			const_iterator end() const;
		};

		/// \brief Ctor from a vector of links
		inline link_list::link_list(link_vec prm_links ///< The links from which to build this link_list
		                            ) : links{ prm_links } {
		}

		/// \brief Return whether this is empty
		inline bool link_list::empty() const {
			return links.empty();
		}

		/// \brief Return the number of links
		inline size_t link_list::size() const {
			return links.size();
		}

		/// \brief Get the link associated with the sequence with the specified index
		inline const link & link_list::operator[](const size_t &prm_index ///< The index of the link to access
		                                          ) const {
			return links[ prm_index ];
		}

		/// \brief Emplace a link to be constructed with the specified arguments
		template <typename... Ts>
		inline void link_list::emplace_back(Ts &&...args ///< The arguments to pass link's ctor
		                                    ) {
			links.emplace_back( std::forward<Ts>( args )... );
		}

		/// \brief Clear the list
		inline void link_list::clear() {
			links.clear();
		}

		/// \brief Shrink the memory usage to fit the current size
		inline void link_list::shrink_to_fit() {
			links.shrink_to_fit();
		}

		/// \brief Standard const begin() method, as part of making this a range over hit links
		inline auto link_list::begin() const -> const_iterator {
			return common::cbegin( links );
		}

		/// \brief Standard const end() method, as part of making this a range over hit links
		inline auto link_list::end() const -> const_iterator {
			return common::cend  ( links );
		}

		std::string link_list_string(const link_list &,
		                             const size_t &);

	} // namespace clust
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_CLUSTAGGLOM_CATH_CLUSTAGGLOM_LINK_LIST_HPP
