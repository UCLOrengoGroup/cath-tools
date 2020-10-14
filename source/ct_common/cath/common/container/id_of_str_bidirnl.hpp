/// \file
/// \brief The id_of_str_bidirnl header

/// \copyright
/// Tony Lewis's Common C++ Library Code (here imported into the CATH Tools project and then tweaked, eg namespaced in cath)
/// Copyright (C) 2007, Tony Lewis
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

#ifndef _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_CONTAINER_ID_OF_STR_BIDIRNL_HPP
#define _CATH_TOOLS_SOURCE_SRC_COMMON_COMMON_CONTAINER_ID_OF_STR_BIDIRNL_HPP

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/functional.hpp>
#include <boost/optional.hpp>
#include <boost/range/empty.hpp>

#include "cath/common/boost_addenda/range/front.hpp"
#include "cath/common/boost_addenda/range/max_proj_element.hpp"
#include "cath/common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "cath/common/container/id_of_string_view.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"
#include "cath/common/type_aliases.hpp"

#include <string>

namespace cath {
	namespace common {

		/// \brief Store cluster names and associate an (non-negative integral type) id with each
		///
		/// This isn't thread-safe (except in the standard way)
		///
		/// The IDs are guaranteed to be stable (in-between calls to clear()) and
		/// sensible for indexing into vectors without wasting much space.
		class id_of_str_bidirnl final {
		private:
			/// \brief A map from the names to the corresponding ID
			id_of_string_view ids_by_name;

			/// \brief A deque of the names where each ID is the index of the corresponding string
			str_deq           names_by_id;

		public:
			/// \brief A const_iterator type alias as part of making this a range over the names
			using const_iterator = str_deq::const_iterator;

			/// \brief Default ctor
			id_of_str_bidirnl() = default;
			inline size_t add_name(const boost::string_ref &);
			inline size_t add_name(const std::string &);
			inline size_t add_name(std::string &&);
			inline const std::string & get_name_of_id(const size_t &) const;
			inline size_t get_id_of_name(const std::string &) const;
			inline size_t get_id_of_name(const boost::string_ref &) const;
			inline bool empty() const;
			inline size_t size() const;
			inline id_of_str_bidirnl & clear();
			inline const_iterator begin() const;
			inline const_iterator end() const;
		};

		/// \brief Function template for doing the name-lookup to pickup to_string() functions
		///        in the arguments' namespace(s) via ADL or fall back to std::to_string() otherwise
		///
		/// \TODO Move this somewhere else and make it a Niebler-esque function-object
		template <typename... Ts>
		auto to_string_hook(Ts &&...args ///< The arguments to pass to std::to_string()
		                    ) -> decltype( auto ) {
			using std::to_string;
			return to_string( std::forward<Ts>( args )... );
		}

		/// \brief Add the specified name and return its ID
		///
		/// Can be used if the name already exists
		inline size_t id_of_str_bidirnl::add_name(const boost::string_ref &prm_name ///< The name to add
		                                          ) {
			if ( ids_by_name.contains( prm_name ) ) {
				return *ids_by_name[ prm_name ];
			}
			names_by_id.push_back( prm_name.to_string() );
			const size_t &id = ids_by_name.emplace( boost::string_ref{ names_by_id.back() } ).second;
			if ( id + 1 != names_by_id.size() ) {
				BOOST_THROW_EXCEPTION(out_of_range_exception(
					"id_of_str_bidirnl tried to add a name with ID "
					+ to_string_hook( id                 )
					+ " when only "
					+ to_string_hook( names_by_id.size() )
					+ " names have been stored."
				));
			}
			return id;
		}

		/// \brief Add the specified name and return its ID
		///
		/// Can be used if the name already exists
		inline size_t id_of_str_bidirnl::add_name(const std::string &prm_name ///< The name to add
		                                          ) {
			return add_name( boost::string_ref{ prm_name } );
		}

		/// \brief Add the specified name and return its ID
		///
		/// Can be used if the name already exists
		inline size_t id_of_str_bidirnl::add_name(std::string &&prm_name ///< The name to add
		                                          ) {
			if ( ids_by_name.contains( prm_name ) ) {
				return *ids_by_name[ prm_name ];
			}
			names_by_id.push_back( std::move( prm_name ) );
			const size_t &id = ids_by_name.emplace( boost::string_ref{ names_by_id.back() } ).second;
			if ( id + 1 != names_by_id.size() ) {
				BOOST_THROW_EXCEPTION(out_of_range_exception(
					"id_of_str_bidirnl tried to add a name with ID "
					+ to_string_hook( id                 )
					+ " when only "
					+ to_string_hook( names_by_id.size() )
					+ " names have been stored."
				));
			}
			return id;
		}

		/// \brief Get the name associated with the specified ID
		///
		/// \pre The ID must be a valid ID else this triggers undefined behaviour
		///      (or something less nasty on a range-checked debug build)
		inline const std::string & id_of_str_bidirnl::get_name_of_id(const size_t &prm_id ///< The ID to query
		                                                             ) const {
			return names_by_id[ prm_id ];
		}

		/// \brief Get the ID associated with the specified name
		inline auto id_of_str_bidirnl::get_id_of_name(const boost::string_ref &prm_name ///< The name to query
		                                              ) const -> size_t {
			return *ids_by_name[ prm_name ];
		}

		/// \brief Get the ID associated with the specified name
		inline auto id_of_str_bidirnl::get_id_of_name(const std::string &prm_name ///< The name to query
		                                              ) const -> size_t {
			return *ids_by_name[ prm_name ];
		}

		/// \brief Return whether the id_of_str_bidirnl is empty
		inline bool id_of_str_bidirnl::empty() const {
			return names_by_id.empty();
		}

		/// \brief The number of names currently stored in the id_of_str_bidirnl
		inline size_t id_of_str_bidirnl::size() const {
			return names_by_id.size();
		}

		/// \brief Clear the id_of_str_bidirnl
		///
		/// Invalidates all IDs and iterators
		inline id_of_str_bidirnl & id_of_str_bidirnl::clear() {
			names_by_id.clear();
			ids_by_name.clear();
			return *this;
		}

		/// \brief Standard const begin() operator to provide range access
		inline auto id_of_str_bidirnl::begin() const -> const_iterator {
			return common::cbegin( names_by_id );
		}

		/// \brief Standard const end() operator to provide range access
		inline auto id_of_str_bidirnl::end() const -> const_iterator {
			return common::cend  ( names_by_id );
		}

		/// \brief Get the largest number if the specified strings are all numbers, or return none otherwise
		///
		/// This rejects strings unless they're: \-?\d+
		/// This is stricter than stol, in the sense that stol will convert "-1.2e2" to 1
		///
		/// \todo Consider just doing a regex instead
		template <typename StrRng>
		boost::optional<ptrdiff_t> largest_number_if_names_all_numeric_integers(const StrRng &prm_strings ///< The range of strings to query
		                                                                        ) {
			using const_str_ref = range_const_reference_t<StrRng>;
			const auto str_is_not_int_like_fn = [] (const_str_ref x) {
				return (
					boost::empty( x )
					||
					( ! std::isdigit( front( x ) ) && ( front( x ) != '-' ) )
					||
					( ( front( x ) == '-' ) && boost::size( x ) == 1 )
					||
					boost::algorithm::any_of(
						next( common::cbegin( x ) ),
						common::cend( x ),
						[] (const auto &y) { return ! boost::algorithm::is_digit()( y ); }
					)
				);
			};
			if ( boost::empty( prm_strings ) || boost::algorithm::any_of( prm_strings, str_is_not_int_like_fn ) ) {
				return boost::none;
			}
			return max_proj(
				prm_strings,
				std::less<>{},
				[] (const_str_ref x) { return std::stol( x ); }
			);
		}

	} // namespace common
} // namespace cath

#endif
