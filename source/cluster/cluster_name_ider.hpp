/// \file
/// \brief The cluster_name_ider class header

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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_CLUSTER_NAME_IDER_H
#define _CATH_TOOLS_SOURCE_CLUSTER_CLUSTER_NAME_IDER_H

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/functional.hpp>
#include <boost/optional.hpp>
#include <boost/range/empty.hpp>

#include "cluster/cluster_type_aliases.hpp"
#include "common/boost_addenda/range/front.hpp"
#include "common/boost_addenda/range/max_proj_element.hpp"
#include "common/boost_addenda/range/range_concept_type_aliases.hpp"
#include "common/container/id_of_string_view.hpp"
#include "common/type_aliases.hpp"
#include "exception/out_of_range_exception.hpp"

#include <string>

namespace cath {
	namespace clust {

		/// \brief Store cluster names and associate an (non-negative integral type) id with each
		///
		/// This isn't thread-safe (except in the standard way)
		///
		/// The IDs are guaranteed to be stable (in-between calls to clear()) and
		/// sensible for indexing into vectors without wasting much space.
		class cluster_name_ider final {
		private:
			/// \brief A map from the names to the corresponding ID
			common::id_of_string_view ids_by_name;

			/// \brief A deque of the names where each ID is the index of the corresponding string
			str_deq                   names_by_id;

		public:
			/// \brief A const_iterator type alias as part of making this a range over the names
			using const_iterator = str_deq::const_iterator;

			/// \brief Default ctor
			cluster_name_ider() = default;
			inline cluster_id_t add_name(const boost::string_ref &);
			inline cluster_id_t add_name(const std::string &);
			inline cluster_id_t add_name(std::string &&);
			inline const std::string & get_name_of_id(const size_t &) const;
			inline cluster_id_t get_id_of_name(const std::string &) const;
			inline cluster_id_t get_id_of_name(const boost::string_ref &) const;
			inline bool empty() const;
			inline size_t size() const;
			inline cluster_name_ider & clear();
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
		inline cluster_id_t cluster_name_ider::add_name(const boost::string_ref &arg_name ///< The name to add
		                                                ) {
			if ( ids_by_name.contains( arg_name ) ) {
				return ids_by_name[ arg_name ];
			}
			names_by_id.push_back( arg_name.to_string() );
			const cluster_id_t &id = ids_by_name.emplace( boost::string_ref{ names_by_id.back() } ).second;
			if ( id + 1 != names_by_id.size() ) {
				BOOST_THROW_EXCEPTION(common::out_of_range_exception(
					"cluster_name_ider tried to add a name with ID "
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
		inline cluster_id_t cluster_name_ider::add_name(const std::string &arg_name ///< The name to add
		                                                ) {
			return add_name( boost::string_ref{ arg_name } );
		}

		/// \brief Add the specified name and return its ID
		///
		/// Can be used if the name already exists
		inline cluster_id_t cluster_name_ider::add_name(std::string &&arg_name ///< The name to add
		                                                ) {
			if ( ids_by_name.contains( arg_name ) ) {
				return ids_by_name[ arg_name ];
			}
			names_by_id.push_back( std::move( arg_name ) );
			const cluster_id_t &id = ids_by_name.emplace( boost::string_ref{ names_by_id.back() } ).second;
			if ( id + 1 != names_by_id.size() ) {
				BOOST_THROW_EXCEPTION(common::out_of_range_exception(
					"cluster_name_ider tried to add a name with ID "
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
		inline const std::string & cluster_name_ider::get_name_of_id(const size_t &arg_id ///< The ID to query
		                                                             ) const {
			return names_by_id[ arg_id ];
		}

		/// \brief Get the ID associated with the specified name
		inline auto cluster_name_ider::get_id_of_name(const boost::string_ref &arg_name ///< The name to query
		                                              ) const -> cluster_id_t {
			return ids_by_name[ arg_name ];
		}

		/// \brief Get the ID associated with the specified name
		inline auto cluster_name_ider::get_id_of_name(const std::string &arg_name ///< The name to query
		                                              ) const -> cluster_id_t {
			return ids_by_name[ arg_name ];
		}

		/// \brief Return whether the cluster_name_ider is empty
		inline bool cluster_name_ider::empty() const {
			return names_by_id.empty();
		}

		/// \brief The number of names currently stored in the cluster_name_ider
		inline size_t cluster_name_ider::size() const {
			return names_by_id.size();
		}

		/// \brief Clear the cluster_name_ider
		///
		/// Invalidates all IDs and iterators
		inline cluster_name_ider & cluster_name_ider::clear() {
			names_by_id.clear();
			ids_by_name.clear();
			return *this;
		}

		/// \brief Standard const begin() operator to provide range access
		inline auto cluster_name_ider::begin() const -> const_iterator {
			return common::cbegin( names_by_id );
		}

		/// \brief Standard const end() operator to provide range access
		inline auto cluster_name_ider::end() const -> const_iterator {
			return common::cend  ( names_by_id );
		}

		/// \brief Get the largest number if the specified strings are all numbers, or return none otherwise
		///
		/// This rejects strings unless they're: \-?\d+
		/// This is stricter than stol, in the sense that stol will convert "-1.2e2" to 1
		///
		/// \todo Consider just doing a regex instead
		template <typename StrRng>
		boost::optional<ptrdiff_t> largest_number_if_names_all_numeric_integers(const StrRng &arg_strings ///< The range of strings to query
		                                                                        ) {
			using const_str_ref = common::range_const_reference_t<StrRng>;
			const auto str_is_not_int_like_fn = [] (const_str_ref x) {
				return (
					boost::empty( x )
					||
					( ! std::isdigit( cath::common::front( x ) ) && ( cath::common::front( x ) != '-' ) )
					||
					( ( cath::common::front( x ) == '-' ) && boost::size( x ) == 1 )
					||
					boost::algorithm::any_of(
						next( common::cbegin( x ) ),
						common::cend( x ),
						[] (const auto &x) { return ! boost::algorithm::is_digit()( x ); }
					)
				);
			};
			if ( boost::empty( arg_strings ) || boost::algorithm::any_of( arg_strings, str_is_not_int_like_fn ) ) {
				return boost::none;
			}
			return cath::common::max_proj(
				arg_strings,
				std::less<>{},
				[] (const_str_ref x) { return std::stol( x ); }
			);
		}

	} // namespace clust
} // namespace cath

#endif
