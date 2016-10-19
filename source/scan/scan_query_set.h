/// \file
/// \brief The scan_query_set class header

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

#ifndef _CATH_TOOLS_SOURCE_SCAN_SCAN_QUERY_SET_H
#define _CATH_TOOLS_SOURCE_SCAN_SCAN_QUERY_SET_H

//#include <boost/log/trivial.hpp> // ***** TEMPORARY *****
#include <boost/numeric/conversion/cast.hpp>
#include <boost/throw_exception.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/systems/information/byte.hpp>

#include "common/boost_addenda/range/back.h"
#include "common/chrono/chrono_type_aliases.h"
#include "common/debug_numeric_cast.h"
#include "exception/invalid_argument_exception.h"
#include "scan/detail/scan_index_store/scan_index_store_helper.h"
#include "scan/detail/scan_index_store/scan_index_vector_store.h"
#include "scan/detail/scan_multi_structure_data.h"
#include "scan/detail/scan_role.h"
#include "scan/detail/scan_type_aliases.h"
#include "scan/detail/stride/roled_scan_stride.h"
#include "scan/res_pair_keyer/res_pair_keyer.h"
#include "scan/scan_index.h"
#include "scan/scan_policy.h"

#include <cassert>
#include <chrono>
#include <utility>

namespace cath { class protein; }
namespace cath { class protein_list; }
namespace cath { namespace scan { namespace detail { class multi_struc_res_rep_pair; } } }

namespace cath {
	namespace scan {

		/// \brief TODOCUMENT
		template <typename... KPs>
		class scan_query_set final {
		private:
			/// \brief TODOCUMENT
			using key_t = typename res_pair_keyer<KPs...>::key_tuple_type;

//			/// \brief TODOCUMENT
//			using store_t = detail::scan_index_vector_store<key_t>;

//			/// \brief TODOCUMENT
//			using store_value_type = common::range_value_t<store_t>;

			/// \brief TODOCUMENT
			std::reference_wrapper<const scan_policy<KPs...>> the_policy;

			/// \brief TODOCUMENT
			detail::scan_multi_structure_data structures_data;

			/// \brief TODOCUMENT
			durn_mem_pair structure_build_durn_and_size = make_pair( hrc_duration::zero(), 0 * boost::units::information::bytes );			

			/// \brief TODOCUMENT
			detail::scan_index_vector_store<key_t> the_store;

			/// \brief TODOCUMENT
			hrc_duration index_build_durn = hrc_duration::zero();

			void add_entry(const detail::multi_struc_res_rep_pair &);

		public:
			scan_query_set(const scan_policy<KPs...> &);

			/// \brief Prevent construction from a temporary scan_policy
			scan_query_set(const scan_policy<KPs...> &&) = delete;

			const scan_policy<KPs...> & get_scan_policy() const;

			void add_structure(const protein &);

			index_type get_num_structures() const;
			index_type get_num_residues_of_structure_of_index(const index_type &) const;

			durn_mem_pair get_structures_build_durn_and_size() const;
			durn_mem_pair get_index_build_durn_and_size() const;

			template <typename FN>
			hrc_duration do_magic(const scan_index<KPs...> &,
			                      FN &) const;

			/// \brief TODOCUMENT
			template <typename FN>
			void act_on_matches(FN &) const;
		};

		/// \brief TODOCUMENT
		template <typename... KPs>
		scan_query_set<KPs...>::scan_query_set(const scan_policy<KPs...> &arg_scan_policy ///< TODOCUMENT
		                                       ) : the_policy( arg_scan_policy ) {
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		const scan_policy<KPs...> & scan_query_set<KPs...>::get_scan_policy() const {
			return the_policy;
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		void scan_query_set<KPs...>::add_structure(const protein &arg_protein ///< TODOCUMENT
		                                           ) {
			// BOOST_LOG_TRIVIAL( warning ) << "About to add structure to structures_data...";

			const auto add_structure_data_starttime = std::chrono::high_resolution_clock::now();
			add_structure_data(
				structures_data,
				arg_protein,
				detail::roled_scan_stride{ detail::scan_role::QUERY, get_scan_policy().get_scan_stride() }
			);

			structure_build_durn_and_size.first  += std::chrono::high_resolution_clock::now() - add_structure_data_starttime;
			structure_build_durn_and_size.second += common::back( structures_data ).get_info_size();

			// BOOST_LOG_TRIVIAL( warning ) << "Added structure to structures_data";
			// BOOST_LOG_TRIVIAL( warning ) << "About to add structure to store...";

			const auto add_structure_to_store_starttime = std::chrono::high_resolution_clock::now();
			detail::add_structure_to_store(
				the_store,
				debug_unwarned_numeric_cast<index_type>( structures_data.size() - 1 ),
				arg_protein,
				get_scan_policy(),
				detail::scan_role::QUERY
			);

			index_build_durn += std::chrono::high_resolution_clock::now() - add_structure_to_store_starttime;
			// BOOST_LOG_TRIVIAL( warning ) << "Added structure to store";
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		index_type scan_query_set<KPs...>::get_num_structures() const {
			return boost::numeric_cast<index_type>( structures_data.size() );
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		index_type scan_query_set<KPs...>::get_num_residues_of_structure_of_index(const index_type &arg_structure_index ///< TODOCUMENT
		                                                                          ) const {
			return boost::numeric_cast<index_type>( structures_data[ arg_structure_index ].get_num_residues() );
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		durn_mem_pair scan_query_set<KPs...>::get_structures_build_durn_and_size() const {
			return structure_build_durn_and_size;
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		durn_mem_pair scan_query_set<KPs...>::get_index_build_durn_and_size() const {
			return make_pair( index_build_durn, the_store.get_info_size() );
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		template <typename FN>
		hrc_duration scan_query_set<KPs...>::do_magic(const scan_index<KPs...> &arg_scan_index, ///< TODOCUMENT
		                                              FN                       &arg_fn          ///< TODOCUMENT
		                                              ) const {
			if ( &( arg_scan_index.get_scan_policy() ) != & ( get_scan_policy() ) ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Unable to scan query_set against index constructed with different policy"));
			}

			const auto scan_starttime = std::chrono::high_resolution_clock::now();
			for (const auto &x : the_store) {
				const auto &key           = x.first;
				const auto &res_pair_list = x.second;
//				BOOST_LOG_TRIVIAL( warning ) << "Searching query key " << detail::output_key( key ) << " list of size " << res_pair_list.size();
				assert( ! res_pair_list.empty() );
				arg_scan_index.act_on_matches( key, structures_data, res_pair_list, arg_fn );
			}
			return std::chrono::high_resolution_clock::now() - scan_starttime;
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		scan_query_set<KPs...> make_scan_query_set(const scan_policy<KPs...> &arg_policy ///< TODOCUMENT
		                                           ) {
			return { arg_policy };
		}

		/// \brief Prevent factory building from a temporary scan_policy
		template <typename... KPs>
		scan_query_set<KPs...> make_scan_query_set(const scan_policy<KPs...> &&) = delete; // Don't try to build a scan_query_set from a temporary scan_policy

		/// \brief TODOCUMENT
		template <typename... KPs>
		scan_query_set<KPs...> make_scan_query_set(const scan_policy<KPs...> &arg_scan_policy, ///< TODOCUMENT
		                                           const protein_list        &arg_protein_list ///< TODOCUMENT
		                                           ) {
			auto new_query_set = make_scan_query_set( arg_scan_policy );
			for (const protein &the_protein : arg_protein_list) {
				new_query_set.add_structure( the_protein );
			}
			return new_query_set;
		}

		/// \brief Prevent factory building from a temporary scan_policy
		template <typename... KPs>
		scan_query_set<KPs...> make_scan_query_set(const scan_policy<KPs...> &&,
		                                           const protein_list &) = delete; // Don't try to build a scan_query_set from a temporary scan_policy

	} // namespace scan
} // namespace cath

#endif
