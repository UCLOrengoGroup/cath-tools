/// \file
/// \brief The scan_index class header

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

#ifndef _CATH_TOOLS_SOURCE_SCAN_SCAN_INDEX_H
#define _CATH_TOOLS_SOURCE_SCAN_SCAN_INDEX_H

#include <boost/numeric/conversion/cast.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/systems/information/byte.hpp>

#include "common/boost_addenda/range/back.h"
#include "common/chrono/chrono_type_aliases.h"
#include "common/debug_numeric_cast.h"
#include "scan/detail/res_pair/multi_struc_res_rep_pair_list.h"
#include "scan/detail/scan_index_store/scan_index_hash_store.h"
#include "scan/detail/scan_index_store/scan_index_store_helper.h"
#include "scan/detail/scan_multi_structure_data.h"
#include "scan/detail/scan_role.h"
#include "scan/detail/scan_type_aliases.h"
#include "scan/detail/stride/roled_scan_stride.h"
#include "scan/res_pair_keyer/res_pair_keyer.h"
#include "scan/scan_index.h"
#include "scan/scan_policy.h"
#include "structure/protein/protein_list.h"

//#include "structure/protein/sec_struc.h"
//#include "structure/protein/sec_struc_planar_angles.h"

#include <chrono> /// ***** TEMPORARY ****
#include <utility>

namespace cath { class protein; }

namespace cath {
	namespace scan {

		/// \brief TODOCUMENT
		template <typename... KPs>
		class scan_index final {
		private:
			/// \brief TODOCUMENT
			using key_t = typename res_pair_keyer<KPs...>::key_tuple_type;

			/// \brief TODOCUMENT
			std::reference_wrapper<const scan_policy<KPs...>> the_policy;

			/// \brief TODOCUMENT
			detail::scan_multi_structure_data structures_data;

			/// \brief TODOCUMENT
			durn_mem_pair structure_build_durn_and_size = make_pair( hrc_duration::zero(), 0 * boost::units::information::bytes );			

			/// \brief TODOCUMENT
			detail::scan_index_hash_store<key_t> the_store;

			/// \brief TODOCUMENT
			hrc_duration index_build_durn = hrc_duration::zero();

			void populate_index_from_structures_data();

		public:
			scan_index(const scan_policy<KPs...> &);

			/// \brief Prevent construction from a temporary scan_policy
			scan_index(const scan_policy<KPs...> &&) = delete;

			const scan_policy<KPs...> & get_scan_policy() const;

			void add_structure(const protein &);

			index_type get_num_structures() const;
			index_type get_num_residues_of_structure_of_index(const index_type &) const;

			durn_mem_pair get_structures_build_durn_and_size() const;
			durn_mem_pair get_index_build_durn_and_size() const;

			template <typename FN>
			void act_on_matches(const key_t &,
			                    const detail::scan_multi_structure_data &,
			                    const detail::multi_struc_res_rep_pair_list &,
			                    FN &) const;
		};


		/// \brief TODOCUMENT
		template <typename... KPs>
		void scan_index<KPs...>::populate_index_from_structures_data() {
			for (const auto &the_structure_data : structures_data) {
//				the_structure_data.get_from_stride();
			}
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		scan_index<KPs...>::scan_index(const scan_policy<KPs...> &arg_policy ///< TODOCUMENT
		                               ) : the_policy( arg_policy ) {
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		const scan_policy<KPs...> & scan_index<KPs...>::get_scan_policy() const {
			return the_policy;
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		void scan_index<KPs...>::add_structure(const protein &arg_protein ///< TODOCUMENT
		                                       ) {
			// BOOST_LOG_TRIVIAL( warning ) << "About to add_structure_data()";

			const auto add_structure_data_starttime = std::chrono::high_resolution_clock::now();
			add_structure_data(
				structures_data,
				arg_protein,
				detail::roled_scan_stride{ detail::scan_role::INDEX, get_scan_policy().get_scan_stride() }
			);
			structure_build_durn_and_size.first  += std::chrono::high_resolution_clock::now() - add_structure_data_starttime;
			structure_build_durn_and_size.second += common::back( structures_data ).get_info_size();

			// BOOST_LOG_TRIVIAL( warning ) << "Finished add_structure_data() - took " << durn_to_seconds_string( add_structure_data_durn );
			// BOOST_LOG_TRIVIAL( warning ) << "About to dense_add_structure_to_store()";
			const auto add_structure_to_store_starttime = std::chrono::high_resolution_clock::now();
			detail::dense_add_structure_to_store(
				the_store,
				debug_unwarned_numeric_cast<index_type>( structures_data.size() - 1 ),
				arg_protein,
				get_scan_policy(),
				detail::scan_role::INDEX
			);
			index_build_durn += std::chrono::high_resolution_clock::now() - add_structure_to_store_starttime;
			// BOOST_LOG_TRIVIAL( warning ) << "Finished dense_add_structure_to_store() - took " << durn_to_seconds_string( dense_add_structure_to_store_durn );
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		index_type scan_index<KPs...>::get_num_structures() const {
			return boost::numeric_cast<index_type>( structures_data.size() );
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		index_type scan_index<KPs...>::get_num_residues_of_structure_of_index(const index_type &arg_structure_index ///< TODOCUMENT
		                                                                      ) const {
			return structures_data[ arg_structure_index ].get_num_residues();
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		durn_mem_pair scan_index<KPs...>::get_structures_build_durn_and_size() const {
			return structure_build_durn_and_size;
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		durn_mem_pair scan_index<KPs...>::get_index_build_durn_and_size() const {
			return make_pair( index_build_durn, the_store.get_info_size() );
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		template <typename FN>
		inline void scan_index<KPs...>::act_on_matches(const key_t                                 &arg_key,                   ///< TODOCUMENT
		                                               const detail::scan_multi_structure_data     &arg_query_structures_data, ///< TODOCUMENT
		                                               const detail::multi_struc_res_rep_pair_list &arg_query_list,            ///< TODOCUMENT
		                                               FN                                          &arg_fn                     ///< TODOCUMENT
		                                               ) const {
			const auto &the_matches = the_store.find_matches( arg_key );

//			static size_t counter     = 0;
//			static size_t num_entries = 0;
//			++counter;
//			num_entries += arg_query_list.size();
//			std::cerr << counter << " " << num_entries << " comparing a query cell of " << std::right << std::setw( 5 ) << arg_query_list.size() << " entries against an index cell of " << the_matches.size() << " multi_res_pairs\n";

			if ( ! the_matches.empty() ) {
				act_on_multi_matches(
					arg_query_list,
					the_matches,
					arg_query_structures_data,
					structures_data,
					get_scan_policy().get_criteria(),
//					get_scan_policy().get_scan_stride(),
					arg_fn
				);
			}
//			else {
//				std::cerr << "Skipping due to matches cell being empty\n";
//			}
		}

		/// \brief TODOCUMENT
		template <typename... KPs>
		scan_index<KPs...> make_scan_index(const scan_policy<KPs...> &arg_policy ///< TODOCUMENT
		                                   ) {
			return { arg_policy };
		}

		/// \brief Prevent factory building from a temporary scan_policy
		template <typename... KPs>
		scan_index<KPs...> make_scan_index(const scan_policy<KPs...> &&) = delete; // Don't try to build a scan_index from a temporary scan_policy

		/// \brief TODOCUMENT
		template <typename... KPs>
		scan_index<KPs...> make_scan_index(const scan_policy<KPs...> &arg_policy,      ///< TODOCUMENT
		                                   const protein_list        &arg_protein_list ///< TODOCUMENT
		                                   ) {
			auto the_scan_index = make_scan_index( arg_policy );
			for (const protein &the_protein : arg_protein_list) {
				the_scan_index.add_structure( the_protein );
			}
			return the_scan_index;
		}

		/// \brief Prevent factory building from a temporary scan_policy
		template <typename... KPs>
		scan_index<KPs...> make_scan_index(const scan_policy<KPs...> &&,
		                                   const protein_list &) = delete; // Don't try to build a scan_index from a temporary scan_policy

	} // namespace scan
} // namespace cath

#endif
