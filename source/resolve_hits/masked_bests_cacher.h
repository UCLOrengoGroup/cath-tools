/// \file
/// \brief The masked_bests_cacher class header

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

#ifndef MASKED_BESTS_CACHER_H_INCLUDED
#define MASKED_BESTS_CACHER_H_INCLUDED

#include <boost/range/sub_range.hpp>

#include "resolve_hits/masked_bests_cache.h"
// #include "resolve_hits/resolve_hits_type_aliases.h"

// #include <functional>
//
// namespace cath { namespace rslv { class best_scan_arches; } }
namespace cath { namespace rslv { class hit_list; } }

namespace cath {
	namespace rslv {

		namespace detail {

			/// \brief TODOCUMENT
			class masked_bests_cacher final {
			private:
				/// \brief TODOCUMENT
				using const_iterator = res_arrow_vec::const_iterator;

				/// \brief TODOCUMENT
				std::reference_wrapper<masked_bests_cache> cache_ref;

				/// \brief TODOCUMENT
				std::reference_wrapper<const hit_vec> masks_ref;

				/// \brief TODOCUMENT
				const res_arrow_vec arrows_to_store;

				/// \brief TODOCUMENT
				const_iterator current_itr;

				void advance_to_itr_with_best_so_far(const const_iterator &,
				                                     const scored_arch_proxy &);

			public:
				masked_bests_cacher(masked_bests_cache &,
				                    const hit_vec &,
				                    const res_arrow_vec &);

				void advance_to_pos_with_best_so_far(const res_arrow &,
				                                     const scored_arch_proxy &);

				void advance_to_end_with_best_so_far(const scored_arch_proxy &);
			};

			/// \brief TODOCUMENT
			inline void masked_bests_cacher::advance_to_itr_with_best_so_far(const const_iterator    &arg_itr,        ///< TODOCUMENT
			                                                                 const scored_arch_proxy &arg_best_so_far ///< TODOCUMENT
			                                                                 ) {
				for (const res_arrow &the_arrow : boost::sub_range<const res_arrow_vec>( current_itr, arg_itr ) ) {
					store_best_for_masks_up_to_arrow(
						cache_ref.get(),
						arg_best_so_far,
						masks_ref.get(),
						the_arrow
					);
				}
				current_itr = arg_itr;
			}

			/// \brief TODOCUMENT
			///
			/// \pre `is_sorted( arg_arrows )`
			inline masked_bests_cacher::masked_bests_cacher(masked_bests_cache  &arg_masked_bests_cache, ///< TODOCUMENT
			                                                const hit_vec       &arg_masks,              ///< TODOCUMENT
			                                                const res_arrow_vec &arg_arrows              ///< TODOCUMENT
			                                                ) : cache_ref       ( arg_masked_bests_cache            ),
			                                                    masks_ref       ( arg_masks                         ),
			                                                    arrows_to_store ( arg_arrows                        ),
			                                                    current_itr     ( common::cbegin( arrows_to_store ) ) {
			}

			/// \brief TODOCUMENT
			inline void masked_bests_cacher::advance_to_pos_with_best_so_far(const res_arrow         &arg_new_position, ///< TODOCUMENT
			                                                                 const scored_arch_proxy &arg_best_so_far   ///< TODOCUMENT
			                                                                 ) {
				// std::cerr << "Number of arrows : " << arrows_to_store.size() << "\n";
				// std::cerr << "New position     : " << arg_new_position << "\n";
				// std::cerr << "current_itr      : " << current_itr << "\n";
				// std::cerr << "current_itr      : " << current_itr << "\n";
				// const auto var = std::find_if(
				// 	current_itr,
				// 	common::cend( arrows_to_store ),
				// 	[&] (const res_arrow &x) { return x >= arg_new_position; }
				// );
				advance_to_itr_with_best_so_far(
					std::find_if(
						current_itr,
						common::cend( arrows_to_store ),
						[&] (const res_arrow &x) { return x >= arg_new_position; }
					),
					arg_best_so_far
				);
			}

			/// \brief TODOCUMENT
			inline void masked_bests_cacher::advance_to_end_with_best_so_far(const scored_arch_proxy &arg_best_so_far ///< TODOCUMENT
			                                                                 ) {
				advance_to_itr_with_best_so_far(
					common::cend( arrows_to_store ),
					arg_best_so_far
				);
			}

			masked_bests_cacher make_masked_bests_cacher(masked_bests_cache &,
			                                             const hit_vec &,
			                                             const hit_list &,
			                                             const res_arrow &);

			res_arrow_vec get_arrows_before_starts_of_doms_right_interspersed_with_all_of(const hit_vec &,
			                                                                              const hit_list &,
			                                                                              const res_arrow &);

		}
	}
}

#endif