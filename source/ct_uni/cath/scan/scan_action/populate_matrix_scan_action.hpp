/// \file
/// \brief The populate_matrix_scan_action class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_SCAN_ACTION_POPULATE_MATRIX_SCAN_ACTION_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_SCAN_ACTION_POPULATE_MATRIX_SCAN_ACTION_HPP

#include <filesystem>

#include "cath/scan/detail/res_pair/single_struc_res_pair.hpp"
#include "cath/scan/scan_index.hpp"
#include "cath/scan/scan_query_set.hpp"

namespace cath { namespace align { namespace detail { class matrix_plotter; } } }

namespace cath {
	namespace scan {

		/// \brief TODOCUMENT
		class populate_matrix_scan_action final {
		private:
			/// \brief TODOCUMENT
			doub_vec_vec the_matrix;

			/// \brief TODOCUMENT
			index_type structure_a;

			/// \brief TODOCUMENT
			index_type structure_b;

			/// \brief const-agnostic implementation of get_entry()
			///
			/// See GSL rule: Pro.Type.3: Don't use const_cast to cast away const (i.e., at all)
			/// (https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md#Pro-type-constcast)
			template <typename Action>
			static auto get_entry_impl(Action           &prm_action,      ///< TODOCUMENT
			                           const index_type &prm_structure_a, ///< TODOCUMENT
			                           const index_type &prm_structure_b  ///< TODOCUMENT
			                           ) -> decltype( prm_action.get_entry( prm_structure_a, prm_structure_b ) ) {
				return prm_action.the_matrix[ prm_structure_a ][ prm_structure_b ];
			}

		public: // ***** TEMPORARILY PUBLIC *****
			double & get_entry(const index_type &, const index_type &);
			const double & get_entry(const index_type &, const index_type &) const;

		public:
			populate_matrix_scan_action(const index_type &prm_num_residues_a,
			                            const index_type &prm_num_residues_b,
			                            const index_type & = 0,
			                            const index_type & = 0);
			void operator()(const detail::single_struc_res_pair &,
			                const detail::single_struc_res_pair &,
			                const index_type &,
			                const index_type &);


			void plot_to_file(const ::std::filesystem::path &,
			                  align::detail::matrix_plotter &) const;

			index_type get_length_a() const;
			index_type get_length_b() const;
//			new_matrix_dyn_prog_score_source
		};

		/// \brief TODOCUMENT
		inline double & populate_matrix_scan_action::get_entry(const index_type &prm_structure_a, ///< TODOCUMENT
		                                                       const index_type &prm_structure_b  ///< TODOCUMENT
		                                                       ) {
			return get_entry_impl( *this, prm_structure_a, prm_structure_b );
		}

		/// \brief TODOCUMENT
		inline const double & populate_matrix_scan_action::get_entry(const index_type &prm_structure_a, ///< TODOCUMENT
		                                                             const index_type &prm_structure_b  ///< TODOCUMENT
		                                                             ) const {
			return get_entry_impl( *this, prm_structure_a, prm_structure_b );
		}

		/// \brief TODOCUMENT
		inline void populate_matrix_scan_action::operator()(const detail::single_struc_res_pair &prm_res_pair_a,  ///< TODOCUMENT
		                                                    const detail::single_struc_res_pair &prm_res_pair_b,  ///< TODOCUMENT
		                                                    const index_type                    &prm_structure_a, ///< TODOCUMENT
		                                                    const index_type                    &prm_structure_b  ///< TODOCUMENT
		                                                    ) {
			if ( prm_structure_a == structure_a && prm_structure_b == structure_b ) {
				const auto distance   = sqrt( detail::squared_distance( prm_res_pair_a, prm_res_pair_b ) );
				// Add half the score to the "from" cell and half to the "to" cell
				const auto half_score = 0.5 - ( distance / 7.0 );
				get_entry( prm_res_pair_a.get_from_res_idx(), prm_res_pair_b.get_from_res_idx() ) += half_score;
				get_entry( prm_res_pair_a.get_to_res_idx(),   prm_res_pair_b.get_to_res_idx()   ) += half_score;
//				std::cerr << "In populate_matrix_scan_action[" << this << "], adding score " << score
//				          << " to "     << prm_res_pair_a.get_from_res_idx()
//						  << ","        << prm_res_pair_b.get_from_res_idx()
//						  << " and to " << prm_res_pair_a.get_to_res_idx()
//						  << ","        << prm_res_pair_b.get_to_res_idx()
//						  << "\n";
			}
//			std::cerr << get_entry( prm_res_pair_a.get_from_res_idx(), prm_res_pair_b.get_from_res_idx() ) << "\n";
//			std::cerr << get_entry( prm_res_pair_a.get_to_res_idx(),   prm_res_pair_b.get_to_res_idx() ) << "\n";
		}

		void gnuplot_to_file(const populate_matrix_scan_action &,
		                     const ::std::filesystem::path &);

		/// \brief TODOCUMENT
		template <typename... KPs>
		populate_matrix_scan_action make_populate_matrix_scan_action(const scan_query_set<KPs...> &prm_query_set,             ///< TODOCUMENT
		                                                             const scan_index<KPs...>     &prm_index,                 ///< TODOCUMENT
		                                                             const index_type             &prm_query_structure_index, ///< TODOCUMENT
		                                                             const index_type             &prm_index_structure_index  ///< TODOCUMENT
		                                                             ) {
			return {
				prm_query_set.get_num_residues_of_structure_of_index( prm_query_structure_index ),
				prm_index.get_num_residues_of_structure_of_index    ( prm_index_structure_index ),
				prm_query_structure_index,
				prm_index_structure_index
			};
		}


	} // namespace scan
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCAN_SCAN_ACTION_POPULATE_MATRIX_SCAN_ACTION_HPP
