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

#ifndef POPULATE_MATRIX_SCAN_ACTION_H_INCLUDED
#define POPULATE_MATRIX_SCAN_ACTION_H_INCLUDED

#include <boost/filesystem/path.hpp>

#include "scan/detail/res_pair/single_struc_res_pair.h"
#include "scan/scan_index.h"
#include "scan/scan_query_set.h"

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

		public: // ***** TEMPORARILY PUBLIC *****
			double & get_entry(const index_type &, const index_type &);
			const double & get_entry(const index_type &, const index_type &) const;

		public:
			populate_matrix_scan_action(const index_type &arg_num_residues_a,
			                            const index_type &arg_num_residues_b,
			                            const index_type &arg_structure_a = 0,
			                            const index_type &arg_structure_b = 0);
			void operator()(const detail::single_struc_res_pair &,
			                const detail::single_struc_res_pair &,
			                const index_type &,
			                const index_type &);


			void plot_to_file(const boost::filesystem::path &,
			                  align::detail::matrix_plotter &) const;

			index_type get_length_a() const;
			index_type get_length_b() const;
//			new_matrix_dyn_prog_score_source
		};

		/// \brief TODOCUMENT
		inline double & populate_matrix_scan_action::get_entry(const index_type &arg_structure_a, ///< TODOCUMENT
		                                                       const index_type &arg_structure_b  ///< TODOCUMENT
		                                                       ) {
			return const_cast<double &>( static_cast<const populate_matrix_scan_action &>( *this ).get_entry( arg_structure_a, arg_structure_b ) );
		}

		/// \brief TODOCUMENT
		inline const double & populate_matrix_scan_action::get_entry(const index_type &arg_structure_a, ///< TODOCUMENT
		                                                             const index_type &arg_structure_b  ///< TODOCUMENT
		                                                             ) const {
			return the_matrix[ arg_structure_a ][ arg_structure_b ];
		}

		/// \brief TODOCUMENT
		inline void populate_matrix_scan_action::operator()(const detail::single_struc_res_pair &arg_res_pair_a,  ///< TODOCUMENT
		                                                    const detail::single_struc_res_pair &arg_res_pair_b,  ///< TODOCUMENT
		                                                    const index_type                    &arg_structure_a, ///< TODOCUMENT
		                                                    const index_type                    &arg_structure_b  ///< TODOCUMENT
		                                                    ) {
			if ( arg_structure_a == structure_a && arg_structure_b == structure_b ) {
				const auto distance   = sqrt( detail::squared_distance( arg_res_pair_a, arg_res_pair_b ) );
				// Add half the score to the "from" cell and half to the "to" cell
				const auto half_score = 0.5 - ( distance / 7.0 );
				get_entry( arg_res_pair_a.get_from_res_idx(), arg_res_pair_b.get_from_res_idx() ) += half_score;
				get_entry( arg_res_pair_a.get_to_res_idx(),   arg_res_pair_b.get_to_res_idx()   ) += half_score;
//				std::cerr << "In populate_matrix_scan_action[" << this << "], adding score " << score
//				          << " to "     << arg_res_pair_a.get_from_res_idx()
//						  << ","        << arg_res_pair_b.get_from_res_idx()
//						  << " and to " << arg_res_pair_a.get_to_res_idx()
//						  << ","        << arg_res_pair_b.get_to_res_idx()
//						  << "\n";
			}
//			std::cerr << get_entry( arg_res_pair_a.get_from_res_idx(), arg_res_pair_b.get_from_res_idx() ) << "\n";
//			std::cerr << get_entry( arg_res_pair_a.get_to_res_idx(),   arg_res_pair_b.get_to_res_idx() ) << "\n";
		}

		void gnuplot_to_file(const populate_matrix_scan_action &,
		                     const boost::filesystem::path &);

		/// \brief TODOCUMENT
		template <typename... KPs>
		populate_matrix_scan_action make_populate_matrix_scan_action(const scan_query_set<KPs...> &arg_query_set,             ///< TODOCUMENT
		                                                             const scan_index<KPs...>     &arg_index,                 ///< TODOCUMENT
		                                                             const index_type             &arg_query_structure_index, ///< TODOCUMENT
		                                                             const index_type             &arg_index_structure_index  ///< TODOCUMENT
		                                                             ) {
			return {
				arg_query_set.get_num_residues_of_structure_of_index( arg_query_structure_index ),
				arg_index.get_num_residues_of_structure_of_index    ( arg_index_structure_index ),
				arg_query_structure_index,
				arg_index_structure_index
			};
		}


	}
}

#endif
