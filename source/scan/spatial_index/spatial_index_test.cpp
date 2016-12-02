/// \file
/// \brief The spatial_index test suite

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

#include <boost/test/auto_unit_test.hpp>

#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/irange.hpp>
// #include <boost/log/trivial.hpp> /// ***** TEMPORARY *****

#include "common/boost_addenda/range/front.hpp"
#include "common/boost_addenda/range/utility/iterator/cross_itr.hpp"
#include "common/size_t_literal.hpp"
#include "common/tuple/tuple_mins_maxs_element.hpp"
#include "scan/detail/scan_index_store/scan_index_lattice_store.hpp"
#include "scan/detail/scan_type_aliases.hpp"
#include "scan/res_pair_keyer/res_pair_keyer.hpp"
#include "scan/res_pair_keyer/res_pair_keyer_part/detail/axis_keyer_part.hpp"
#include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_view_x_keyer_part.hpp"
#include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_view_y_keyer_part.hpp"
#include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_view_z_keyer_part.hpp"
#include "ssap/context_res.hpp"
#include "structure/geometry/coord.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/protein_io.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"
#include "test/global_test_constants.hpp"
// #include "common/chrono/duration_to_seconds_string.hpp"
// #include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_from_phi_keyer_part.hpp"
// #include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_from_psi_keyer_part.hpp"
// #include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_index_dirn_keyer_part.hpp"
// #include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_to_phi_keyer_part.hpp"
// #include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_to_psi_keyer_part.hpp"
// #include "scan/scan_action/populate_matrix_scan_action.hpp" // ***** TEMPORARY? *****
// #include "scan/scan_action/record_scores_scan_action.hpp"   // ***** TEMPORARY? *****
// #include "scan/scan_query_set.hpp"
// #include "scan/spatial_index.hpp"
// #include "structure/geometry/angle.hpp"
// #include "structure/protein/protein_source_file_set/protein_source_from_pdb.hpp"
// //#include "scan/detail/res_pair_dirn/res_pair_dirn.hpp"
// //#include "scan/res_pair_keyer/res_pair_keyer_part/res_pair_orient_keyer_part.hpp"

// #include <chrono> // ***** TEMPORARY? *****

namespace cath { namespace test { } }

using namespace cath;
using namespace cath::common;
using namespace cath::geom;
using namespace cath::scan;
using namespace cath::scan::detail;
using namespace cath::test;
// using namespace std::chrono; // ***** TEMPORARY? *****

using boost::adaptors::transformed;
using boost::algorithm::join;
using boost::irange;
using std::make_tuple;

// template <typename T> class TD;
// template <size_t T> class SD;
// template <int    T> class ID;

namespace cath {
	namespace scan {

		/// \brief TODOCUMENT
		struct simple_locn_index final {
			/// \brief TODOCUMENT
			view_base_type view_x;

			/// \brief TODOCUMENT
			view_base_type view_y;

			/// \brief TODOCUMENT
			view_base_type view_z;

			/// \brief TODOCUMENT
			unsigned int index;
		};

		/// \brief Convenience function to get the x component of the view in the specified simple_locn_index
		///
		/// \relates simple_locn_index
		inline constexpr const view_base_type & get_view_x(const simple_locn_index &arg_res_pair ///< The simple_locn_index to query
		                                                   ) {
		        return arg_res_pair.view_x;
		}

		/// \brief Convenience function to get the y component of the view in the specified simple_locn_index
		///
		/// \relates simple_locn_index
		inline constexpr const view_base_type & get_view_y(const simple_locn_index &arg_res_pair ///< The simple_locn_index to query
		                                                   ) {
		        return arg_res_pair.view_y;
		}

		/// \brief Convenience function to get the z component of the view in the specified simple_locn_index
		///
		/// \relates simple_locn_index
		inline constexpr const view_base_type & get_view_z(const simple_locn_index &arg_res_pair ///< The simple_locn_index to query
		                                                   ) {
		        return arg_res_pair.view_z;
		}

		/// \brief TODOCUMENT
		simple_locn_index make_simple_locn_index(const coord        &arg_coord, ///< TODOCUMENT
		                                         const unsigned int &arg_index  ///< TODOCUMENT
		                                         ) {
			return {
				debug_numeric_cast< view_base_type>( arg_coord.get_x() ),
				debug_numeric_cast< view_base_type>( arg_coord.get_y() ),
				debug_numeric_cast< view_base_type>( arg_coord.get_z() ),
				arg_index
			};
		}

		/// \brief TODOCUMENT
		simple_locn_index make_simple_locn_index_of_ca(const residue      &arg_res,   ///< TODOCUMENT
		                                               const unsigned int &arg_index  ///< TODOCUMENT
		                                               ) {
			return make_simple_locn_index(
				arg_res.get_carbon_alpha_coord(),
				arg_index
			);
		}

		/// \brief TODOCUMENT
		struct simple_locn_crit final {

			/// \brief TODOCUMENT
			detail::view_base_type maximum_squared_distance;
		};

		std::string to_string(const simple_locn_index &arg_simple_locn_index
		                      ) {
			return "simple_locn_index[ ("
				+ ::std::to_string( arg_simple_locn_index.view_x )
				+ ", "
				+ ::std::to_string( arg_simple_locn_index.view_y )
				+ ", "
				+ ::std::to_string( arg_simple_locn_index.view_z )
				+ "), "
				+ ::std::to_string( arg_simple_locn_index.index )
				+ "]";
		}

			/// \brief TODOCUMENT
		detail::view_base_type get_maximum_distance(const simple_locn_crit &arg_crit ///< TODOCUMENT
		                                            ) {
			return std::sqrt( arg_crit.maximum_squared_distance );
		}

		using simple_locn_x_keyer_part = detail::axis_keyer_part<detail::res_pair_view_x_keyer_part_spec< simple_locn_index, simple_locn_crit > >;
		using simple_locn_y_keyer_part = detail::axis_keyer_part<detail::res_pair_view_y_keyer_part_spec< simple_locn_index, simple_locn_crit > >;
		using simple_locn_z_keyer_part = detail::axis_keyer_part<detail::res_pair_view_z_keyer_part_spec< simple_locn_index, simple_locn_crit > >;
		// SD< sizeof( simple_locn_index ) > fred;

	} // namespace scan
} // namespace cath

// using doub_doub_doub_tpl = std::tuple<double, double, double>;
// BOOST_TEST_DONT_PRINT_LOG_VALUE( doub_doub_doub_tpl )

// using float_float_float_tpl = std::tuple<float, float, float>;
// BOOST_TEST_DONT_PRINT_LOG_VALUE( float_float_float_tpl )

namespace cath {
	namespace test {

		using keyer_type = res_pair_keyer<simple_locn_x_keyer_part, simple_locn_y_keyer_part, simple_locn_z_keyer_part>;

		/// \brief The spatial_index_test_suite_fixture to assist in testing spatial_index
		struct spatial_index_test_suite_fixture : protected global_test_constants {
		protected:
			~spatial_index_test_suite_fixture() noexcept = default;

			const protein protein_a = read_protein_from_dssp_and_pdb( EXAMPLE_A_DSSP_FILENAME(), EXAMPLE_A_PDB_FILENAME(), true, EXAMPLE_A_PDB_STEMNAME() );
			const protein protein_b = read_protein_from_dssp_and_pdb( EXAMPLE_B_DSSP_FILENAME(), EXAMPLE_B_PDB_FILENAME(), true, EXAMPLE_B_PDB_STEMNAME() );

			static constexpr float CELL_SIZE = 10.0;
			static constexpr float RADIUS    =  7.0;

			// static constexpr simple_locn_x_keyer_part x_keyer{ CELL_SIZE };
			// static constexpr simple_locn_y_keyer_part y_keyer{ CELL_SIZE };
			// static constexpr simple_locn_z_keyer_part z_keyer{ CELL_SIZE };
			static constexpr keyer_type the_keyer = make_res_pair_keyer(
				simple_locn_x_keyer_part{ CELL_SIZE },
				simple_locn_y_keyer_part{ CELL_SIZE },
				simple_locn_z_keyer_part{ CELL_SIZE }
			);
		};

		constexpr float                    spatial_index_test_suite_fixture::CELL_SIZE;
		constexpr float                    spatial_index_test_suite_fixture::RADIUS;
		// constexpr simple_locn_x_keyer_part spatial_index_test_suite_fixture::x_keyer;
		// constexpr simple_locn_y_keyer_part spatial_index_test_suite_fixture::y_keyer;
		// constexpr simple_locn_z_keyer_part spatial_index_test_suite_fixture::z_keyer;
		constexpr keyer_type               spatial_index_test_suite_fixture::the_keyer;

	} // namespace test
} // namespace cath

namespace cath {
	namespace scan {

		/// \brief TODOCUMENT
		template <typename Rng, typename Keyer>
		auto make_empty_lattice_store(const Rng   &arg_rng,  ///< TODOCUMENT
		                              const Keyer &arg_keyer ///< TODOCUMENT
		                              ) {
			using value_t = range_value_t< Rng >;

			const auto mins_maxs = tuple_mins_maxs_element(
				arg_rng
					| transformed( [&] (const value_t &value) { return arg_keyer.make_key( value ); } )
			);

			return make_scan_index_lattice_store< std::vector<value_t> >( mins_maxs.first, mins_maxs.second );
		}

		/// \brief TODOCUMENT
		template <typename Rng, typename Keyer>
		auto make_sparse_lattice_store(const Rng   &arg_rng,  ///< TODOCUMENT
		                               const Keyer &arg_keyer ///< TODOCUMENT
		                               ) {
			auto the_store = make_empty_lattice_store( arg_rng, arg_keyer );
			for (const auto &value : arg_rng) {
				the_store.push_back_entry_to_cell( arg_keyer.make_key( value ), value );
			}
			return the_store;
		}

		/// \brief TODOCUMENT
		template <typename Rng, typename Keyer, typename Crit>
		auto make_dense_lattice_store(const Rng   &arg_rng,   ///< TODOCUMENT
		                              const Keyer &arg_keyer, ///< TODOCUMENT
		                              const Crit  &arg_crit   ///< TODOCUMENT
		                              ) {
			auto the_store = make_empty_lattice_store( arg_rng, arg_keyer );
			for (const auto &value : arg_rng) {
				const auto close_keys = arg_keyer.make_close_keys( value, arg_crit );
				for (const auto &the_key : common::cross( close_keys ) ) {
					the_store.push_back_entry_to_cell( the_key, value );
				}
			}
			return the_store;
		}

	} // namespace scan
} // namespace cath


BOOST_FIXTURE_TEST_SUITE(spatial_index_test_suite, spatial_index_test_suite_fixture)


BOOST_AUTO_TEST_CASE(basic) {
	// // ID< std::get<2>( the_keyer.make_key  ( simple_locn_index{  1.0, -9.0, 31.0, 49 } ) ) > bob;

	// constexpr simple_locn_crit the_crit{ RADIUS * RADIUS };

	// static_assert( the_keyer.make_value( simple_locn_index{   1.0,  -9.0,  31.0, 49 } ) == make_tuple(  1.0, -9.0, 31.0 ), "" );
	// static_assert( the_keyer.make_key  ( simple_locn_index{   1.0,  -9.0,  31.0, 49 } ) == make_tuple(    0,   -1,    3 ), "" );
	// static_assert( the_keyer.make_key  ( simple_locn_index{   0.0, -10.0,  10.0, 49 } ) == make_tuple(    0,   -1,    1 ), "" );

	// const auto the_close_keys = the_keyer.make_close_keys( simple_locn_index{   1.0,  -9.0,  31.0, 49 }, the_crit );
	// for (const auto &crossed_val : cross( the_close_keys ) ) {
	// 	std::cerr
	// 		<< std::get<0>( crossed_val )
	// 		<< " "
	// 		<< std::get<1>( crossed_val )
	// 		<< " "
	// 		<< std::get<2>( crossed_val )
	// 		<< "\n";
	// }

	// constexpr simple_locn_x_keyer_part x_keyer{ CELL_SIZE };
	// constexpr simple_locn_y_keyer_part y_keyer{ CELL_SIZE };
	// constexpr simple_locn_z_keyer_part z_keyer{ CELL_SIZE };

	// // constexpr auto the_keyer = make_res_pair_keyer( x_keyer, y_keyer, z_keyer );

	// // key_tuple_type make_key(const detail::multi_struc_res_rep_pair &) const;
	// // key_ranges_tuple_type make_close_keys(const detail::multi_struc_res_rep_pair &, const quad_criteria &) const;

	// const auto sls_a = make_sparse_lattice_store(
	// 	irange( 0_z, protein_a.get_length() )
	// 		| transformed( [&] (const size_t &x) {
	// 			return make_simple_locn_index_of_ca( protein_a.get_residue_ref_of_index( x ), debug_numeric_cast<unsigned int>( x ) );
	// 		} ),
	// 	the_keyer
	// );
	// std::cerr << "protein_a.get_length() : " << protein_a.get_length() << "\n";
	// std::cerr << "sls_a.get_info_size()  : " << ::std::to_string( sls_a.get_info_size().value() ) << "b" << "\n";

	// size_t full_count = 0; 
	// for (const size_t &the_res_idx : irange( 0_z, protein_a.get_length() ) ) {
		
	// 	const auto &the_res = protein_a.get_residue_ref_of_index( the_res_idx );

	// 	const auto  data    = make_simple_locn_index_of_ca( the_res, debug_numeric_cast<unsigned int>( the_res_idx ) );

	// 	std::cerr << "Searching for residues close to " << to_string( data ) << "\n";
	// 	// const auto  key     = the_keyer.make_key( data );
	// 	size_t count = 0;
	// 	for (const auto &key : common::cross( the_keyer.make_close_keys( data, the_crit ) ) ) {
	// 		// std::cerr << "key is " << std::get<0>( key ) << ", " << std::get<1>( key ) << ", " << std::get<2>( key ) << "\n";
	// 		if ( sls_a.has_matches( key ) ) {
	// 			for (const simple_locn_index &eg : sls_a.find_matches( key ) ) {
	// 				boost::ignore_unused( eg );
	// 				// std::cerr << to_string( eg ) << "\n";
	// 				++count;
	// 				++full_count;
	// 			}
	// 		}
	// 	}
	// 	std::cerr << "count is " << count << "\n";
	// }
	// std::cerr << "full_count is " << full_count << "\n";


	// std::cerr << the_keyer.parts_names() << "\n";;

	// const simple_locn_y_keyer_part fred{ CELL_SIZE };

	// std::cerr << "Name : " << fred.get_name() << "\n";

	// const auto search_radius = fred.get_search_radius( the_crit );
	// std::cerr << "Crit : " << search_radius << "\n";

	// const auto &first_res_a = front( protein_a );
	// for (const size_t &other_res_a_idx : irange( 0_z, protein_a.get_length() ) ) {
	// 	const auto &other_res_a = protein_a.get_residue_ref_of_index( other_res_a_idx );
	// 	const simple_locn_index bob = make_simple_locn_index(
	// 		view_vector_of_residue_pair( first_res_a, other_res_a ),
	// 		static_cast<unsigned int>( other_res_a_idx )
	// 	);
	// 	const auto value = fred.get_value( bob );

	// 	std::cerr << bob.index << "\t" << coord{ bob.view_x, bob.view_y, bob.view_z };

	// 	std::cerr << ", value is : " << value;
	// 	std::cerr << ", key_part is : " << fred.key_part( value );
	// 	std::cerr << ", close_key_parts is :";
	// 	for (const auto &x : fred.close_key_parts( value, search_radius ) ) {
	// 		std::cerr << " " << x;
	// 	}
	// 	std::cerr << "\n";
	// 	std::cerr << " " << std::get<0>( the_keyer.make_value( bob ) );
	// 	std::cerr << " " << std::get<1>( the_keyer.make_value( bob ) );
	// 	std::cerr << " " << std::get<2>( the_keyer.make_value( bob ) );
	// 	std::cerr << "\n";
	// 	std::cerr << " " << std::get<0>( the_keyer.make_key( bob ) );
	// 	std::cerr << " " << std::get<1>( the_keyer.make_key( bob ) );
	// 	std::cerr << " " << std::get<2>( the_keyer.make_key( bob ) );
	// 	std::cerr << "\n\n";
	// }

	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
