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

#include "spatial_index.hpp"

#include <boost/test/unit_test.hpp>

#include "test/global_test_constants.hpp"

namespace cath { namespace test { } }

using namespace cath::scan;
using namespace cath::test;

using std::make_tuple;
using std::tuple;

using short_short_short_tuple = tuple<short int, short int, short int>;

BOOST_TEST_DONT_PRINT_LOG_VALUE( short_short_short_tuple )

namespace cath {
	namespace test {

		/// \brief The spatial_index_test_suite_fixture to assist in testing spatial_index
		struct spatial_index_test_suite_fixture : protected global_test_constants {
		protected:
			~spatial_index_test_suite_fixture() noexcept = default;
		};

	} // namespace test
} // namespace cath


BOOST_FIXTURE_TEST_SUITE(spatial_index_test_suite, spatial_index_test_suite_fixture)

BOOST_AUTO_TEST_CASE(basic) {
	constexpr auto index_x_y_z_keyer = make_res_pair_keyer(
		simple_locn_x_keyer_part{ 10.0 },
		simple_locn_y_keyer_part{ 10.0 },
		simple_locn_z_keyer_part{ 10.0 },
		res_pair_from_to_index_keyer_part{}
	);

	static_assert(
		index_x_y_z_keyer.make_value( simple_locn_index{ 12.0, 34.0, -22.0, 4u } ) == make_tuple( 12.0, 34.0, -22.0, 4u ),
		""
	);

	// Making this test a static_assert means that a floor() has to be switched for a constexpr_floor() which
	// appears to be quite a bit slower
	// static_assert(
	// 	index_x_y_z_keyer.make_key  ( simple_locn_index{ 12.0, 34.0, -22.0, 4u } ) == make_tuple(  1,    3,    -3       ),
	// 	""
	// );

	BOOST_CHECK_EQUAL( index_x_y_z_keyer.make_key  ( simple_locn_index{ 12.0, 34.0, -22.0, 4u } ), short_short_short_tuple(  1,    3,    -3       ) );

	// for (const auto &spatial_try_val : { spatial_try::SPARSE_FROM_PROT,
	//                                      // spatial_try::SPARSE_FROM_SCAN,
	//                                      // spatial_try::DENSE_FROM_PROT,
	//                                      // spatial_try::DENSE_FROM_SCAN
	//                                  } ) {
	// 	// 18 81023 308584 0.000922 seconds 0.006032 seconds

	// 	// for (const float cell_size : { 1.25f, 2.5f, 5.0f, 10.0f, 15.0f, 20.0f } ) {
	// 	for (const float cell_size : irange( 5.0f, 30.0f ) ) {
	// 	// for (const float cell_size : irange( 18.0f, 19.0f ) ) {

	// 		const auto make_start_time = high_resolution_clock::now();
	// 		// const auto lattice         = make_sparse_lattice( protein_a, cell_size );
	// 		const auto lattice = [&] {
	// 			if ( spatial_try_val == spatial_try::SPARSE_FROM_PROT || spatial_try_val == spatial_try::SPARSE_FROM_SCAN ) {
	// 				return make_sparse_lattice ( protein_a, cell_size );
	// 			}
	// 			else {
	// 				return make_dense_lattice ( protein_a, cell_size, MAX_DIST );
	// 			}
	// 		} ();
	// 		const auto make_durn       = high_resolution_clock::now() - make_start_time;

	// 		size_t count = 0;

	// 		const auto scan_start_time = high_resolution_clock::now();

	// 		switch ( spatial_try_val ) {
	// 			case ( spatial_try::SPARSE_FROM_PROT ) : { scan_sparse_lattice        ( lattice, protein_a, cell_size, MAX_DIST, [&] (const simple_locn_index &x, const simple_locn_index &y) { if ( x.index < y.index ) { ++count; } } ) ; break; }
	// 			case ( spatial_try::SPARSE_FROM_SCAN ) : { scan_sparse_lattice_locally( lattice, protein_a, cell_size, MAX_DIST, [&] (const simple_locn_index &x, const simple_locn_index &y) { if ( x.index < y.index ) { ++count; } } ) ; break; }
	// 			case ( spatial_try::DENSE_FROM_PROT  ) : { scan_dense_lattice         ( lattice, protein_a, cell_size, MAX_DIST, [&] (const simple_locn_index &x, const simple_locn_index &y) { if ( x.index < y.index ) { ++count; } } ) ; break; }
	// 			case ( spatial_try::DENSE_FROM_SCAN  ) : { scan_dense_lattice_locally ( lattice, protein_a, cell_size, MAX_DIST, [&] (const simple_locn_index &x, const simple_locn_index &y) { if ( x.index < y.index ) { ++count; } } ) ; break; }
	// 		}

	// 		const auto scan_durn = high_resolution_clock::now() - scan_start_time;

	// 		// std::cerr << "cell_size    : " << cell_size << "\n";
	// 		// std::cerr << "count is     : " << count << "\n";
	// 		// std::cerr << "lattice size : " << lattice.get_info_size().value() << "b\n";
	// 		// std::cerr << "make_time    : " << durn_to_seconds_string( make_durn ) << "\n";
	// 		// std::cerr << "scan_durn    : " << durn_to_seconds_string( scan_durn ) << "\n";
	// 		// std::cerr << "\n";
	// 		switch ( spatial_try_val ) {
	// 			case ( spatial_try::SPARSE_FROM_PROT ) : { std::cerr << "SPARSE_FROM_PROT " ; break; }
	// 			case ( spatial_try::SPARSE_FROM_SCAN ) : { std::cerr << "SPARSE_FROM_SCAN " ; break; }
	// 			case ( spatial_try::DENSE_FROM_PROT  ) : { std::cerr << "DENSE_FROM_PROT  " ; break; }
	// 			case ( spatial_try::DENSE_FROM_SCAN  ) : { std::cerr << "DENSE_FROM_SCAN  " ; break; }
	// 		}

	// 		std::cerr        << cell_size;
	// 		std::cerr << " " << count;
	// 		std::cerr << " " << lattice.get_info_size().value();
	// 		std::cerr << " " << durn_to_seconds_string ( make_durn );
	// 		std::cerr << " " << durn_to_seconds_string ( scan_durn );
	// 		std::cerr << " " << durn_to_seconds_string ( make_durn + scan_durn );
	// 		std::cerr << " " << durn_to_rate_per_second( make_durn + scan_durn );
	// 		std::cerr << "\n";
	// 	}
	// }


	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
