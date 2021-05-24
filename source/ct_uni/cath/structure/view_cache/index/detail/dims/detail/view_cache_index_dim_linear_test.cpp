/// \file
/// \brief The view_cache_index_dim_linear test suite

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

#include <boost/test/unit_test.hpp>

// #include "cath/alignment/alignment.hpp"
// #include "cath/alignment/io/alignment_io.hpp"
// #include "cath/file/pdb/pdb.hpp"
// #include "cath/file/pdb/pdb_atom.hpp"
// #include "cath/file/pdb/pdb_list.hpp"
// #include "cath/file/pdb/pdb_residue.hpp"
// #include "cath/ssap/ssap.hpp"
// #include "cath/structure/geometry/angle.hpp"
// #include "cath/structure/protein/protein.hpp"
// #include "cath/structure/protein/protein_source_file_set/protein_from_pdb.hpp"
// #include "cath/structure/protein/sec_struc.hpp"
// #include "cath/structure/protein/sec_struc_planar_angles.hpp"
// #include "cath/structure/view_cache/index/detail/vcie_match_criteria.hpp"
// #include "cath/structure/view_cache/index/quad_find_action.hpp"
// #include "cath/structure/view_cache/index/view_cache_index_dim_linear.hpp"
// #include "cath/test/global_test_constants.hpp"

#include "cath/structure/view_cache/index/detail/dims/detail/view_cache_index_dim_linear.hpp"

// using namespace ::cath;
// using namespace ::cath::align;
// using namespace ::cath::index;
// using namespace ::cath::index::detail;
// using namespace ::cath::index::detail::detail;
using namespace ::cath::index::detail::detail::detail;
// using namespace ::cath::geom;
using namespace ::std;

// using ::boost::numeric_cast;

/// \todonow Tidy up this file

/// \brief TODOCUMENT
BOOST_AUTO_TEST_SUITE(view_cache_index_dim_linear_test_suite)

// template <size_t N> class TNS;
// template <int    N> class TNI;

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
	/// \todonow Put more sensible comments in the static asserts
	static_assert( clamped_cell_index_of_value_in_current( 20.0, 2, 10, 0.0 ) == 0, "This should work" );

	static_assert( search_cell_ranges( 20.0, 0, 10, 29.0,   31.0 ) == make_tuple( 1,  2, 0, 0 ), "This will not work yet" );
	static_assert( search_cell_ranges( 20.0, 0, 10, 30.0,   50.0 ) == make_tuple( 1,  3, 0, 0 ), "This will not work yet" );
	static_assert( search_cell_ranges( 20.0, 2, 10, 69.0,   71.0 ) == make_tuple( 1,  2, 0, 0 ), "This will not work yet" );

	static_assert( search_cell_ranges( 20.0, 2, 10,  0.0,   71.0 ) == make_tuple( 0,  2, 0, 0 ), "This will not work yet" );
	static_assert( search_cell_ranges( 20.0, 2, 10, 69.0, 1000.0 ) == make_tuple( 1, 10, 0, 0 ), "This will not work yet" );
	static_assert( search_cell_ranges( 20.0, 2, 10,  0.0, 1000.0 ) == make_tuple( 0, 10, 0, 0 ), "This will not work yet" );
	// static_assert( search_cell_ranges( 20.0, 0, 10, 78.0, 122.0 ) == make_tuple( 1, 2, 0, 0 ), "This will not work yet" );
	// static_assert( search_cell_ranges( 20.0, 0, 10, 78.0, 122.0 ) == make_tuple( 1, 2, 0, 0 ), "This will not work yet" );
	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()

