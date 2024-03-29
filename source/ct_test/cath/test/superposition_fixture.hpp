/// \file
/// \brief Test fixture for test suites involving superpositions

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

#ifndef CATH_TOOLS_SOURCE_CT_TEST_CATH_TEST_SUPERPOSITION_FIXTURE_HPP
#define CATH_TOOLS_SOURCE_CT_TEST_CATH_TEST_SUPERPOSITION_FIXTURE_HPP

#include <string>

#include "cath/chopping/chopping_type_aliases.hpp"
#include "cath/file/pdb/pdb.hpp"
#include "cath/file/pdb/pdb_atom.hpp"
#include "cath/file/pdb/pdb_list.hpp"
#include "cath/file/pdb/pdb_residue.hpp"
#include "cath/structure/geometry/coord_list.hpp"
#include "cath/superposition/superposition_context.hpp"
#include "cath/test/global_test_constants.hpp"

namespace cath::sup {

	using ::std::literals::string_literals::operator""s;

	/// \brief A test fixture for superposition tests that adds some extras to global_test_constants
	class superposition_fixture : protected global_test_constants {
	protected:
		static constexpr size_t NUM_ENTRIES = 2;

		const geom::coord_list      coord_list_1{ { geom::coord{  0.0,  2.0, -1.0 }, geom::coord{  0.0,  2.0,  1.0 }, geom::coord{  0.0, -4.0,  0.0 } } };
		const geom::coord_list      coord_list_2{ { geom::coord{ -1.0,  0.0, -2.0 }, geom::coord{  1.0,  0.0, -2.0 }, geom::coord{  0.0,  0.0,  4.0 } } };

		const file::pdb_list        pdbs{ file::pdb_vec{ NUM_ENTRIES, file::pdb{} } };
		const file::name_set_list   names{ file::name_set_vec{ { file::name_set{ "1c0pA01"s }, file::name_set{ "1hdoA00"s } } } };
		const superposition         the_sup{ create_pairwise_superposition( coord_list_1, coord_list_2 ) };
		const superposition_context the_sup_con{ the_sup, pdbs, names, cath::chop::region_vec_opt_vec( NUM_ENTRIES ) };

		const std::string sup_json_str = R"({"transformations":[{"translation":)"
			R"({"x":"0","y":"0","z":"0"},)"
			R"("rotation":)"
			R"([["1","0","0"],)"
			R"(["0","1","0"],)"
			R"(["0","0","1"])"
			R"(]},{"translation":)"
			R"({"x":"0","y":"0","z":"0"},)"
			R"("rotation":)"
			R"([["0","-1","0"],)"
			R"(["0","0","-1"],)"
			R"(["1","0","0"])"
			R"(]}]})"
			"\n";
		const std::string sup_context_json_str = R"({"entries":[{"name":"1c0pA01","transformation":{"translation":)"
			R"({"x":"0","y":"0","z":"0"},)"
			R"("rotation":)"
			R"([["1","0","0"],)"
			R"(["0","1","0"],)"
			R"(["0","0","1"])"
			R"(]}},{"name":"1hdoA00","transformation":{"translation":)"
			R"({"x":"0","y":"0","z":"0"},)"
			R"("rotation":)"
			R"([["0","-1","0"],)"
			R"(["0","0","-1"],)"
			R"(["1","0","0"])"
			R"(]}}]})" "\n";
	};

} // namespace cath::sup

#endif // CATH_TOOLS_SOURCE_CT_TEST_CATH_TEST_SUPERPOSITION_FIXTURE_HPP
