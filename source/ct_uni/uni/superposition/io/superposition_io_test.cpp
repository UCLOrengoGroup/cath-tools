///// \file
///// \brief The superposition I/O test suite

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

#include <boost/filesystem/path.hpp>
#include <boost/test/unit_test.hpp>

#include "chopping/region/region.hpp"
#include "common/boost_addenda/log/stringstream_log_sink.hpp"
#include "common/file/temp_file.hpp"
#include "common/property_tree/from_json_string.hpp"
#include "common/property_tree/to_json_string.hpp"
#include "superposition/io/superposition_io.hpp"
#include "test/superposition_fixture.hpp"

#include <regex>
#include <sstream>

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::geom;
using namespace cath::sup;

using ::boost::filesystem::path;
using ::std::ostringstream;
using ::std::regex;

/// \brief The superposition_io_test_suite_fixture to assist in testing superposition I/O
struct superposition_io_test_suite_fixture : protected superposition_fixture{
protected:
	~superposition_io_test_suite_fixture() noexcept = default;
};

BOOST_FIXTURE_TEST_SUITE(superposition_io_test_suite, superposition_io_test_suite_fixture)

BOOST_AUTO_TEST_SUITE(json)

BOOST_AUTO_TEST_SUITE(write)

BOOST_AUTO_TEST_CASE(to_json_string_works_for_example_sup) {
	BOOST_CHECK_EQUAL( to_json_string( the_sup, json_style::COMPACT ), sup_json_str );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(write)

BOOST_AUTO_TEST_CASE(from_json_string_works_for_identity) {
	BOOST_CHECK_EQUAL( from_json_string<superposition>( sup_json_str ), the_sup );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE(warns_if_overwriting_chain_label_in_pdb_with_multiple_chains) {
	const stringstream_log_sink log_sink;
	const temp_file temp_file("cath_tools_test_temp_file.superposition_chain_munge_warning.%%%%");
	const auto the_pdbs = make_pdb_list( { read_pdb_file( TEST_SOURCE_DATA_DIR() / "supn_content" / "1bdh" ) } );

	write_superposed_pdb_to_file(
		make_identity_superposition_of( the_pdbs ),
		get_filename( temp_file ),
		the_pdbs,
		sup_pdbs_script_policy::LEAVE_RAW_PDBS,
		chain_relabel_policy::RELABEL
	);

	BOOST_CHECK( regex_search(
		log_sink.str(),
		regex{ R"(Overwriting chain labels.*contains multiple chain labels)" }
	) );
}

BOOST_AUTO_TEST_CASE(does_not_warn_if_overwriting_chain_label_in_pdb_with_single_chain) {
	const stringstream_log_sink log_sink;
	const temp_file temp_file("cath_tools_test_temp_file.superposition_chain_munge_non_warning.%%%%");
	const auto the_pdbs = make_pdb_list( { read_pdb_file( TEST_SOURCE_DATA_DIR() / "1c0pA01" ) } );

	write_superposed_pdb_to_file(
		make_identity_superposition_of( the_pdbs ),
		get_filename( temp_file ),
		the_pdbs,
		sup_pdbs_script_policy::LEAVE_RAW_PDBS,
		chain_relabel_policy::RELABEL
	);

	BOOST_CHECK( log_sink.str_is_empty() );
}


BOOST_AUTO_TEST_SUITE_END()
