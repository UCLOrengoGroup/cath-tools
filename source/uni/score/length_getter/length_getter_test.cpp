/// \file
/// \brief The length_getter test suite

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

#include "score/length_getter/geometric_mean_length_getter.hpp"
#include "score/length_getter/length_getter_enum.hpp"
#include "score/length_getter/length_getter_types.hpp"
#include "score/length_getter/length_of_first_getter.hpp"
#include "score/length_getter/length_of_longer_getter.hpp"
#include "score/length_getter/length_of_second_getter.hpp"
#include "score/length_getter/length_of_shorter_getter.hpp"
#include "score/length_getter/mean_length_getter.hpp"
#include "score/length_getter/num_aligned_length_getter.hpp"
#include "structure/protein/protein.hpp"
#include "structure/protein/residue.hpp"
#include "structure/protein/sec_struc.hpp"
#include "structure/protein/sec_struc_planar_angles.hpp"
#include "test/boost_addenda/boost_check_no_throw_diag.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::score;
using namespace cath::score::detail;

using std::is_same;

namespace cath {
	namespace test {

		/// \brief The length_getter_test_suite_fixture to assist in testing length_getter
		struct length_getter_test_suite_fixture {
		protected:
			~length_getter_test_suite_fixture() noexcept = default;

			protein make_dummy_protein_of_n_residues(const size_t &);

			protein dummy_protein_10 = make_dummy_protein_of_n_residues( 10 );
			protein dummy_protein_90 = make_dummy_protein_of_n_residues( 90 );
		};

		/// \brief TODOCUMENT
		protein length_getter_test_suite_fixture::make_dummy_protein_of_n_residues(const size_t &prm_index ///< TODOCUMENT
		                                                                           ) {
			return build_protein( residue_vec{ prm_index, residue::NULL_RESIDUE } );
		}

//		/// \brief TODOCUMENT
//		template <size_t I>
//		class length_getter_enum_tester final {
//		public:
//			/// \brief TODOCUMENT
//			void operator()() {
//				constexpr length_getter_enum the_enum = get<I>( all_length_getter_enums );
//				using the_length_getter = typename length_getter_of_length_getter_enum<the_enum>::type;
//				static_assert( boost::mpl::contains< length_getter_types, the_length_getter>::value, "length_getter_types does not contain one of the length_getter types" );
//				constexpr length_getter_enum orig_enum = the_length_getter::enum_val;
//				static_assert( orig_enum == the_enum, "Length_getter enum round trip produced different enum value" );
//			}
//		};

	}  // namespace test
}  // namespace cath

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(length_getter_test_suite, cath::test::length_getter_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE_TEMPLATE(ctor_works, length_getter_type, length_getter_types) {
	BOOST_CHECK_NO_THROW_DIAG( length_getter_type() );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE_TEMPLATE(prot_only_get_length_does_not_throw, length_getter_type, protein_only_length_getter_types) {
	BOOST_CHECK_NO_THROW_DIAG( length_getter_type().get_prot_only_length( dummy_protein_10, dummy_protein_90 ) );
}

///// \brief TODOCUMENT
//BOOST_AUTO_TEST_CASE(static__length_getter_type_checks) {
//	BOOST_CHECK( true );
//	constexpr_for_n<cath::test::length_getter_enum_tester, num_length_getter_enums>();
//}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(example_numbers) {
	BOOST_CHECK_EQUAL(           mean_length_getter().get_prot_only_length( dummy_protein_10, dummy_protein_90 ), 50 );
	BOOST_CHECK_EQUAL(           mean_length_getter().get_prot_only_length( dummy_protein_90, dummy_protein_10 ), 50 );

	BOOST_CHECK_EQUAL( geometric_mean_length_getter().get_prot_only_length( dummy_protein_10, dummy_protein_90 ), 30 );
	BOOST_CHECK_EQUAL( geometric_mean_length_getter().get_prot_only_length( dummy_protein_90, dummy_protein_10 ), 30 );

	BOOST_CHECK_EQUAL(      length_of_longer_getter().get_prot_only_length( dummy_protein_10, dummy_protein_90 ), 90 );
	BOOST_CHECK_EQUAL(      length_of_longer_getter().get_prot_only_length( dummy_protein_90, dummy_protein_10 ), 90 );

	BOOST_CHECK_EQUAL(     length_of_shorter_getter().get_prot_only_length( dummy_protein_10, dummy_protein_90 ), 10 );
	BOOST_CHECK_EQUAL(     length_of_shorter_getter().get_prot_only_length( dummy_protein_90, dummy_protein_10 ), 10 );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE_TEMPLATE(names_work, length_getter_type, length_getter_types) {
	const auto the_getter = length_getter_type{};

	BOOST_CHECK_GT( the_getter.human_friendly_short_name().length(),   0 );
	BOOST_CHECK_GT( the_getter.id_name().length(),                     0 );
	BOOST_CHECK_GT( the_getter.long_name().length(),                   0 );
	BOOST_CHECK_GT( the_getter.description().length(),                 0 );

	constexpr length_getter_enum my_length_getter_enum = enum_of_length_getter_impl<length_getter_type>::value;
	using orig_type = length_getter_of_enum_t<my_length_getter_enum>;
	static_assert( is_same<length_getter_type, orig_type>::value, "Length_getter enum round trip produced different type" );
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE_TEMPLATE(short_name_suffix_works, length_getter_type, length_getter_types) {
	const auto short_name_suffixes = length_getter_as_short_name_suffixes( length_getter_type(), { length_getter_category::LONGER, length_getter_category::SHORTER } );
	BOOST_CHECK_GT( short_name_suffixes.size(), 0 );
}

BOOST_AUTO_TEST_SUITE_END()
