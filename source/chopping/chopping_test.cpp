/// \file
/// \brief The chopping test suite

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

#include "chopping/chopping.hpp"
#include "chopping/region/region.hpp"
#include "test/global_test_constants.hpp"

using namespace cath;
using namespace cath::chop;

// Chopping types: domall, new_to_jmol_selection, new_to_string, scop, simple

// Update:
// "2p16 D16-27[A]+121-232[A] D28-120[A]+233-244[A]"
// "1hvc D1(B)-99(B)[A] D1(A)-99(A)[A] F200-204[A]"
//
// Domall
// "2p16A D02 F00  2  A   16 - A   27 -  A  121 - A  232 -  2  A   28 - A  120 -  A  233 - A  244 -"
// "1hvcA D02 F01  1  A    1 B A   99 B  1  A    1 A A   99 A  A  200 - A  204 - (5)"
//
// Domain seqres chopping:
// 2p16A01 : "1-12,107-222"
// 2p16A02 : "13-106,223-234"
// 1hvcA01 : "1-99"
// 1hvcA02 : "105-203"
//
// SCOP insert codes:
// d1jqga2 1jqg    A:4P-100P       d.58.3.1        71792   cl=53931,cf=54861,sf=54897,fa=54898,dm=54899,sp=75429,px=71792
// d1kwma2 1kwm    A:1A-95A        d.58.3.1        73080   cl=53931,cf=54861,sf=54897,fa=54898,dm=54903,sp=75430,px=73080
// d1kwmb2 1kwm    B:1A-95A        d.58.3.1        73082   cl=53931,cf=54861,sf=54897,fa=54898,dm=54903,sp=75430,px=73082
// d1pcaa1 1pca    A:4A-99A        d.58.3.1        39063   cl=53931,cf=54861,sf=54897,fa=54898,dm=54899,sp=54900,px=39063
// d1qdma1 1qdm    A:1S-104S       a.64.1.2        18141   cl=46456,cf=47861,sf=47862,fa=47866,dm=47867,sp=47868,px=18141
// d1qdmb1 1qdm    B:1S-104S       a.64.1.2        18142   cl=46456,cf=47861,sf=47862,fa=47866,dm=47867,sp=47868,px=18142
// d1qdmc1 1qdm    C:1S-104S       a.64.1.2        18143   cl=46456,cf=47861,sf=47862,fa=47866,dm=47867,sp=47868,px=18143
//
// SCOP multi-chain-domains:
// d1a0h.1 1a0h    A:271-320,B:    b.47.1.2        26229   cl=48724,cf=50493,sf=50494,fa=50514,dm=50531,sp=50533,px=26229
// d1a0h.2 1a0h    D:271-320,E:    b.47.1.2        26230   cl=48724,cf=50493,sf=50494,fa=50514,dm=50531,sp=50533,px=26230
// d1bi6.2 1bi6    L:,H:1-7,H:32-41        g.3.12.1        44341   cl=56992,cf=57015,sf=57243,fa=57244,dm=57245,sp=57246,px=44341
// d1bpl.1 1bpl    A:,B:193-393    c.1.8.1 28703   cl=51349,cf=51350,sf=51445,fa=51446,dm=51447,sp=51448,px=28703
// d1cl7.1 1cl7    H:114-128,I:    b.1.1.2 21389   cl=48724,cf=48725,sf=48726,fa=48942,dm=88574,sp=88576,px=21389
// d1dan.1 1dan    T:,U:91-106     b.1.2.1 21953   cl=48724,cf=48725,sf=49265,fa=49266,dm=49267,sp=49268,px=21953
// d1qrj.1 1qrj    A:,B:16-130     a.73.1.1        18332   cl=46456,cf=47942,sf=47943,fa=47944,dm=47949,sp=47950,px=18332
// d2bi6.2 2bi6    L:,H:1-7,H:32-41        g.3.12.1        44339   cl=56992,cf=57015,sf=57243,fa=57244,dm=57245,sp=57246,px=44339
// d2pjr.1 2pjr    A:319-548,B:    c.37.1.19       32396   cl=51349,cf=52539,sf=52540,fa=81268,dm=52701,sp=52702,px=32396
// d2pjr.2 2pjr    F:1019-1247,G:  c.37.1.19       32398   cl=51349,cf=52539,sf=52540,fa=81268,dm=52701,sp=52702,px=32398

namespace cath {
	namespace test {

		/// \brief The chopping_test_suite_fixture to assist in testing chopping
		struct chopping_test_suite_fixture: protected global_test_constants {
		protected:
			~chopping_test_suite_fixture() noexcept = default;
		};

	}
}

/// \brief TODOCUMENT
BOOST_FIXTURE_TEST_SUITE(chopping_test_suite, cath::test::chopping_test_suite_fixture)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(basic) {
//	2p16A01:
//	new_to_string '2p16A01' '2p16 D16-27[A]+121-232[A]',
//	simple        '2p16A01' '2p16 16-27[A]+121-232[A]',
//	new_to_string '1hvcA01' '1hvc D1B-99B[A]'
//	simple        '1hvcA01' '1hvc 1(B)-99(B)[A]'

	BOOST_CHECK( true );
}

BOOST_AUTO_TEST_SUITE_END()
