/// \file
/// \brief The alignment_breaks test suite

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

#include "cath/alignment/alignment.hpp"
#include "cath/alignment/io/align_scaffold.hpp"
#include "cath/alignment/tools/alignment_breaks.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/test/boost_addenda/boost_check_equal_ranges.hpp"
#include "cath/test/boost_test_print_type.hpp"
//#include "cath/alignment/test/alignment_fixture.hpp"
//#include "cath/common/exception/invalid_argument_exception.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::common;
using namespace ::std;

namespace {

	/// \brief The alignment_breaks_test_suite_fixture to assist in testing alignment_breaks
	struct alignment_breaks_test_suite_fixture {
	protected:
		~alignment_breaks_test_suite_fixture() noexcept = default;

		/// \brief An example alignment for testing alignment breaks code
		const alignment the_alignment = alignment_of_scaffold_lines( {
			// 0         1          < These two lines count up from 01 to 17 (one less than the alignment length)
			//  12345678901234567   <
			//   //     /    ////   < This line connects the between-index breaks (think: fence posts between fences) with their corresponding indices
			  "XX   XXXX  XXX   X",
			  "   XXXX  XXXX X X ",
			  "  X  XX    XX  X  "
			//   \\     \    \\\\   < This line helps re-emphasise the between-index breaks
		} );

		/// \brief The correct alignment breaks of the_alignment
		const size_vec expected_breaks{  2,  3,  9, 14, 15, 16, 17 };

		/// \brief The correct alignment break pairs of the_alignment
		const size_size_pair_vec expected_break_pairs{
			{  2,  2 },
			{  2,  3 },
			{  3,  3 },
			{  9,  9 },
			{ 14, 14 },
			{ 14, 15 },
			{ 15, 15 },
			{ 16, 16 },
			{ 16, 17 },
			{ 17, 17 }
		};

		/// \brief Example FunFam 22219 in 3.40.710.10
		const alignment big_alignment = alignment_of_scaffold_lines( {
			// 0                                                                                                   1                                                                                                   2                                                                                                   3
			// 0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8         9         0         1         2         3         4         5         6         7         8
			//  123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901
			  "                         X   XXX     XX XXXXXXX X    X      XX XXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXX XX XXXXX  XXXXXXXXXXXXXX X XXXXX XXXXX XXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXX      X               ",
			  "                           XXXXX     XX XXXXXXX X    X      XX XXXXXXXXXXXXX  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXX XX XXX    XXXXXXXXX  XXX X XXXXX XXXXX XXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXX                      ",
			  "                        XX   XXX     X XXXXXXX XX    X      X XXXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXXXXX XXXX X    XXXXXXXXX  XXX X XXXX XXXXXX XXXXXXXXXXXXXXXXXXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXXX XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXX          X           ",
			  "                  XX    XX   XXX     X XXXXXXX XX    X      X XXXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXXXXX XXXX X    XXXXXXXXX  XXX X XXXX XXXXXX XXXXXXXXXXXXXXXXXXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXXX XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXX          X           ",
			  "           X              XXXXXX     XX XXXXXXX X    X      XX XXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXX XX XXX    XXXXXXXXX  XXX X XXXXX XXXXX XXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXX  XX                  ",
			  "                         X   XXX     XX XXXXXXX X    X      XX XXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXX XX XXX    XXXXXXXXX  XXX X XXXXX XXXXX XXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXX       XXX            ",
			  "                           XXXXX     XX XXXXXXX X    X      XX XXXXXXXXXXXXX  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXX  XX XXX    XXXXXXXXX  XX  X XXXXX XXXX  XXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXX  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXX                      ",
			  "                           XXXXX     XX XXXXXXX X    X      XX XXXXXXXXXXXXX X XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXX XX XXX    XXXXXXXXX  XXX X XXXXX XXXXX XXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXX            XX        ",
			  " XXXXXXXXXX                  XXX     XX XXXXXXX X    X      XX XXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXX XX XXX    XXXXXXXXX  XXX X XXXXX XXXXX XXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXXXXX XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXXX                     ",
			  "X                            XXX     XX XXXXXXX XXXXXX      XX XXXXXXXXXXXXXX  XXXXXXXXXXXXXXXXXXXXXXXXXXXX      XXXX                                XXXX XXXXX XXXXXXXXXXXXXXX XXXXXXXXXXXXXX    XXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXXXXXX   XX  XXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXXXX XXXXXXXX  X XXXXXXXXXXXXX  XXXXXXXXXXXX X    X XXXXXXXXXXXXXXXXXXX                       ",
			  "                         X   XXX     XX XXXXXXX X    X      XX XXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXX XX XXX    XXXXXXXXX  XXX X XXXXX XXXXX XXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXX  XXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXX  XXXXXXXXXXXXXXXXXXXXXX X                    ",
			  "                X        X   XXX     XX XXXXXXX X    X      XX XXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXX XX XXX    XXXXXXXXX  XXX   XXXXX XXXXX XXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XX XXXXXXXXXXXX XXXXX                       ",
			  "            XXXX             XXX     XX XXXXXXX X    X      XX XXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXX XX XXX    XXXXXXXXX  XXX XXXXXXX XXXXX XXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXX X  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXX    XX                ",
			  "                            XXXX     XX XXXXXXX X    XXXXXXXXX XXXXXXXXXXXXX  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXX XX XXX    XXXXXXXXX  XXX X XXXXX XXXXX XXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXX              XX      ",
			  "                                     XX XXXXXXX X    X      XX XXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXX  XXXX XXXXXXX XX XXX    XXXXXXXXX  XXX X XXXXX XXXXX XXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXX                      ",
			  "                         X   XXX     XX XXXXXXX X    X      XX XXXXX  XXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXX X XX XXX    XX XXXXXX  XXX X XXXXX XXXXX XXXXXXXXXXXXXXX XXX               XX  XXXXX XXXXX   XXXXXXXXXXX XXXXXXXXXXXXXXX XXXX XXXXX XXXX   XXXXXXXX        XXXXX XXXXXXXX XXXXXXXXXX X  XXXXXXXXX XXXXX  XXXXXXXXXXXX X    X XXXXXX  XXXXXXXXXXXX                      ",
			  "                 X       X   XXX     XX XXXXXXX X    X      XX XXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXX XX XXX    XXXXXXXXX  XXX X XXXXX XXXXX XXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXX                      ",
			  "                         X   XXX     XX XXXXXXX X    X      XX XXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXX XX XXX    XXXXXXXXX  XXX X XXXXX XXXXX XXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXX                      ",
			  "                        XX   XXX     XX XXXXXXX X    X      XX XXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXX XX XXX    XXXXXXXXX  XXX X XXXXX XXXXX XXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXXX XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXX          X           ",
			  "                                XX   XX XXXXXXX X    X      XX XXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXX  XXXX XXXXXXX XX XXX    XXXXXXXXX  XXX X XXXXX XXXXX XXXXXXXXXXXXXX  XXX               XXXXXX  XXXXXX   XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXX XXXX XXXXXXX XXXXXXX        XXXXXXXXXXXX X XXXXXXXXXXXX  XXXXXXXXXXXXXXX XXXXXXXXX XXX X    XXXXXXXXXXXXXXXXXXXXXX                 XX   ",
			  "                                  XXXXX XXXXXXX X    X      XX XXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXX  XXXX XXXXXXX XX XXX    XXXXXXXXX  XXX X XXXXX XXXXX XXXXXXXXXXXXXXX XXX           XXXXXXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXXXX  XXXXXXXXXXXX X  XXXXXXXXXXXXXXXXXXXXXXXX                   X  ",
			  "                             XXX     XX XXXXXXX X    X      XX XXXXXXXX    X   XXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXX XXXXXXX XX XXX  XXXXXXXXXXX  XXX X XXXXX XXXXX XXX XXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXX X XXXXXXXXXXXX XXXXXXXXXXXXXXXX  XXXXXXXXXXXXXX    XXXXXXXXXXXXXXXXXXXXXX       XXX          XX",
			  "                              XX     XX XXXXXXX X    X      XX XXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXX  XXXX XXXXXXX XX XXX    XXXXXXXXX  XXX X XXXXX XXXXX XXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXX                X     ",
			  "                                                X    X      XX XXXXXXXXXXXXX  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXX XXX XXX    XXXXXXXXX  XX XX XXXXX XXXX XXXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXX                       ",
			  "                          XXXXXX     XX XXXXXXX X    X      XX XXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXX XX XXX    XXXXXXXXX  XXX X XXXXX XXXXX XXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXX                      ",
			  "                             XXX     XX XXXXXXX X    X      XX XXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XXXXXXX XX XXX    XXXXXXXXX  XXX X XXXXX XXXXX XXXXXXXXXXXXXXX XXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXX  XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXX                      ",
			  "                    XXXXXX   XXX     X  XXXXXX  X    X      X  XXXXXXXXXXXXX   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX XX XXXXXXX XX X X    XXXXXXXXX  XXX X XXXX  XXXXX XXXXXXXXXXXXXXXXXXX               XXXXXXXXXXXXXXX  XXXXXXXXXXXX XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX        XXXXXXXXXXXXXX XXXXXXXXXXXXX XXXXXXXXXXXXXXX  XXXXXXXXXXXX X    XXXXXXXXXXXXXXXXXXXXXX          XX          "
		} );
	};

} // namespace

BOOST_FIXTURE_TEST_SUITE(alignment_breaks_test_suite, alignment_breaks_test_suite_fixture)

BOOST_AUTO_TEST_CASE(alignment_breaks_are_correct_for_simple_alignment) {
	BOOST_CHECK_EQUAL_RANGES( get_alignment_breaks( the_alignment ), expected_breaks );
}

BOOST_AUTO_TEST_CASE(alignment_break_pairs_are_correct_for_simple_alignment) {
	BOOST_TEST( get_alignment_break_pairs( the_alignment ) == expected_break_pairs );
}

BOOST_AUTO_TEST_CASE(alignment_break_pairs_are_correct_for_big_alignment) {
	BOOST_CHECK_EQUAL  ( big_alignment.length(), 382 );
	BOOST_REQUIRE_EQUAL( get_alignment_break_pairs( big_alignment ).size(), 121_z );
	BOOST_TEST         ( get_alignment_break_pairs( big_alignment ).front() == make_pair(   1_z,   1_z ) );
	BOOST_TEST         ( get_alignment_break_pairs( big_alignment ).back()  == make_pair( 380_z, 380_z ) );
}


BOOST_AUTO_TEST_SUITE_END()
