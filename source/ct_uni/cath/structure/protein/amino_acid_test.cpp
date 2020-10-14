/// \file
/// \brief The amino_acid test suite

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

#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/cpp20/make_array.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/structure/protein/amino_acid.hpp"
#include "cath/test/boost_addenda/boost_check_no_throw_diag.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace std;

/// \brief The sequence of 1cjn0, 1cjo0, 1roeA, 2cjnA and 2cjoA (prepended with HMWX) in names
const     auto NAMES   = make_array( "Histidine"s, "Methionine"s, "Tryptophan"s, "Unknown"s, "Alanine"s, "Threonine"s, "Tyrosine"s, "Lysine"s, "Valine"s, "Threonine"s, "Leucine"s, "Valine"s, "Arginine"s, "Proline"s, "Aspartic Acid"s, "Glycine"s, "Serine"s, "Glutamic Acid"s, "Threonine"s, "Threonine"s, "Isoleucine"s, "Aspartic Acid"s, "Valine"s, "Proline"s, "Glutamic Acid"s, "Aspartic Acid"s, "Glutamic Acid"s, "Tyrosine"s, "Isoleucine"s, "Leucine"s, "Aspartic Acid"s, "Valine"s, "Alanine"s, "Glutamic Acid"s, "Glutamic Acid"s, "Glutamine"s, "Glycine"s, "Leucine"s, "Aspartic Acid"s, "Leucine"s, "Proline"s, "Phenylalanine"s, "Serine"s, "Cysteine"s, "Arginine"s, "Alanine"s, "Glycine"s, "Alanine"s, "Cysteine"s, "Serine"s, "Threonine"s, "Cysteine"s, "Alanine"s, "Glycine"s, "Lysine"s, "Leucine"s, "Leucine"s, "Glutamic Acid"s, "Glycine"s, "Glutamic Acid"s, "Valine"s, "Aspartic Acid"s, "Glutamine"s, "Serine"s, "Aspartic Acid"s, "Glutamine"s, "Serine"s, "Phenylalanine"s, "Leucine"s, "Aspartic Acid"s, "Aspartic Acid"s, "Aspartic Acid"s, "Glutamine"s, "Isoleucine"s, "Glutamic Acid"s, "Lysine"s, "Glycine"s, "Phenylalanine"s, "Valine"s, "Leucine"s, "Threonine"s, "Cysteine"s, "Valine"s, "Alanine"s, "Tyrosine"s, "Proline"s, "Arginine"s, "Serine"s, "Aspartic Acid"s, "Cysteine"s, "Lysine"s, "Isoleucine"s, "Leucine"s, "Threonine"s, "Asparagine"s, "Glutamine"s, "Glutamic Acid"s, "Glutamic Acid"s, "Glutamic Acid"s, "Leucine"s, "Tyrosine"s );

/// \brief The sequence of 1cjn0, 1cjo0, 1roeA, 2cjnA and 2cjoA (prepended with HMWX) in codes
const     auto CODES   = make_array( "HIS"s, "MET"s, "TRP"s, "UNK"s, "ALA"s, "THR"s, "TYR"s, "LYS"s, "VAL"s, "THR"s, "LEU"s, "VAL"s, "ARG"s, "PRO"s, "ASP"s, "GLY"s, "SER"s, "GLU"s, "THR"s, "THR"s, "ILE"s, "ASP"s, "VAL"s, "PRO"s, "GLU"s, "ASP"s, "GLU"s, "TYR"s, "ILE"s, "LEU"s, "ASP"s, "VAL"s, "ALA"s, "GLU"s, "GLU"s, "GLN"s, "GLY"s, "LEU"s, "ASP"s, "LEU"s, "PRO"s, "PHE"s, "SER"s, "CYS"s, "ARG"s, "ALA"s, "GLY"s, "ALA"s, "CYS"s, "SER"s, "THR"s, "CYS"s, "ALA"s, "GLY"s, "LYS"s, "LEU"s, "LEU"s, "GLU"s, "GLY"s, "GLU"s, "VAL"s, "ASP"s, "GLN"s, "SER"s, "ASP"s, "GLN"s, "SER"s, "PHE"s, "LEU"s, "ASP"s, "ASP"s, "ASP"s, "GLN"s, "ILE"s, "GLU"s, "LYS"s, "GLY"s, "PHE"s, "VAL"s, "LEU"s, "THR"s, "CYS"s, "VAL"s, "ALA"s, "TYR"s, "PRO"s, "ARG"s, "SER"s, "ASP"s, "CYS"s, "LYS"s, "ILE"s, "LEU"s, "THR"s, "ASN"s, "GLN"s, "GLU"s, "GLU"s, "GLU"s, "LEU"s, "TYR"s );

/// \brief The sequence of 1cjn0, 1cjo0, 1roeA, 2cjnA and 2cjoA (prepended with HMWX) in letters
constexpr auto LETTERS = make_array( 'H',   'M',   'W',   'X',   'A',   'T',   'Y',   'K',   'V',   'T',   'L',   'V',   'R',   'P',   'D',   'G',   'S',   'E',   'T',   'T',   'I',   'D',   'V',   'P',   'E',   'D',   'E',   'Y',   'I',   'L',   'D',   'V',   'A',   'E',   'E',   'Q',   'G',   'L',   'D',   'L',   'P',   'F',   'S',   'C',   'R',   'A',   'G',   'A',   'C',   'S',   'T',   'C',   'A',   'G',   'K',   'L',   'L',   'E',   'G',   'E',   'V',   'D',   'Q',   'S',   'D',   'Q',   'S',   'F',   'L',   'D',   'D',   'D',   'Q',   'I',   'E',   'K',   'G',   'F',   'V',   'L',   'T',   'C',   'V',   'A',   'Y',   'P',   'R',   'S',   'D',   'C',   'K',   'I',   'L',   'T',   'N',   'Q',   'E',   'E',   'E',   'L',   'Y'   );

BOOST_AUTO_TEST_SUITE(amino_acid_test_suite)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(simple_conversion) {
	BOOST_REQUIRE_EQUAL(CODES.size(), NAMES.size());
	BOOST_REQUIRE_EQUAL(CODES.size(), LETTERS.size());
	for (const size_t &letter_ctr : indices( CODES.size() ) ) {
		const string &code   = CODES[letter_ctr];
		const string &name   = NAMES[letter_ctr];
		const char   &letter = LETTERS[letter_ctr];

		// Check non-member, non-friend 1->3 and 3->1 converters
		BOOST_CHECK_EQUAL( code,   get_code_str_of_amino_acid_letter( letter ) );
		BOOST_CHECK_EQUAL( letter, get_letter_of_amino_acid_code    ( code   ) );

		// Check that the correct amino acid is constructed from any of the three labels
		const str_vec all_names_and_codes = { code, string{ letter }, name} ;
		for (const string &name_or_code : all_names_and_codes) {
			const amino_acid the_amino_acid(name_or_code);
			BOOST_CHECK_EQUAL( the_amino_acid.get_name(),                  name   );
			BOOST_CHECK_EQUAL( *the_amino_acid.get_letter_if_amino_acid(), letter );
			BOOST_CHECK_EQUAL( get_code_string( the_amino_acid ),          code   );
		}
	}
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(throws_on_invalid_conversion) {
	const str_str_pair_vec valid_and_invalid_pairs = {
		{ "D",             "#"                  },
		{ "ASP",           "BOB"                },
		{ "Aspartic Acid", "yessiree_billy_bob" }
	};

	for (const str_str_pair &valid_and_invalid : valid_and_invalid_pairs) {
		BOOST_CHECK_THROW        ( amino_acid the_amino_acid( valid_and_invalid.second ), invalid_argument_exception );
		BOOST_CHECK_NO_THROW_DIAG( amino_acid the_amino_acid( valid_and_invalid.first  )                             );
	}
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(dna_and_rna_pseudo_residue_names_do_not_throw) {
	const str_vec dna_and_rna_pseudo_residue_names = { "  A", "  C", "  G", "  N", "  U", " DA", " DC", " DG", " DT", " DU" };
	for (const string &dna_and_rna_pseudo_residue_name : dna_and_rna_pseudo_residue_names) {
		BOOST_CHECK_NO_THROW_DIAG( amino_acid the_amino_acid(dna_and_rna_pseudo_residue_name) );
	}
}

BOOST_AUTO_TEST_CASE(less_than_works) {
	// AOP in 1b3a
	// MSP in 1pfy
	// PCP in 1ap5
	// SIS in 4f8u
	// TIY in 3tiy
	BOOST_CHECK_LT   ( amino_acid( 'F'                       ), amino_acid( 'G'                       ) );
	BOOST_CHECK_LT   ( amino_acid( 'G'                       ), amino_acid( "AOP", pdb_record::HETATM ) );
	BOOST_CHECK_LT   ( amino_acid( "AOP", pdb_record::HETATM ), amino_acid( "MSP", pdb_record::HETATM ) );
	BOOST_CHECK_LT   ( amino_acid( "MSP", pdb_record::HETATM ), amino_acid( "PCP", pdb_record::HETATM ) );
	BOOST_CHECK_LT   ( amino_acid( "PCP", pdb_record::HETATM ), amino_acid( "SIS", pdb_record::HETATM ) );
	BOOST_CHECK_LT   ( amino_acid( "SIS", pdb_record::HETATM ), amino_acid( "TIY", pdb_record::HETATM ) );

	BOOST_CHECK_EQUAL( amino_acid( "GLU", pdb_record::HETATM ), amino_acid( "GLU"                     ) );

	BOOST_CHECK_EQUAL( amino_acid( 'G'                       ), amino_acid( 'G'                       ) );
	BOOST_CHECK_EQUAL( amino_acid( "MSP", pdb_record::HETATM ), amino_acid( "MSP", pdb_record::HETATM ) );
}

BOOST_AUTO_TEST_CASE(is_water_works) {
	BOOST_CHECK(   is_water( amino_acid( "HOH", pdb_record::HETATM ) ) );
	BOOST_CHECK( ! is_water( amino_acid( "MSP", pdb_record::HETATM ) ) );
}

BOOST_AUTO_TEST_SUITE_END()

