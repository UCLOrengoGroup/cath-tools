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

#include <boost/test/auto_unit_test.hpp>

#include "common/boost_addenda/range/indices.hpp"
#include "common/exception/invalid_argument_exception.hpp"
#include "common/type_aliases.hpp"
#include "structure/protein/amino_acid.hpp"
#include "test/boost_addenda/boost_check_no_throw_diag.hpp"

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace std;

/// \brief The sequence of 1cjn0, 1cjo0, 1roeA, 2cjnA and 2cjoA (prepended with HMWX) in names
const str_vec  NAMES   = { "Histidine", "Methionine", "Tryptophan", "Unknown", "Alanine", "Threonine", "Tyrosine", "Lysine", "Valine", "Threonine", "Leucine", "Valine", "Arginine", "Proline", "Aspartic Acid", "Glycine", "Serine", "Glutamic Acid", "Threonine", "Threonine", "Isoleucine", "Aspartic Acid", "Valine", "Proline", "Glutamic Acid", "Aspartic Acid", "Glutamic Acid", "Tyrosine", "Isoleucine", "Leucine", "Aspartic Acid", "Valine", "Alanine", "Glutamic Acid", "Glutamic Acid", "Glutamine", "Glycine", "Leucine", "Aspartic Acid", "Leucine", "Proline", "Phenylalanine", "Serine", "Cysteine", "Arginine", "Alanine", "Glycine", "Alanine", "Cysteine", "Serine", "Threonine", "Cysteine", "Alanine", "Glycine", "Lysine", "Leucine", "Leucine", "Glutamic Acid", "Glycine", "Glutamic Acid", "Valine", "Aspartic Acid", "Glutamine", "Serine", "Aspartic Acid", "Glutamine", "Serine", "Phenylalanine", "Leucine", "Aspartic Acid", "Aspartic Acid", "Aspartic Acid", "Glutamine", "Isoleucine", "Glutamic Acid", "Lysine", "Glycine", "Phenylalanine", "Valine", "Leucine", "Threonine", "Cysteine", "Valine", "Alanine", "Tyrosine", "Proline", "Arginine", "Serine", "Aspartic Acid", "Cysteine", "Lysine", "Isoleucine", "Leucine", "Threonine", "Asparagine", "Glutamine", "Glutamic Acid", "Glutamic Acid", "Glutamic Acid", "Leucine", "Tyrosine" };

/// \brief The sequence of 1cjn0, 1cjo0, 1roeA, 2cjnA and 2cjoA (prepended with HMWX) in codes
const str_vec  CODES   = { "HIS", "MET", "TRP", "UNK", "ALA", "THR", "TYR", "LYS", "VAL", "THR", "LEU", "VAL", "ARG", "PRO", "ASP", "GLY", "SER", "GLU", "THR", "THR", "ILE", "ASP", "VAL", "PRO", "GLU", "ASP", "GLU", "TYR", "ILE", "LEU", "ASP", "VAL", "ALA", "GLU", "GLU", "GLN", "GLY", "LEU", "ASP", "LEU", "PRO", "PHE", "SER", "CYS", "ARG", "ALA", "GLY", "ALA", "CYS", "SER", "THR", "CYS", "ALA", "GLY", "LYS", "LEU", "LEU", "GLU", "GLY", "GLU", "VAL", "ASP", "GLN", "SER", "ASP", "GLN", "SER", "PHE", "LEU", "ASP", "ASP", "ASP", "GLN", "ILE", "GLU", "LYS", "GLY", "PHE", "VAL", "LEU", "THR", "CYS", "VAL", "ALA", "TYR", "PRO", "ARG", "SER", "ASP", "CYS", "LYS", "ILE", "LEU", "THR", "ASN", "GLN", "GLU", "GLU", "GLU", "LEU", "TYR" };

/// \brief The sequence of 1cjn0, 1cjo0, 1roeA, 2cjnA and 2cjoA (prepended with HMWX) in letters
const char_vec LETTERS = { 'H',   'M',   'W',   'X',   'A',   'T',   'Y',   'K',   'V',   'T',   'L',   'V',   'R',   'P',   'D',   'G',   'S',   'E',   'T',   'T',   'I',   'D',   'V',   'P',   'E',   'D',   'E',   'Y',   'I',   'L',   'D',   'V',   'A',   'E',   'E',   'Q',   'G',   'L',   'D',   'L',   'P',   'F',   'S',   'C',   'R',   'A',   'G',   'A',   'C',   'S',   'T',   'C',   'A',   'G',   'K',   'L',   'L',   'E',   'G',   'E',   'V',   'D',   'Q',   'S',   'D',   'Q',   'S',   'F',   'L',   'D',   'D',   'D',   'Q',   'I',   'E',   'K',   'G',   'F',   'V',   'L',   'T',   'C',   'V',   'A',   'Y',   'P',   'R',   'S',   'D',   'C',   'K',   'I',   'L',   'T',   'N',   'Q',   'E',   'E',   'E',   'L',   'Y'   };

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

