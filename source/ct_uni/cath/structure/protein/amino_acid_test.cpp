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

#include <array>
#include <optional>
#include <string>
#include <string_view>

#include <boost/test/unit_test.hpp>

#include "cath/common/boost_addenda/range/indices.hpp"
#include "cath/common/cpp20/make_array.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/string/string_view_of_char_arr.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/structure/protein/amino_acid.hpp"
#include "cath/test/boost_addenda/boost_check_no_throw_diag.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::file;

using ::std::array;
using ::std::literals::string_view_literals::operator""sv;
using ::std::make_optional;
using ::std::nullopt;
using ::std::string;
using ::std::string_view;

/// \brief The sequence of 1cjn0, 1cjo0, 1roeA, 2cjnA and 2cjoA (prepended with HMWX) in names
constexpr     auto NAMES   = make_array( "Histidine"sv, "Methionine"sv, "Tryptophan"sv, "Unknown"sv, "Alanine"sv, "Threonine"sv, "Tyrosine"sv, "Lysine"sv, "Valine"sv, "Threonine"sv, "Leucine"sv, "Valine"sv, "Arginine"sv, "Proline"sv, "Aspartic Acid"sv, "Glycine"sv, "Serine"sv, "Glutamic Acid"sv, "Threonine"sv, "Threonine"sv, "Isoleucine"sv, "Aspartic Acid"sv, "Valine"sv, "Proline"sv, "Glutamic Acid"sv, "Aspartic Acid"sv, "Glutamic Acid"sv, "Tyrosine"sv, "Isoleucine"sv, "Leucine"sv, "Aspartic Acid"sv, "Valine"sv, "Alanine"sv, "Glutamic Acid"sv, "Glutamic Acid"sv, "Glutamine"sv, "Glycine"sv, "Leucine"sv, "Aspartic Acid"sv, "Leucine"sv, "Proline"sv, "Phenylalanine"sv, "Serine"sv, "Cysteine"sv, "Arginine"sv, "Alanine"sv, "Glycine"sv, "Alanine"sv, "Cysteine"sv, "Serine"sv, "Threonine"sv, "Cysteine"sv, "Alanine"sv, "Glycine"sv, "Lysine"sv, "Leucine"sv, "Leucine"sv, "Glutamic Acid"sv, "Glycine"sv, "Glutamic Acid"sv, "Valine"sv, "Aspartic Acid"sv, "Glutamine"sv, "Serine"sv, "Aspartic Acid"sv, "Glutamine"sv, "Serine"sv, "Phenylalanine"sv, "Leucine"sv, "Aspartic Acid"sv, "Aspartic Acid"sv, "Aspartic Acid"sv, "Glutamine"sv, "Isoleucine"sv, "Glutamic Acid"sv, "Lysine"sv, "Glycine"sv, "Phenylalanine"sv, "Valine"sv, "Leucine"sv, "Threonine"sv, "Cysteine"sv, "Valine"sv, "Alanine"sv, "Tyrosine"sv, "Proline"sv, "Arginine"sv, "Serine"sv, "Aspartic Acid"sv, "Cysteine"sv, "Lysine"sv, "Isoleucine"sv, "Leucine"sv, "Threonine"sv, "Asparagine"sv, "Glutamine"sv, "Glutamic Acid"sv, "Glutamic Acid"sv, "Glutamic Acid"sv, "Leucine"sv, "Tyrosine"sv );

/// \brief The sequence of 1cjn0, 1cjo0, 1roeA, 2cjnA and 2cjoA (prepended with HMWX) in codes
constexpr     auto CODES   = make_array( "HIS"sv, "MET"sv, "TRP"sv, "UNK"sv, "ALA"sv, "THR"sv, "TYR"sv, "LYS"sv, "VAL"sv, "THR"sv, "LEU"sv, "VAL"sv, "ARG"sv, "PRO"sv, "ASP"sv, "GLY"sv, "SER"sv, "GLU"sv, "THR"sv, "THR"sv, "ILE"sv, "ASP"sv, "VAL"sv, "PRO"sv, "GLU"sv, "ASP"sv, "GLU"sv, "TYR"sv, "ILE"sv, "LEU"sv, "ASP"sv, "VAL"sv, "ALA"sv, "GLU"sv, "GLU"sv, "GLN"sv, "GLY"sv, "LEU"sv, "ASP"sv, "LEU"sv, "PRO"sv, "PHE"sv, "SER"sv, "CYS"sv, "ARG"sv, "ALA"sv, "GLY"sv, "ALA"sv, "CYS"sv, "SER"sv, "THR"sv, "CYS"sv, "ALA"sv, "GLY"sv, "LYS"sv, "LEU"sv, "LEU"sv, "GLU"sv, "GLY"sv, "GLU"sv, "VAL"sv, "ASP"sv, "GLN"sv, "SER"sv, "ASP"sv, "GLN"sv, "SER"sv, "PHE"sv, "LEU"sv, "ASP"sv, "ASP"sv, "ASP"sv, "GLN"sv, "ILE"sv, "GLU"sv, "LYS"sv, "GLY"sv, "PHE"sv, "VAL"sv, "LEU"sv, "THR"sv, "CYS"sv, "VAL"sv, "ALA"sv, "TYR"sv, "PRO"sv, "ARG"sv, "SER"sv, "ASP"sv, "CYS"sv, "LYS"sv, "ILE"sv, "LEU"sv, "THR"sv, "ASN"sv, "GLN"sv, "GLU"sv, "GLU"sv, "GLU"sv, "LEU"sv, "TYR"sv );

/// \brief The sequence of 1cjn0, 1cjo0, 1roeA, 2cjnA and 2cjoA (prepended with HMWX) in letters
constexpr auto LETTERS = make_array( 'H',   'M',   'W',   'X',   'A',   'T',   'Y',   'K',   'V',   'T',   'L',   'V',   'R',   'P',   'D',   'G',   'S',   'E',   'T',   'T',   'I',   'D',   'V',   'P',   'E',   'D',   'E',   'Y',   'I',   'L',   'D',   'V',   'A',   'E',   'E',   'Q',   'G',   'L',   'D',   'L',   'P',   'F',   'S',   'C',   'R',   'A',   'G',   'A',   'C',   'S',   'T',   'C',   'A',   'G',   'K',   'L',   'L',   'E',   'G',   'E',   'V',   'D',   'Q',   'S',   'D',   'Q',   'S',   'F',   'L',   'D',   'D',   'D',   'Q',   'I',   'E',   'K',   'G',   'F',   'V',   'L',   'T',   'C',   'V',   'A',   'Y',   'P',   'R',   'S',   'D',   'C',   'K',   'I',   'L',   'T',   'N',   'Q',   'E',   'E',   'E',   'L',   'Y'   );

BOOST_AUTO_TEST_SUITE(amino_acid_test_suite)

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE( static_tests ) {
	static_assert( dna_atom_of_code( make_char_arr( "  T" ) ) == make_optional( dna_atom::T ) );
	static_assert( dna_atom_of_code( make_char_arr( " +C" ) ) == make_optional( dna_atom::PLUS_C ) );
	static_assert( dna_atom_of_code( make_char_arr( " DI" ) ) == make_optional( dna_atom::DI ) );
	static_assert( dna_atom_of_code( make_char_arr( " G " ) ) == make_optional( dna_atom::G_SPACE ) );
	static_assert( dna_atom_of_code( make_char_arr( "AAA" ) ) == nullopt );

	static_assert( get_letter_code_or_name<char>( index_of_letter( 'F' ) ) == 'F' );
	static_assert( get_letter_code_or_name<char>( index_of_letter( 'L' ) ) == 'L' );
	static_assert( get_letter_code_or_name<char>( index_of_letter( 'T' ) ) == 'T' );

	static_assert( get_letter_code_or_name<char>( index_of_code( make_char_arr( "PHE" ) ) ) == 'F' );
	static_assert( get_letter_code_or_name<char>( index_of_code( make_char_arr( "LEU" ) ) ) == 'L' );
	static_assert( get_letter_code_or_name<char>( index_of_code( make_char_arr( "THR" ) ) ) == 'T' );

	static_assert( get_letter_code_or_name<char>( index_of_name( "Phenylalanine"sv ) ) == 'F' );
	static_assert( get_letter_code_or_name<char>( index_of_name( "Leucine"sv ) ) == 'L' );
	static_assert( get_letter_code_or_name<char>( index_of_name( "Threonine"sv ) ) == 'T' );

	/////

	static_assert( amino_acid( make_char_arr( "  T" ) ).get_type() == amino_acid_type::DNA );
	static_assert( amino_acid( make_char_arr( " +C" ) ).get_type() == amino_acid_type::DNA );
	static_assert( amino_acid( make_char_arr( " DI" ) ).get_type() == amino_acid_type::DNA );
	static_assert( amino_acid( make_char_arr( " G " ) ).get_type() == amino_acid_type::DNA );

	static_assert( amino_acid( make_char_arr( "PHE" ) ).get_type() == amino_acid_type::AA );
	static_assert( amino_acid( make_char_arr( "LEU" ) ).get_type() == amino_acid_type::AA );
	static_assert( amino_acid( make_char_arr( "THR" ) ).get_type() == amino_acid_type::AA );

	static_assert( amino_acid( make_char_arr( "MLY" ), pdb_record::HETATM ).get_type() == amino_acid_type::HETATOM );
	static_assert( amino_acid( make_char_arr( "MSE" ), pdb_record::HETATM ).get_type() == amino_acid_type::HETATOM );
	static_assert( amino_acid( make_char_arr( "NCX" ), pdb_record::HETATM ).get_type() == amino_acid_type::HETATOM );

	/////

	static_assert( string_view_of_char_arr( amino_acid( make_char_arr( "  T" ) ).get_code() ) == "  T"sv );
	static_assert( string_view_of_char_arr( amino_acid( make_char_arr( " +C" ) ).get_code() ) == " +C"sv );
	static_assert( string_view_of_char_arr( amino_acid( make_char_arr( " DI" ) ).get_code() ) == " DI"sv );
	static_assert( string_view_of_char_arr( amino_acid( make_char_arr( " G " ) ).get_code() ) == " G "sv );

	static_assert( string_view_of_char_arr( amino_acid( make_char_arr( "PHE" ) ).get_code() ) == "PHE"sv );
	static_assert( string_view_of_char_arr( amino_acid( make_char_arr( "LEU" ) ).get_code() ) == "LEU"sv );
	static_assert( string_view_of_char_arr( amino_acid( make_char_arr( "THR" ) ).get_code() ) == "THR"sv );

	static_assert( string_view_of_char_arr( amino_acid( make_char_arr( "MLY" ), pdb_record::HETATM ).get_code() ) == "MLY"sv );
	static_assert( string_view_of_char_arr( amino_acid( make_char_arr( "MSE" ), pdb_record::HETATM ).get_code() ) == "MSE"sv );
	static_assert( string_view_of_char_arr( amino_acid( make_char_arr( "NCX" ), pdb_record::HETATM ).get_code() ) == "NCX"sv );

	BOOST_TEST(true);
}

/// \brief TODOCUMENT
BOOST_AUTO_TEST_CASE(simple_conversion) {
	BOOST_REQUIRE_EQUAL(CODES.size(), NAMES.size());
	BOOST_REQUIRE_EQUAL(CODES.size(), LETTERS.size());
	for ( const size_t &letter_ctr : indices( CODES.size() ) ) {
		const string_view &code   = CODES.at( letter_ctr );
		const string_view &name   = NAMES.at( letter_ctr );
		const char &       letter = LETTERS.at( letter_ctr );

		// Check non-member, non-friend 1->3 and 3->1 converters
		BOOST_CHECK_EQUAL( code,   get_code_str_of_amino_acid_letter( letter ) );
		BOOST_CHECK_EQUAL( letter, get_letter_of_amino_acid_code    ( code   ) );

		// Check that the correct amino acid is constructed from any of the three labels
		const string letter_string { letter };
		const array all_names_and_codes = { code, string_view{ letter_string }, name} ;
		for (const string_view &name_or_code : all_names_and_codes) {
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

