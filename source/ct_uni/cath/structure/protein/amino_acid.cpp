/// \file
/// \brief The amino_acid class definitions

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

#include "amino_acid.hpp"

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "cath/common/algorithm/contains.hpp"
#include "cath/common/algorithm/transform_build.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"

#include <set>

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::std;

constexpr size_t amino_acid::NUM_HETATM_CHARS;

/// \brief TODOCUMENT
///
/// Examples:
///  * PDBs 1dgo, 1ii1, 1js7, 1qe7, 3ga6 and 4fm0 all contain
///    ATOM records with each of the residue names: DA, DC, DG, DT and DU
///  * PDBs 1eg0, 1fjf, 2i82 and 3u5f all contain
///    ATOM records with each of the residue names: A, C, G, N and U
const map<string, dna_atom> DNA_RNA_RESIDUE_NAMES = {
	{ "  A", dna_atom::A       },
	{ "  C", dna_atom::C       },
	{ "  G", dna_atom::G       },
	{ "  I", dna_atom::I       }, // eg in 1xnr
	{ "  N", dna_atom::N       },
	{ "  T", dna_atom::T       }, // eg in 3dpv
	{ "  U", dna_atom::U       },

	{ " DA", dna_atom::DA      },
	{ " DC", dna_atom::DC      },
	{ " DG", dna_atom::DG      },
	{ " DI", dna_atom::DI      }, // eg in 3rzl
	{ " DT", dna_atom::DT      },
	{ " DU", dna_atom::DU      },

	{ " +A", dna_atom::PLUS_A  }, // eg in 1hp6
	{ " +C", dna_atom::PLUS_C  }, // eg in 356d
	{ " +G", dna_atom::PLUS_G  }, // eg in 1gpg
	{ " +U", dna_atom::PLUS_U  }, // eg in 1hp6

	{ " A ", dna_atom::A_SPACE }, // eg in 3zvp
	{ " C ", dna_atom::C_SPACE }, // eg in 3zvp
	{ " G ", dna_atom::G_SPACE }, // eg in 4a1d
	{ " U ", dna_atom::U_SPACE }  // eg in 4a1d
};

/// \brief TODOCUMENT
//const string amino_acid::UNKNOWN_AMINO_ACID_NAME("Unknown");

///// \brief TODOCUMENT
//const string & global_test_constants::EXAMPLE_B_PDB_STEMNAME() {
//	static const string example_b_pdb_stemname( "1hdoA00" );
//	return example_b_pdb_stemname;
//}

/// \brief Generate a string describing the specified amino_acid_type
///
/// \relates amino_acid_type
string cath::to_string(const amino_acid_type &prm_amino_acid_type ///< The amino_acid_type to describe
                       ) {
	switch ( prm_amino_acid_type ) {
		case ( amino_acid_type::AA      ) : { return "AA"      ; };
		case ( amino_acid_type::HETATOM ) : { return "HETATOM" ; };
		case ( amino_acid_type::DNA     ) : { return "DNA"     ; };
	}
	BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Value of amino_acid_type not recognised whilst converting to_string()"));
}

/// \brief TODOCUMENT
auto amino_acid::INDEX_OF_LETTER() -> const char_size_unordered_map & {
	static const char_size_unordered_map index_of_letter(amino_acid::build_index_unordered_map<char,   0>());
	return index_of_letter;
}

/// \brief TODOCUMENT
auto amino_acid::INDEX_OF_CODE() -> const string_size_unordered_map & {
	static const string_size_unordered_map index_of_code(amino_acid::build_index_unordered_map<string, 1>(
		[] (const char_3_arr &x) { return char_arr_to_string( x ); }
	));
	return index_of_code;
}

/// \brief TODOCUMENT
auto amino_acid::INDEX_OF_NAME() -> const string_size_unordered_map & {
	static const string_size_unordered_map index_of_name(amino_acid::build_index_unordered_map<string, 2>());
	return index_of_name;
}

/// \brief TODOCUMENT
void amino_acid::set_letter_code_or_name(const string &prm_letter_code_or_name ///< The 1 letter code, three letter code or name to which the amino_acid should be set
                                         ) {
	// If the argument matches any names, then set to the corresponding index and return
	if ( prm_letter_code_or_name.length() > 3 ) {
		const auto index_itr = INDEX_OF_NAME().find( prm_letter_code_or_name );
		if ( index_itr != common::cend( INDEX_OF_NAME() ) ) {
			data = index_itr->second;
			return;
		}
	}

	// If the argument has one character and matches any of the LETTERS, then set to the corresponding index and return
	if ( prm_letter_code_or_name.length() == 1 ) {
		data = get_letter_index( prm_letter_code_or_name.front() );
		return;
	}

	// If the argument has three characters and matches any of the codes,  then set to the corresponding index and return
	if ( prm_letter_code_or_name.length() == 3 ) {
		const auto index_itr = INDEX_OF_CODE().find( prm_letter_code_or_name );
		if ( index_itr != common::cend( INDEX_OF_CODE() ) ) {
			data = index_itr->second;
			return;
		}
	}

	// Check whether the string represents a DNA/RNA base
	if ( prm_letter_code_or_name.length() == 3 ) {
		if ( contains( DNA_RNA_RESIDUE_NAMES, prm_letter_code_or_name ) ) {
			data = DNA_RNA_RESIDUE_NAMES.at( prm_letter_code_or_name );
		}
		// Hacks to handle erroneous ATOM AA in versions of 2yjd from before it was fixed in January 2017
		//
		// From http://deposit.rcsb.org/format-faq-v1.html :
		//
		// > Noteworthy exceptions to the above treatment of modified residues are the cases of
		// > acetylation of the N-terminus (residue ACE) and amidation of the C-terminus (residue NH2).
		else if ( prm_letter_code_or_name == "ACE" ) { data = hetatm_variant_t{ { 'A', 'C', 'E' } }; }
		else if ( prm_letter_code_or_name == "NH2" ) { data = hetatm_variant_t{ { 'N', 'H', '2' } }; }
		else {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Amino acid string \"" + prm_letter_code_or_name + "\" has three characters but is not a recognised code (currently case-sensitive)"));
		}
		return;
	}

	// No code or name has been recognised so throw a wobbly
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Amino acid name not recognised"));
}

/// \brief Ctor for amino_acid
///
/// This is deliberately not explicit to allow implicit conversions from strings to amino_acid object
amino_acid::amino_acid(const string &prm_code_or_name
                       ) {
	set_letter_code_or_name(prm_code_or_name);
}

/// \brief Make a vector of amino acids from a vector of amino acid chars
///
/// \relates amino_acid
amino_acid_vec cath::make_amino_acids_of_chars(const char_vec &prm_amino_acid_chars ///< A vector of chars for amino acids
                                               ) {
	return transform_build<amino_acid_vec>(
		prm_amino_acid_chars,
		[] (const char &x) { return amino_acid{ x }; }
	);
}

/// \brief Get the three-letter-code char_3_arr associated with the specified one letter
///
/// eg 'A' -> "ALA"
///
/// \relates amino_acid
char_3_arr cath::get_code_of_amino_acid_letter(const char &prm_one_letter_aa ///< The single-letter amino acid (eg 'A' for alanine)
                                               ) {
	return amino_acid( prm_one_letter_aa ).get_code();
}

/// \brief Get the three-letter-code string associated with the specified one letter
///
/// eg 'A' -> "ALA"
///
/// \relates amino_acid
string cath::get_code_str_of_amino_acid_letter(const char &prm_one_letter_aa ///< The single-letter amino acid (eg 'A' for alanine)
                                               ) {
	return char_arr_to_string( get_code_of_amino_acid_letter( prm_one_letter_aa ) );
}

/// \brief TODOCUMENT
///
/// \relates amino_acid
char cath::get_letter_of_amino_acid_code(const string &prm_three_letter_aa ///< TODOCUMENT
                                         ) {
	return *( amino_acid( prm_three_letter_aa ).get_letter_if_amino_acid() );
}

/// \brief Insert a description of the specified amino_acid into the specified ostream
///
/// \relates amino_acid
ostream & cath::operator<<(ostream          &prm_os,        ///< The ostream into which the description should be inserted
                           const amino_acid &prm_amino_acid ///< The amino_acid to describe
                           ) {
	prm_os << get_code_string( prm_amino_acid );
	return prm_os;
}

/// \brief TODOCUMENT
///
/// \relates amino_acid
istream & cath::operator>>(istream    &is,            ///< The istream from which to parse the amino_acid
                           amino_acid &prm_amino_acid ///< The amino acid to populate
                           ) {
	string input_string;
	is >> input_string;

	try {
		prm_amino_acid = amino_acid(input_string);
	}
	catch (const invalid_argument_exception &) {
		BOOST_THROW_EXCEPTION(invalid_argument("invalid_evaluator_argument"));
		exit(1);
	}
	return is;
}
