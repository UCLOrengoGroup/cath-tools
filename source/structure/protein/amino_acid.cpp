/// \file
/// \brief The amino_acid class definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "amino_acid.h"

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "common/algorithm/contains.h"
#include "common/algorithm/transform_build.h"
#include "exception/invalid_argument_exception.h"
#include "exception/out_of_range_exception.h"

#include <set>

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace std;

using boost::algorithm::all;
using boost::algorithm::is_print;
using boost::algorithm::to_upper_copy;
using boost::algorithm::trim_copy;

/// \brief TODOCUMENT
///
/// Examples:
///  * PDBs 1dgo, 1ii1, 1js7, 1qe7, 3ga6 and 4fm0 all contain
///    ATOM records with each of the residue names: DA, DC, DG, DT and DU
///  * PDBs 1eg0, 1fjf, 2i82 and 3u5f all contain
///    ATOM records with each of the residue names: A, C, G, N and U
const set<string> DNA_RNA_RESIDUE_NAMES = { "A", "C", "G", "N", "U", "DA", "DC", "DG", "DT", "DU" };

/// \brief TODOCUMENT
//const string amino_acid::UNKNOWN_AMINO_ACID_NAME("Unknown");

/// \todo What about ACE and NH2 (eg PDB 2yjd)
///       From http://deposit.rcsb.org/format-faq-v1.html :
///
/// > Noteworthy exceptions to the above treatment of modified residues are the cases of
/// > acetylation of the N-terminus (residue ACE) and amidation of the C-terminus (residue NH2).

///// \brief TODOCUMENT
//const string & global_test_constants::EXAMPLE_B_PDB_STEMNAME() {
//	static const string example_b_pdb_stemname( "1hdoA00" );
//	return example_b_pdb_stemname;
//}

/// \brief TODOCUMENT
const amino_acid::char_size_map & amino_acid::INDEX_OF_LETTER() {
	static const amino_acid::char_size_map index_of_letter(amino_acid::build_index_map<char,   0>());
	return index_of_letter;
}

/// \brief TODOCUMENT
const amino_acid::string_size_map & amino_acid::INDEX_OF_CODE() {
	static const amino_acid::string_size_map index_of_code(amino_acid::build_index_map<string, 1>());
	return index_of_code;
}

/// \brief TODOCUMENT
const amino_acid::string_size_map & amino_acid::INDEX_OF_NAME() {
	static const amino_acid::string_size_map index_of_name(amino_acid::build_index_map<string, 2>());
	return index_of_name;
}

/// \brief TODOCUMENT
void amino_acid::check_is_proper_amino_acid() const {
	if ( ! is_proper_amino_acid() ) {
		BOOST_THROW_EXCEPTION(out_of_range_exception(
			"Cannot use a generic amino_acid \""
			+ ( raw_string ? *raw_string : string("") )
			+ "\" as a proper, ATOM-record amino acid"));
	}
}


/// \brief TODOCUMENT
void amino_acid::set_letter_code_or_name(const string &arg_letter_code_or_name ///< The 1 letter code, three letter code or name to which the amino_acid should be set
                                         ) {
	// If the argument matches any names, then set to the corresponding index and return
	if ( arg_letter_code_or_name.length() > 3 && contains( INDEX_OF_NAME(), arg_letter_code_or_name ) ) {
		index = INDEX_OF_NAME().find( arg_letter_code_or_name )->second;
		return;
	}

	// If the argument has one character and matches any of the LETTERS, then set to the corresponding index and return
	if ( arg_letter_code_or_name.length() == 1 ) {
//		const string upper_letter_string = to_upper_copy(arg_letter_code_or_name);
		const char &letter = arg_letter_code_or_name.front();
		if ( contains( INDEX_OF_LETTER(), letter ) ) {
			index = INDEX_OF_LETTER().find( letter )->second;
			return;
		}
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Amino acid string \"" + arg_letter_code_or_name + "\" has one character but is not a recognised letter (currently case-sensitive)"));
	}

	// If the argument has three characters and matches any of the codes,  then set to the corresponding index and return
	if (arg_letter_code_or_name.length() == 3) {
//		const string upper_code = to_upper_copy(arg_letter_code_or_name);
		if ( contains( INDEX_OF_CODE(), arg_letter_code_or_name ) ) {
			index = INDEX_OF_CODE().find(arg_letter_code_or_name)->second;
			return;
		}
	}

	// Check whether a trimmed, upper-cased version of this string is actually a residue name storing a DNA/RNA base
	//
	// Note: this has to come after the 1-letter check so that, for example, this doesn't reject "A"
	//        before it gets handled as a valid alanine
	//
	// Note: this has to come before the throw at the end of the 3-letter-code check, so that this more informative
	//       exception gets thrown for three letter strings like " DG"
	const string upper_trimmed_copy = to_upper_copy(trim_copy(arg_letter_code_or_name));
	if ( contains( DNA_RNA_RESIDUE_NAMES, upper_trimmed_copy ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception(
			"Amino acid name is actually a DNA/RNA base (\""
			+ upper_trimmed_copy
			+ "\") rather than a valid amino acid"
		));
	}

	if (arg_letter_code_or_name.length() == 3) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Amino acid string \"" + arg_letter_code_or_name + "\" has three characters but is not a recognised code (currently case-sensitive)"));
	}


	// No code or name has been recognised so throw a wobbly
	BOOST_THROW_EXCEPTION(invalid_argument_exception("Amino acid name not recognised"));
}

/// \brief Ctor for amino_acid
amino_acid::amino_acid(const string     &arg_string,    ///< TODOCUMENT
                       const pdb_record &arg_pdb_record ///< TODOCUMENT
					   ) {
	if ( ! all( arg_string, is_print() ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot construct an amino acid from a string with non-printing characters"));
	}
	switch ( arg_pdb_record ) {
		case ( pdb_record::ATOM   ) : {
			set_letter_code_or_name( arg_string );
			break;
		}
		case ( pdb_record::HETATM ) : {
			if ( arg_string.length() != 3 ) {
				BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot create a HETATM amino acid from a string that is not 3 characters long"));
			}
			raw_string = arg_string;
			break;
		}
		default : {
			BOOST_THROW_EXCEPTION(invalid_argument_exception("Value of arg_pdb_record not recognised whilst constructing an amino_acid"));
		}
	}
	assert(   raw_string ||   index );
	assert( ! raw_string || ! index );
}

/// \brief Ctor for amino_acid
///
/// This is deliberately not explicit to allow implicit conversions from strings to amino_acid object
amino_acid::amino_acid(const string &arg_code_or_name
                       ) {
	set_letter_code_or_name(arg_code_or_name);
}

/// \brief Ctor for amino_acid
///
/// This is deliberately not explicit to allow implicit conversions from chars to amino_acid object
amino_acid::amino_acid(const char &arg_code
                       ) {
	set_letter_code_or_name(string(1, arg_code));
}

/// \brief Make a vector of amino acids from a vector of amino acid chars
///
/// \relates amino_acid
amino_acid_vec cath::make_amino_acids_of_chars(const char_vec &arg_amino_acid_chars ///< A vector of chars for amino acids
                                               ) {
	return transform_build<amino_acid_vec>(
		arg_amino_acid_chars,
		[] (const char &x) { return amino_acid{ x }; }
	);
}

/// \brief TODOCUMENT
///
/// \relates amino_acid
string cath::get_code_of_amino_acid_letter(const char &arg_one_letter_aa ///< TODOCUMENT
                                           ) {
	return amino_acid(arg_one_letter_aa).get_code();
}

/// \brief TODOCUMENT
///
/// \relates amino_acid
char cath::get_letter_of_amino_acid_code(const string &arg_three_letter_aa ///< TODOCUMENT
                                         ) {
	return amino_acid(arg_three_letter_aa).get_letter();
}

/// \brief TODOCUMENT
istream & cath::operator>>(istream    &is,            ///< The istream from which to parse the amino_acid
                           amino_acid &arg_amino_acid ///< The amino acid to populate
                           ) {
	string input_string;
	is >> input_string;

	try {
		arg_amino_acid = amino_acid(input_string);
	}
	catch (const invalid_argument_exception &) {
		BOOST_THROW_EXCEPTION(invalid_argument("invalid_evaluator_argument"));
		exit(1);
	}
	return is;
}
