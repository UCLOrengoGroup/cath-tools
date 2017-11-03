/// \file
/// \brief The ssap_scores_entry class definitions

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

#include "ssap_scores_entry.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/finder.hpp>

#include "common/boost_addenda/string_algorithm/split_build.hpp"
#include "common/exception/runtime_error_exception.hpp"
#include "common/type_aliases.hpp"

#include <string>

using namespace cath::common;
using namespace cath::file;
using namespace std;

using boost::algorithm::is_space;
using boost::token_compress_on;

/// \brief Ctor from all of the required pieces of information
ssap_scores_entry::ssap_scores_entry(string        arg_prot1,      ///< Name of protein 1
                                     string        arg_prot2,      ///< Name of protein 2
                                     const size_t &arg_length1,    ///< Length of protein 1
                                     const size_t &arg_length2,    ///< Length of protein 2
                                     const double &arg_ssap_score, ///< SSAP score for structural comparison
                                     const size_t &arg_num_equivs, ///< Number of equivalent/aligned residues
                                     const double &arg_overlap_pc, ///< Percentage overlap  (100% x overlap /length of largest)
                                     const double &arg_seq_id_pc,  ///< Percentage identity (100% x identity/length of smallest)
                                     const double &arg_rmsd        ///< RMSD of superposed structures
                                     ) : name_1     { std::move( arg_prot1 ) },
                                         name_2     { std::move( arg_prot2 ) },
                                         length_1   { arg_length1            },
                                         length_2   { arg_length2            },
                                         ssap_score { arg_ssap_score         },
                                         num_equivs { arg_num_equivs         },
                                         overlap_pc { arg_overlap_pc         },
                                         seq_id_pc  { arg_seq_id_pc          },
                                         rmsd       { arg_rmsd               } {
}

/// \brief Getter for the name of protein 1
const string & ssap_scores_entry::get_name_1() const {
	return name_1;
}

/// \brief Getter for the name of protein 2
const string & ssap_scores_entry::get_name_2() const {
	return name_2;
}

/// \brief Getter for the length of protein 1
const size_t & ssap_scores_entry::get_length_1() const {
	return length_1;
}

/// \brief Getter for the length of protein 2
const size_t & ssap_scores_entry::get_length_2() const {
	return length_2;
}

/// \brief Getter for the SSAP score for structural comparison
const double & ssap_scores_entry::get_ssap_score() const {
	return ssap_score;
}

/// \brief Getter for the number of equivalent/aligned residues
const size_t & ssap_scores_entry::get_num_equivs() const {
	return num_equivs;
}

/// \brief Getter for the percentage overlap  (100% x overlap /length of largest)
const double & ssap_scores_entry::get_overlap_pc() const {
	return overlap_pc;
}

/// \brief Getter for the Percentage identity (100% x identity/length of smallest)
const double & ssap_scores_entry::get_seq_id_pc() const {
	return seq_id_pc;
}

/// \brief Getter for the RMSD of superposed structures
const double & ssap_scores_entry::get_rmsd() const {
	return rmsd;
}

/// \brief Non-member equality operator for ssap_scores_entry
///
/// \relates ssap_scores_entry
bool cath::file::operator==(const ssap_scores_entry &arg_entry_a, ///< The first  ssap_scores_entry to compare
                            const ssap_scores_entry &arg_entry_b  ///< The second ssap_scores_entry to compare
                            ) {
	return (
		arg_entry_a.get_name_1()     == arg_entry_b.get_name_1()
		&&
		arg_entry_a.get_name_2()     == arg_entry_b.get_name_2()
		&&
		arg_entry_a.get_length_1()   == arg_entry_b.get_length_1()
		&&
		arg_entry_a.get_length_2()   == arg_entry_b.get_length_2()
		&&
		arg_entry_a.get_ssap_score() == arg_entry_b.get_ssap_score()
		&&
		arg_entry_a.get_num_equivs() == arg_entry_b.get_num_equivs()
		&&
		arg_entry_a.get_overlap_pc() == arg_entry_b.get_overlap_pc()
		&&
		arg_entry_a.get_seq_id_pc()  == arg_entry_b.get_seq_id_pc()
		&&
		arg_entry_a.get_rmsd()       == arg_entry_b.get_rmsd()
	);
}

/// \brief Parse a ssap_scores_entry from a SSAP scores line as output by SSAP
///
/// An example line:
///
///     1cukA03  1hjpA03   48   44  94.92   44   91   97   0.71
///
/// \relates ssap_scores_entry
ssap_scores_entry cath::file::ssap_scores_entry_from_line(const string &arg_ssap_line_entry ///< The line from which to parse the data
                                                          ) {
	const auto line_parts = split_build<str_vec>( arg_ssap_line_entry, is_space(), token_compress_on );
	if ( line_parts.size() != 9 ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to parse ssap_scores_entry from line that doesn't contain 9 parts"));
	}

	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return ssap_scores_entry{
		       line_parts[ 0 ],
		       line_parts[ 1 ],
		stoul( line_parts[ 2 ] ),
		stoul( line_parts[ 3 ] ),
		stod ( line_parts[ 4 ] ),
		stoul( line_parts[ 5 ] ),
		stod ( line_parts[ 6 ] ),
		stod ( line_parts[ 7 ] ),
		stod ( line_parts[ 8 ] )
	};
}

/// \brief Simple to_string() overload for ssap_scores_entry
///
/// \relates ssap_scores_entry
string cath::file::to_string(const ssap_scores_entry &arg_ssap_scores_entry ///< The ssap_scores_entry to be output as a string
                             ) {
	return "ssap_scores_entry[" +                   arg_ssap_scores_entry.get_name_1()
	     + ", "                 +                   arg_ssap_scores_entry.get_name_2()
	     + ", "                 + ::std::to_string( arg_ssap_scores_entry.get_length_1()    )
	     + ", "                 + ::std::to_string( arg_ssap_scores_entry.get_length_2()    )
	     + ", "                 + ::std::to_string( arg_ssap_scores_entry.get_ssap_score() )
	     + ", "                 + ::std::to_string( arg_ssap_scores_entry.get_num_equivs() )
	     + ", "                 + ::std::to_string( arg_ssap_scores_entry.get_overlap_pc() )
	     + ", "                 + ::std::to_string( arg_ssap_scores_entry.get_seq_id_pc()  )
	     + ", "                 + ::std::to_string( arg_ssap_scores_entry.get_rmsd()       )
	     + "]";
}

/// \brief Simple insertion operator for ssap_scores_entry
///
/// \relates ssap_scores_entry
ostream & cath::file::operator<<(ostream                 &arg_os,               ///< The ostream to which the ssap_scores_entry should be output
                                 const ssap_scores_entry &arg_ssap_scores_entry ///< The ssap_scores_entry to output
                                 ) {
	arg_os << to_string( arg_ssap_scores_entry );
	return arg_os;
}
