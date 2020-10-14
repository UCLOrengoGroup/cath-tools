/// \file
/// \brief The prc_scores_entry class definitions

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

#include "common/type_aliases.hpp"
#include "file/prc_scores_file/detail/prc_scores_line_parser.hpp"
#include "file/prc_scores_file/prc_scores_entry.hpp"

using namespace cath::file;
using namespace std;

/// \brief Ctor from all of the required pieces of information
prc_scores_entry::prc_scores_entry(string        prm_name_1,   ///< Name of protein 1
                                   const size_t &prm_start_1,  ///< Start on protein 1
                                   const size_t &prm_end_1,    ///< End on protein 1
                                   const size_t &prm_length_1, ///< Length of protein 1
                                   const size_t &prm_hit_num,  ///< Number of this particular hit
                                   string        prm_name_2,   ///< Name of protein 2
                                   const size_t &prm_start_2,  ///< Start on protein 2
                                   const size_t &prm_end_2,    ///< End on protein 2
                                   const size_t &prm_length_2, ///< Length of protein 2
                                   const double &prm_simple,   ///< Simple score
                                   const double &prm_reverse,  ///< Reverse score
                                   const double &prm_evalue    ///< E-value
                                   ) : name_1   { std::move( prm_name_1 ) },
                                       start_1  { prm_start_1             },
                                       end_1    { prm_end_1               },
                                       length_1 { prm_length_1            },
                                       hit_num  { prm_hit_num             },
                                       name_2   { std::move( prm_name_2 ) },
                                       start_2  { prm_start_2             },
                                       end_2    { prm_end_2               },
                                       length_2 { prm_length_2            },
                                       simple   { prm_simple              },
                                       reverse  { prm_reverse             },
                                       evalue   { prm_evalue              } {
}

/// \brief Getter for the name of protein 1
const string & prc_scores_entry::get_name_1() const {
	return name_1;
}

/// \brief Start on protein 1
const size_t & prc_scores_entry::get_start_1() const {
	return start_1;
}

/// \brief End on protein 1
const size_t & prc_scores_entry::get_end_1() const {
	return end_1;
}

/// \brief Length of protein 1
const size_t & prc_scores_entry::get_length_1() const {
	return length_1;
}

/// \brief Number of this particular hit
const size_t & prc_scores_entry::get_hit_num() const {
	return hit_num;
}

/// \brief Name of protein 2
const string & prc_scores_entry::get_name_2() const {
	return name_2;
}

/// \brief Start on protein 2
const size_t & prc_scores_entry::get_start_2() const {
	return start_2;
}

/// \brief End on protein 2
const size_t & prc_scores_entry::get_end_2() const {
	return end_2;
}

/// \brief Length of protein 2
const size_t & prc_scores_entry::get_length_2() const {
	return length_2;
}

/// \brief Simple score
const double & prc_scores_entry::get_simple() const {
	return simple;
}

/// \brief Reverse score
const double & prc_scores_entry::get_reverse() const {
	return reverse;
}

/// \brief E-value
const double & prc_scores_entry::get_evalue() const {
	return evalue;
}

/// \brief Non-member equality operator for prc_scores_entry
///
/// \relates prc_scores_entry
bool cath::file::operator==(const prc_scores_entry &prm_entry_a, ///< The first  prc_scores_entry to compare
                            const prc_scores_entry &prm_entry_b  ///< The second prc_scores_entry to compare
                            ) {
	return (
		prm_entry_a.get_name_1()   == prm_entry_b.get_name_1()
		&&
		prm_entry_a.get_name_2()   == prm_entry_b.get_name_2()
		&&
		prm_entry_a.get_start_1()  == prm_entry_b.get_start_1()
		&&
		prm_entry_a.get_start_2()  == prm_entry_b.get_start_2()
		&&
		prm_entry_a.get_end_1()    == prm_entry_b.get_end_1()
		&&
		prm_entry_a.get_end_2()    == prm_entry_b.get_end_2()
		&&
		prm_entry_a.get_length_1() == prm_entry_b.get_length_1()
		&&
		prm_entry_a.get_length_2() == prm_entry_b.get_length_2()
		&&
		prm_entry_a.get_hit_num()  == prm_entry_b.get_hit_num()
		&&
		prm_entry_a.get_simple()   == prm_entry_b.get_simple()
		&&
		prm_entry_a.get_reverse()  == prm_entry_b.get_reverse()
		&&
		prm_entry_a.get_evalue()   == prm_entry_b.get_evalue()
	);
}

/// \brief Parse a prc_scores_entry from a PRC scores line as output by PRC
///
/// An example line:
///
///     1i4dA00 4       199     201     1       3cazA00 15      209     219       25.3    16.3   1.6e-11
///
/// \relates prc_scores_entry
prc_scores_entry cath::file::prc_scores_entry_from_line(const string &prm_prc_line ///< The line from which to parse the data
                                                        ) {
	return prc_scores_line_parser{}.parse_line( prm_prc_line );
}

/// \brief Simple to_string() overload for prc_scores_entry
///
/// \relates prc_scores_entry
string cath::file::to_string(const prc_scores_entry &prm_prc_scores_entry ///< The prc_scores_entry to be output as a string
                             ) {
	return "prc_scores_entry[" +                 prm_prc_scores_entry.get_name_1()
	     + ", "                + std::to_string( prm_prc_scores_entry.get_start_1()  )
	     + ", "                + std::to_string( prm_prc_scores_entry.get_end_1()    )
	     + ", "                + std::to_string( prm_prc_scores_entry.get_length_1() )
	     + ", "                + std::to_string( prm_prc_scores_entry.get_hit_num()  )
	     + ", "                +                 prm_prc_scores_entry.get_name_2()
	     + ", "                + std::to_string( prm_prc_scores_entry.get_start_2()  )
	     + ", "                + std::to_string( prm_prc_scores_entry.get_end_2()    )
	     + ", "                + std::to_string( prm_prc_scores_entry.get_length_2() )
	     + ", "                + std::to_string( prm_prc_scores_entry.get_simple()   )
	     + ", "                + std::to_string( prm_prc_scores_entry.get_reverse()  )
	     + ", "                + std::to_string( prm_prc_scores_entry.get_evalue()   )
	     + "]";
}

/// \brief Simple insertion operator for prc_scores_entry
///
/// \relates prc_scores_entry
ostream & cath::file::operator<<(ostream                &prm_os,   ///< The ostream to which the prc_scores_entry should be output
                                 const prc_scores_entry &prm_entry ///< The prc_scores_entry to output
                                 ) {
	prm_os << to_string( prm_entry );
	return prm_os;
}
