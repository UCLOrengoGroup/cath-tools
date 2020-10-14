/// \file
/// \brief The hmmer_scores_entry class definitions

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

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "cath/common/boost_addenda/string_algorithm/split_build.hpp"
#include "cath/common/exception/runtime_error_exception.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/file/hmmer_scores_file/hmmer_scores_entry.hpp"

#include <iostream> /// ***** TEMPORARY *****
#include <regex>

using namespace cath::common;
using namespace cath::file;
using namespace std;

using boost::algorithm::is_space;
using boost::algorithm::starts_with;
using boost::numeric_cast;
using boost::token_compress_on;

/// \brief TODOCUMENT
string cath::file::detail::strip_header_name(const string &prm_string ///< TODOCUMENT
                                             ) {
	const regex the_regex( R"(^cath\|\d_\d_\d\|(\d\w{3}\S\d{2})/)" );
	smatch the_smatch;
	if ( regex_search( prm_string, the_smatch, the_regex ) ) {
		return prm_string.substr(
			numeric_cast<size_t>( the_smatch.position( 1 ) ),
			numeric_cast<size_t>( the_smatch.length  ( 1 ) )
		);
	}
	return prm_string;
}

/// \brief Ctor from all of the required pieces of information
hmmer_scores_entry::hmmer_scores_entry(string        prm_name_1,               ///< TODOCUMENT
                                       string        prm_accession_1,          ///< TODOCUMENT
                                       string        prm_name_2,               ///< TODOCUMENT
                                       string        prm_accession_2,          ///< TODOCUMENT
                                       const double &prm_full_sequence_evalue, ///< TODOCUMENT
                                       const double &prm_full_sequence_score,  ///< TODOCUMENT
                                       const double &prm_full_sequence_bias,   ///< TODOCUMENT
                                       const double &prm_best_1_domain_evalue, ///< TODOCUMENT
                                       const double &prm_best_1_domain_score,  ///< TODOCUMENT
                                       const double &prm_best_1_domain_bias,   ///< TODOCUMENT
                                       const double &prm_expected_num_doms,    ///< TODOCUMENT
                                       const size_t &prm_reg,                  ///< TODOCUMENT
                                       const size_t &prm_clu,                  ///< TODOCUMENT
                                       const size_t &prm_ov,                   ///< TODOCUMENT
                                       const size_t &prm_env,                  ///< TODOCUMENT
                                       const size_t &prm_dom,                  ///< TODOCUMENT
                                       const size_t &prm_rep,                  ///< TODOCUMENT
                                       const size_t &prm_inc,                  ///< TODOCUMENT
                                       string        prm_description           ///< TODOCUMENT
                                       ) : name_1               { std::move( prm_name_1      ) },
                                           accession_1          { std::move( prm_accession_1 ) },
                                           name_2               { std::move( prm_name_2      ) },
                                           accession_2          { std::move( prm_accession_2 ) },
                                           full_sequence_evalue { prm_full_sequence_evalue     },
                                           full_sequence_score  { prm_full_sequence_score      },
                                           full_sequence_bias   { prm_full_sequence_bias       },
                                           best_1_domain_evalue { prm_best_1_domain_evalue     },
                                           best_1_domain_score  { prm_best_1_domain_score      },
                                           best_1_domain_bias   { prm_best_1_domain_bias       },
                                           expected_num_doms    { prm_expected_num_doms        },
                                           reg                  { prm_reg                      },
                                           clu                  { prm_clu                      },
                                           ov                   { prm_ov                       },
                                           env                  { prm_env                      },
                                           dom                  { prm_dom                      },
                                           rep                  { prm_rep                      },
                                           inc                  { prm_inc                      },
                                           description          { std::move( prm_description ) } {
}

/// \brief TODOCUMENT
string hmmer_scores_entry::get_name_1() const {
	return name_1;
}

/// \brief TODOCUMENT
string hmmer_scores_entry::get_accession_1() const {
	return accession_1;
}

/// \brief TODOCUMENT
string hmmer_scores_entry::get_name_2() const {
	return name_2;
}

/// \brief TODOCUMENT
string hmmer_scores_entry::get_accession_2() const {
	return accession_2;
}

/// \brief TODOCUMENT
double hmmer_scores_entry::get_full_sequence_evalue() const {
	return full_sequence_evalue;
}

/// \brief TODOCUMENT
double hmmer_scores_entry::get_full_sequence_score() const {
	return full_sequence_score;
}

/// \brief TODOCUMENT
double hmmer_scores_entry::get_full_sequence_bias() const {
	return full_sequence_bias;
}

/// \brief TODOCUMENT
double hmmer_scores_entry::get_best_1_domain_evalue() const {
	return best_1_domain_evalue;
}

/// \brief TODOCUMENT
double hmmer_scores_entry::get_best_1_domain_score() const {
	return best_1_domain_score;
}

/// \brief TODOCUMENT
double hmmer_scores_entry::get_best_1_domain_bias() const {
	return best_1_domain_bias;
}

/// \brief TODOCUMENT
double hmmer_scores_entry::get_expected_num_doms() const {
	return expected_num_doms;
}

/// \brief TODOCUMENT
size_t hmmer_scores_entry::get_reg() const {
	return reg;
}

/// \brief TODOCUMENT
size_t hmmer_scores_entry::get_clu() const {
	return clu;
}

/// \brief TODOCUMENT
size_t hmmer_scores_entry::get_ov() const {
	return ov;
}

/// \brief TODOCUMENT
size_t hmmer_scores_entry::get_env() const {
	return env;
}

/// \brief TODOCUMENT
size_t hmmer_scores_entry::get_dom() const {
	return dom;
}

/// \brief TODOCUMENT
size_t hmmer_scores_entry::get_rep() const {
	return rep;
}

/// \brief TODOCUMENT
size_t hmmer_scores_entry::get_inc() const {
	return inc;
}

/// \brief TODOCUMENT
string hmmer_scores_entry::get_description() const {
	return description;
}

/// \brief Non-member equality operator for hmmer_scores_entry
///
/// \relates hmmer_scores_entry
bool cath::file::operator==(const hmmer_scores_entry &prm_entry_a, ///< The first  hmmer_scores_entry to compare
                            const hmmer_scores_entry &prm_entry_b  ///< The second hmmer_scores_entry to compare
                            ) {
	return (
		prm_entry_a.get_name_1()               == prm_entry_b.get_name_1()
		&&
		prm_entry_a.get_accession_1()          == prm_entry_b.get_accession_1()
		&&
		prm_entry_a.get_name_2()               == prm_entry_b.get_name_2()
		&&
		prm_entry_a.get_accession_2()          == prm_entry_b.get_accession_2()
		&&
		prm_entry_a.get_full_sequence_evalue() == prm_entry_b.get_full_sequence_evalue()
		&&
		prm_entry_a.get_full_sequence_score()  == prm_entry_b.get_full_sequence_score()
		&&
		prm_entry_a.get_full_sequence_bias()   == prm_entry_b.get_full_sequence_bias()
		&&
		prm_entry_a.get_best_1_domain_evalue() == prm_entry_b.get_best_1_domain_evalue()
		&&
		prm_entry_a.get_best_1_domain_score()  == prm_entry_b.get_best_1_domain_score()
		&&
		prm_entry_a.get_best_1_domain_bias()   == prm_entry_b.get_best_1_domain_bias()
		&&
		prm_entry_a.get_expected_num_doms()    == prm_entry_b.get_expected_num_doms()
		&&
		prm_entry_a.get_reg()                  == prm_entry_b.get_reg()
		&&
		prm_entry_a.get_clu()                  == prm_entry_b.get_clu()
		&&
		prm_entry_a.get_ov()                   == prm_entry_b.get_ov()
		&&
		prm_entry_a.get_env()                  == prm_entry_b.get_env()
		&&
		prm_entry_a.get_dom()                  == prm_entry_b.get_dom()
		&&
		prm_entry_a.get_rep()                  == prm_entry_b.get_rep()
		&&
		prm_entry_a.get_inc()                  == prm_entry_b.get_inc()
		&&
		prm_entry_a.get_description()          == prm_entry_b.get_description()
	);
}

/// \brief Parse a hmmer_scores_entry from a PRC scores line as output by PRC
///
/// An example line:
///
///     cath|4_0_0|102mA00/0-153-i5                                   -          cath|4_0_0|3ixfA00/1-137 -            5.4e-09   34.1   0.1   1.1e-08   33.1   0.1   1.4   1   1   0   1   1   1   1 -
///
/// \relates hmmer_scores_entry
hmmer_scores_entry cath::file::hmmer_scores_entry_from_line(const string              &prm_hmmer_line,         ///< The line from which to parse the data
                                                            const hmmer_name_handling &prm_hmmer_name_handling ///< TODOCUMENT
                                                            ) {
	const auto line_parts = split_build<str_vec>( prm_hmmer_line, is_space(), token_compress_on );
	if ( line_parts.size() != 19 ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to parse hmmer_scores_entry from line that doesn't contain 19 parts"));
	}

	const auto &raw_name1     = line_parts[ 0 ];
	const auto &raw_name2     = line_parts[ 2 ];
	const bool  strip_headers = ( prm_hmmer_name_handling == hmmer_name_handling::STRIP );
	const auto  name1         = strip_headers ? detail::strip_header_name( raw_name1 ) : raw_name1;
	const auto  name2         = strip_headers ? detail::strip_header_name( raw_name2 ) : raw_name2;

	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return hmmer_scores_entry{
		       name1,
		       line_parts[  1 ],
		       name2,
		       line_parts[  3 ],
		stod ( line_parts[  4 ] ),
		stod ( line_parts[  5 ] ),
		stod ( line_parts[  6 ] ),
		stod ( line_parts[  7 ] ),
		stod ( line_parts[  8 ] ),
		stod ( line_parts[  9 ] ),
		stod ( line_parts[ 10 ] ),
		stoul( line_parts[ 11 ] ),
		stoul( line_parts[ 12 ] ),
		stoul( line_parts[ 13 ] ),
		stoul( line_parts[ 14 ] ),
		stoul( line_parts[ 15 ] ),
		stoul( line_parts[ 16 ] ),
		stoul( line_parts[ 17 ] ),
		       line_parts[ 18 ]
	};
}

/// \brief Simple to_string() overload for hmmer_scores_entry
///
/// \relates hmmer_scores_entry
string cath::file::to_string(const hmmer_scores_entry &prm_hmmer_scores_entry ///< The hmmer_scores_entry to be output as a string
                             ) {
	return "hmmer_scores_entry[" +                 prm_hmmer_scores_entry.get_name_1()
	     + ", "                  +                 prm_hmmer_scores_entry.get_accession_1()
	     + ", "                  +                 prm_hmmer_scores_entry.get_name_2()
	     + ", "                  +                 prm_hmmer_scores_entry.get_accession_2()
	     + ", "                  + std::to_string( prm_hmmer_scores_entry.get_full_sequence_evalue()  )
	     + ", "                  + std::to_string( prm_hmmer_scores_entry.get_full_sequence_score()   )
	     + ", "                  + std::to_string( prm_hmmer_scores_entry.get_full_sequence_bias()    )
	     + ", "                  + std::to_string( prm_hmmer_scores_entry.get_best_1_domain_evalue()  )
	     + ", "                  + std::to_string( prm_hmmer_scores_entry.get_best_1_domain_score()   )
	     + ", "                  + std::to_string( prm_hmmer_scores_entry.get_full_sequence_bias()    )
	     + ", "                  + std::to_string( prm_hmmer_scores_entry.get_expected_num_doms()     )
	     + ", "                  + std::to_string( prm_hmmer_scores_entry.get_reg()                   )
	     + ", "                  + std::to_string( prm_hmmer_scores_entry.get_clu()                   )
	     + ", "                  + std::to_string( prm_hmmer_scores_entry.get_ov()                    )
	     + ", "                  + std::to_string( prm_hmmer_scores_entry.get_env()                   )
	     + ", "                  + std::to_string( prm_hmmer_scores_entry.get_dom()                   )
	     + ", "                  + std::to_string( prm_hmmer_scores_entry.get_rep()                   )
	     + ", "                  + std::to_string( prm_hmmer_scores_entry.get_inc()                   )
	     + ", "                  +                 prm_hmmer_scores_entry.get_description()
	     + "]";
}

/// \brief Simple insertion operator for hmmer_scores_entry
///
/// \relates hmmer_scores_entry
ostream & cath::file::operator<<(ostream                &prm_os,   ///< The ostream to which the hmmer_scores_entry should be output
                                 const hmmer_scores_entry &prm_entry ///< The hmmer_scores_entry to output
                                 ) {
	prm_os << to_string( prm_entry );
	return prm_os;
}
