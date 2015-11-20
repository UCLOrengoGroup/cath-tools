/// \file
/// \brief The hmmer_scores_entry class definitions

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "common/boost_addenda/string_algorithm/split_build.h"
#include "common/type_aliases.h"
#include "exception/runtime_error_exception.h"
#include "file/hmmer_scores_file/hmmer_scores_entry.h"

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
string cath::file::detail::strip_header_name(const string &arg_string ///< TODOCUMENT
                                             ) {
	const regex the_regex( R"(^cath\|\d_\d_\d\|(\d\w{3}\S\d{2})/)" );
	smatch the_smatch;
	if ( regex_search( arg_string, the_smatch, the_regex ) ) {
		return arg_string.substr(
			numeric_cast<size_t>( the_smatch.position( 1 ) ),
			numeric_cast<size_t>( the_smatch.length  ( 1 ) )
		);
	}
	return arg_string;
}

/// \brief Ctor from all of the required pieces of information
hmmer_scores_entry::hmmer_scores_entry(const string &arg_name_1,               ///< TODOCUMENT
                                       const string &arg_accession_1,          ///< TODOCUMENT
                                       const string &arg_name_2,               ///< TODOCUMENT
                                       const string &arg_accession_2,          ///< TODOCUMENT
                                       const double &arg_full_sequence_evalue, ///< TODOCUMENT
                                       const double &arg_full_sequence_score,  ///< TODOCUMENT
                                       const double &arg_full_sequence_bias,   ///< TODOCUMENT
                                       const double &arg_best_1_domain_evalue, ///< TODOCUMENT
                                       const double &arg_best_1_domain_score,  ///< TODOCUMENT
                                       const double &arg_best_1_domain_bias,   ///< TODOCUMENT
                                       const double &arg_expected_num_doms,    ///< TODOCUMENT
                                       const size_t &arg_reg,                  ///< TODOCUMENT
                                       const size_t &arg_clu,                  ///< TODOCUMENT
                                       const size_t &arg_ov,                   ///< TODOCUMENT
                                       const size_t &arg_env,                  ///< TODOCUMENT
                                       const size_t &arg_dom,                  ///< TODOCUMENT
                                       const size_t &arg_rep,                  ///< TODOCUMENT
                                       const size_t &arg_inc,                  ///< TODOCUMENT
                                       const string &arg_description           ///< TODOCUMENT
                                       ) : name_1               ( arg_name_1                ),
                                           accession_1          ( arg_accession_1           ),
                                           name_2               ( arg_name_2                ),
                                           accession_2          ( arg_accession_2           ),
                                           full_sequence_evalue ( arg_full_sequence_evalue  ),
                                           full_sequence_score  ( arg_full_sequence_score   ),
                                           full_sequence_bias   ( arg_full_sequence_bias    ),
                                           best_1_domain_evalue ( arg_best_1_domain_evalue  ),
                                           best_1_domain_score  ( arg_best_1_domain_score   ),
                                           best_1_domain_bias   ( arg_best_1_domain_bias    ),
                                           expected_num_doms    ( arg_expected_num_doms     ),
                                           reg                  ( arg_reg                   ),
                                           clu                  ( arg_clu                   ),
                                           ov                   ( arg_ov                    ),
                                           env                  ( arg_env                   ),
                                           dom                  ( arg_dom                   ),
                                           rep                  ( arg_rep                   ),
                                           inc                  ( arg_inc                   ),
                                           description          ( arg_description           ) {
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
bool cath::file::operator==(const hmmer_scores_entry &arg_entry_a, ///< The first  hmmer_scores_entry to compare
                            const hmmer_scores_entry &arg_entry_b  ///< The second hmmer_scores_entry to compare
                            ) {
	return (
		arg_entry_a.get_name_1()               == arg_entry_b.get_name_1()
		&&
		arg_entry_a.get_accession_1()          == arg_entry_b.get_accession_1()
		&&
		arg_entry_a.get_name_2()               == arg_entry_b.get_name_2()
		&&
		arg_entry_a.get_accession_2()          == arg_entry_b.get_accession_2()
		&&
		arg_entry_a.get_full_sequence_evalue() == arg_entry_b.get_full_sequence_evalue()
		&&
		arg_entry_a.get_full_sequence_score()  == arg_entry_b.get_full_sequence_score()
		&&
		arg_entry_a.get_full_sequence_bias()   == arg_entry_b.get_full_sequence_bias()
		&&
		arg_entry_a.get_best_1_domain_evalue() == arg_entry_b.get_best_1_domain_evalue()
		&&
		arg_entry_a.get_best_1_domain_score()  == arg_entry_b.get_best_1_domain_score()
		&&
		arg_entry_a.get_best_1_domain_bias()   == arg_entry_b.get_best_1_domain_bias()
		&&
		arg_entry_a.get_expected_num_doms()    == arg_entry_b.get_expected_num_doms()
		&&
		arg_entry_a.get_reg()                  == arg_entry_b.get_reg()
		&&
		arg_entry_a.get_clu()                  == arg_entry_b.get_clu()
		&&
		arg_entry_a.get_ov()                   == arg_entry_b.get_ov()
		&&
		arg_entry_a.get_env()                  == arg_entry_b.get_env()
		&&
		arg_entry_a.get_dom()                  == arg_entry_b.get_dom()
		&&
		arg_entry_a.get_rep()                  == arg_entry_b.get_rep()
		&&
		arg_entry_a.get_inc()                  == arg_entry_b.get_inc()
		&&
		arg_entry_a.get_description()          == arg_entry_b.get_description()
	);
}

/// \brief Parse a hmmer_scores_entry from a PRC scores line as output by PRC
///
/// An example line:
///
///     cath|4_0_0|102mA00/0-153-i5                                   -          cath|4_0_0|3ixfA00/1-137 -            5.4e-09   34.1   0.1   1.1e-08   33.1   0.1   1.4   1   1   0   1   1   1   1 -
///
/// \relates hmmer_scores_entry
hmmer_scores_entry cath::file::hmmer_scores_entry_from_line(const string              &arg_hmmer_line,         ///< The line from which to parse the data
                                                            const hmmer_name_handling &arg_hmmer_name_handling ///< TODOCUMENT
                                                            ) {
	const auto line_parts = split_build<str_vec>( arg_hmmer_line, is_space(), token_compress_on );
	if ( line_parts.size() != 19 ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to parse hmmer_scores_entry from line that doesn't contain 19 parts"));
	}

	const auto &raw_name1     = line_parts[ 0 ];
	const auto &raw_name2     = line_parts[ 2 ];
	const bool  strip_headers = ( arg_hmmer_name_handling == hmmer_name_handling::STRIP );
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
string cath::file::to_string(const hmmer_scores_entry &arg_hmmer_scores_entry ///< The hmmer_scores_entry to be output as a string
                             ) {
	return "hmmer_scores_entry[" +                 arg_hmmer_scores_entry.get_name_1()
	     + ", "                  +                 arg_hmmer_scores_entry.get_accession_1()
	     + ", "                  +                 arg_hmmer_scores_entry.get_name_2()
	     + ", "                  +                 arg_hmmer_scores_entry.get_accession_2()
	     + ", "                  + std::to_string( arg_hmmer_scores_entry.get_full_sequence_evalue()  )
	     + ", "                  + std::to_string( arg_hmmer_scores_entry.get_full_sequence_score()   )
	     + ", "                  + std::to_string( arg_hmmer_scores_entry.get_full_sequence_bias()    )
	     + ", "                  + std::to_string( arg_hmmer_scores_entry.get_best_1_domain_evalue()  )
	     + ", "                  + std::to_string( arg_hmmer_scores_entry.get_best_1_domain_score()   )
	     + ", "                  + std::to_string( arg_hmmer_scores_entry.get_full_sequence_bias()    )
	     + ", "                  + std::to_string( arg_hmmer_scores_entry.get_expected_num_doms()     )
	     + ", "                  + std::to_string( arg_hmmer_scores_entry.get_reg()                   )
	     + ", "                  + std::to_string( arg_hmmer_scores_entry.get_clu()                   )
	     + ", "                  + std::to_string( arg_hmmer_scores_entry.get_ov()                    )
	     + ", "                  + std::to_string( arg_hmmer_scores_entry.get_env()                   )
	     + ", "                  + std::to_string( arg_hmmer_scores_entry.get_dom()                   )
	     + ", "                  + std::to_string( arg_hmmer_scores_entry.get_rep()                   )
	     + ", "                  + std::to_string( arg_hmmer_scores_entry.get_inc()                   )
	     + ", "                  +                 arg_hmmer_scores_entry.get_description()
	     + "]";
}

/// \brief Simple insertion operator for hmmer_scores_entry
///
/// \relates hmmer_scores_entry
ostream & cath::file::operator<<(ostream                &arg_os,   ///< The ostream to which the hmmer_scores_entry should be output
                                 const hmmer_scores_entry &arg_entry ///< The hmmer_scores_entry to output
                                 ) {
	arg_os << to_string( arg_entry );
	return arg_os;
}
