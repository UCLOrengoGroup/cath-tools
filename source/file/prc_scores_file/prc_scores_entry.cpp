/// \file
/// \brief The prc_scores_entry class definitions

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/finder.hpp>

#include "common/boost_addenda/string_algorithm/split_build.h"
#include "common/type_aliases.h"
#include "exception/runtime_error_exception.h"
#include "file/prc_scores_file/prc_scores_entry.h"

using namespace cath::common;
using namespace cath::file;
using namespace std;

using boost::algorithm::is_space;
using boost::token_compress_on;

/// \brief Ctor from all of the required pieces of information
prc_scores_entry::prc_scores_entry(const string &arg_name_1,   ///< Name of protein 1
                                   const size_t &arg_start_1,  ///< Start on protein 1
                                   const size_t &arg_end_1,    ///< End on protein 1
                                   const size_t &arg_length_1, ///< Length of protein 1
                                   const size_t &arg_hit_num,  ///< Number of this particular hit
                                   const string &arg_name_2,   ///< Name of protein 2
                                   const size_t &arg_start_2,  ///< Start on protein 2
                                   const size_t &arg_end_2,    ///< End on protein 2
                                   const size_t &arg_length_2, ///< Length of protein 2
                                   const double &arg_simple,   ///< Simple score
                                   const double &arg_reverse,  ///< Reverse score
                                   const double &arg_evalue    ///< E-value
                                   ) : name_1   ( arg_name_1   ),
                                       start_1  ( arg_start_1  ),
                                       end_1    ( arg_end_1    ),
                                       length_1 ( arg_length_1 ),
                                       hit_num  ( arg_hit_num  ),
                                       name_2   ( arg_name_2   ),
                                       start_2  ( arg_start_2  ),
                                       end_2    ( arg_end_2    ),
                                       length_2 ( arg_length_2 ),
                                       simple   ( arg_simple   ),
                                       reverse  ( arg_reverse  ),
                                       evalue   ( arg_evalue   ) {
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
bool cath::file::operator==(const prc_scores_entry &arg_entry_a, ///< The first  prc_scores_entry to compare
                            const prc_scores_entry &arg_entry_b  ///< The second prc_scores_entry to compare
                            ) {
	return (
		arg_entry_a.get_name_1()   == arg_entry_b.get_name_1()
		&&
		arg_entry_a.get_name_2()   == arg_entry_b.get_name_2()
		&&
		arg_entry_a.get_start_1()  == arg_entry_b.get_start_1()
		&&
		arg_entry_a.get_start_2()  == arg_entry_b.get_start_2()
		&&
		arg_entry_a.get_end_1()    == arg_entry_b.get_end_1()
		&&
		arg_entry_a.get_end_2()    == arg_entry_b.get_end_2()
		&&
		arg_entry_a.get_length_1() == arg_entry_b.get_length_1()
		&&
		arg_entry_a.get_length_2() == arg_entry_b.get_length_2()
		&&
		arg_entry_a.get_hit_num()  == arg_entry_b.get_hit_num()
		&&
		arg_entry_a.get_simple()   == arg_entry_b.get_simple()
		&&
		arg_entry_a.get_reverse()  == arg_entry_b.get_reverse()
		&&
		arg_entry_a.get_evalue()   == arg_entry_b.get_evalue()
	);
}

/// \brief Parse a prc_scores_entry from a PRC scores line as output by PRC
///
/// An example line:
///
///     1i4dA00 4       199     201     1       3cazA00 15      209     219       25.3    16.3   1.6e-11
///
/// \relates prc_scores_entry
prc_scores_entry cath::file::prc_scores_entry_from_line(const string &arg_prc_line ///< The line from which to parse the data
                                                        ) {
	const auto line_parts = split_build<str_vec>( arg_prc_line, is_space(), token_compress_on );
	if ( line_parts.size() != 12 ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to parse prc_scores_entry from line that doesn't contain 12 parts"));
	}

	/// \todo Come C++17, if Herb Sutter has gotten his way (n4029), just use braced list here
	return prc_scores_entry{
		       line_parts[  0 ],
		stoul( line_parts[  1 ] ),
		stoul( line_parts[  2 ] ),
		stoul( line_parts[  3 ] ),
		stoul( line_parts[  4 ] ),
		       line_parts[  5 ],
		stoul( line_parts[  6 ] ),
		stoul( line_parts[  7 ] ),
		stoul( line_parts[  8 ] ),
		stod ( line_parts[  9 ] ),
		stod ( line_parts[ 10 ] ),
		stod ( line_parts[ 11 ] )
	};
}

/// \brief Simple to_string() overload for prc_scores_entry
///
/// \relates prc_scores_entry
string cath::file::to_string(const prc_scores_entry &arg_prc_scores_entry ///< The prc_scores_entry to be output as a string
                             ) {
	return "prc_scores_entry[" +                 arg_prc_scores_entry.get_name_1()
	     + ", "                + std::to_string( arg_prc_scores_entry.get_start_1()  )
	     + ", "                + std::to_string( arg_prc_scores_entry.get_end_1()    )
	     + ", "                + std::to_string( arg_prc_scores_entry.get_length_1() )
	     + ", "                + std::to_string( arg_prc_scores_entry.get_hit_num()  )
	     + ", "                +                 arg_prc_scores_entry.get_name_2()
	     + ", "                + std::to_string( arg_prc_scores_entry.get_start_2()  )
	     + ", "                + std::to_string( arg_prc_scores_entry.get_end_2()    )
	     + ", "                + std::to_string( arg_prc_scores_entry.get_length_2() )
	     + ", "                + std::to_string( arg_prc_scores_entry.get_simple()   )
	     + ", "                + std::to_string( arg_prc_scores_entry.get_reverse()  )
	     + ", "                + std::to_string( arg_prc_scores_entry.get_evalue()   )
	     + "]";
}

/// \brief Simple insertion operator for prc_scores_entry
///
/// \relates prc_scores_entry
ostream & cath::file::operator<<(ostream                &arg_os,   ///< The ostream to which the prc_scores_entry should be output
                                 const prc_scores_entry &arg_entry ///< The prc_scores_entry to output
                                 ) {
	arg_os << to_string( arg_entry );
	return arg_os;
}
