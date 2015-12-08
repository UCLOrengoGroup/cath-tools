/// \file
/// \brief The ssap_and_prc class definitions

#include "ssap_and_prc.h"

#include "exception/invalid_argument_exception.h"
#include "file/prc_scores_file/prc_scores_entry.h"
#include "file/ssap_scores_file/ssap_scores_entry.h"

#include <cmath>
#include <string>

using namespace cath::common;
using namespace cath::file;
using namespace cath::homcheck;
using namespace std;

/// \brief Ctor from a ssap_scores_entry and a prc_scores_entry
///
/// \pre arg_ssap and arg_prc must have matching name_1 and name_2 values,
///      else an invalid_argument_exception will be thrown
ssap_and_prc::ssap_and_prc(const ssap_scores_entry &arg_ssap, ///< The SSAP scores from which this should be constructed
                           const prc_scores_entry  &arg_prc   ///< The PRC scores from which this should be constructed
                           ) : the_ssap            ( arg_ssap ),
                               the_prc             ( arg_prc  ),
                               magic_function_score( magic_function( the_ssap, the_prc ) ) {
	if (
		( the_ssap.get_name_1() != the_prc.get_name_1() )
		||
		( the_ssap.get_name_2() != the_prc.get_name_2() )
		) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot construct a ssap_and_prc from a ssap_scores_entry and prc_scores_entry with mismatching names"));
	}
}

/// \brief Getter for the query_id (name_1) shared by the SSAP and PRC results
const string & ssap_and_prc::get_query_id() const {
	return the_ssap.get_name_1();
}

/// \brief Getter for the query_id (name_1) shared by the SSAP and PRC results
const string & ssap_and_prc::get_match_id() const {
	return the_ssap.get_name_2();
}

/// \brief Get (a const reference to) the SSAP result
const ssap_scores_entry & ssap_and_prc::get_ssap() const {
	return the_ssap;
}

/// \brief Get (a const reference to) the PRC result
const prc_scores_entry & ssap_and_prc::get_prc() const {
	return the_prc;
}

/// \brief Getter for the the magic function score (calculated during construction)
///
/// \seealso magic_function()
const double & ssap_and_prc::get_magic_function_score() const {
	return magic_function_score;
}

/// \brief Calculate the magic function associated with the specified SSAP and PRC results
///        (working on the assumption that their results reflect the same comparison)
///
/// At present, this does not check that the two results have matching IDs
///
/// The magic function is defined as: `ssap_score - log10( prc_evalue )`
double cath::homcheck::magic_function(const ssap_scores_entry &arg_ssap, ///< The SSAP result with which the magic function should be computed
                                      const prc_scores_entry  &arg_prc   ///< The PRC  result with which the magic function should be computed
					                  ) {
	return arg_ssap.get_ssap_score() - std::log10( arg_prc.get_evalue() );
}

/// \brief Getter for the SSAP length_1 of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_ssap_length_1(const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc result to query
                                                 ) {
	return arg_ssap_and_prc.get_ssap().get_length_1();
}

/// \brief Getter for the SSAP length_2 of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_ssap_length_2(const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc result to query
                                                 ) {
	return arg_ssap_and_prc.get_ssap().get_length_2();
}

/// \brief Getter for the SSAP score of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const double & cath::homcheck::get_ssap_score(const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc result to query
                                              ) {
	return arg_ssap_and_prc.get_ssap().get_ssap_score();
}

/// \brief Getter for the SSAP num_equivs of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_ssap_num_equivs(const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc result to query
                                                   ) {
	return arg_ssap_and_prc.get_ssap().get_num_equivs();
}

/// \brief Getter for the SSAP overlap_pc of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const double & cath::homcheck::get_ssap_overlap_pc(const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc result to query
                                                   ) {
	return arg_ssap_and_prc.get_ssap().get_overlap_pc();
}

/// \brief Getter for the SSAP seq_id_pc of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const double & cath::homcheck::get_ssap_seq_id_pc(const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc result to query
                                                  ) {
	return arg_ssap_and_prc.get_ssap().get_seq_id_pc();
}

/// \brief Getter for the SSAP RMSD of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const double & cath::homcheck::get_ssap_rmsd(const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc result to query
                                             ) {
	return arg_ssap_and_prc.get_ssap().get_rmsd();
}

/// \brief Getter for the PRC start_1 of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_prc_start_1(const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc result to query
                                               ) {
	return arg_ssap_and_prc.get_prc().get_start_1();
}

/// \brief Getter for the PRC end_1 of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_prc_end_1(const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc result to query
                                             ) {
	return arg_ssap_and_prc.get_prc().get_end_1();
}

/// \brief Getter for the PRC length_1 of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_prc_length_1(const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc result to query
                                                ) {
	return arg_ssap_and_prc.get_prc().get_length_1();
}

/// \brief Getter for the PRC hit_num of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_prc_hit_num(const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc result to query
                                               ) {
	return arg_ssap_and_prc.get_prc().get_hit_num();
}

/// \brief Getter for the PRC start_2 of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_prc_start_2(const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc result to query
                                               ) {
	return arg_ssap_and_prc.get_prc().get_start_2();
}

/// \brief Getter for the PRC end_2 of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_prc_end_2(const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc result to query
                                             ) {
	return arg_ssap_and_prc.get_prc().get_end_2();
}

/// \brief Getter for the PRC length_2 of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_prc_length_2(const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc result to query
                                                ) {
	return arg_ssap_and_prc.get_prc().get_length_2();
}

/// \brief Getter for the PRC simple of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const double & cath::homcheck::get_prc_simple(const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc result to query
                                              ) {
	return arg_ssap_and_prc.get_prc().get_simple();
}

/// \brief Getter for the PRC reverse of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const double & cath::homcheck::get_prc_reverse(const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc result to query
                                               ) {
	return arg_ssap_and_prc.get_prc().get_reverse();
}

/// \brief Getter for the PRC evalue of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const double & cath::homcheck::get_prc_evalue(const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc result to query
                                              ) {
	return arg_ssap_and_prc.get_prc().get_evalue();
}

/// \brief Simple to_string() overload for ssap_and_prc
///
/// \relates ssap_and_prc
string cath::homcheck::to_string(const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc to be output as a string
                                 ) {
	return "ssap_and_prc[ssap_score:"
		+ std::to_string( get_ssap_score( arg_ssap_and_prc ) )
		+ "; prc_evalue:"
		+ std::to_string( get_ssap_score( arg_ssap_and_prc ) )
		+ "; magic_function:"
		+ std::to_string( arg_ssap_and_prc.get_magic_function_score() )
		+ "]";
}

/// \brief Simple insertion operator for ssap_and_prc
///
/// \relates ssap_and_prc
ostream & cath::homcheck::operator<<(ostream            &arg_ostream,     ///< The ostream to which the ssap_and_prc should be output
                                     const ssap_and_prc &arg_ssap_and_prc ///< The ssap_and_prc to output
                                     ) {
	arg_ostream << to_string( arg_ssap_and_prc );
	return arg_ostream;
}

