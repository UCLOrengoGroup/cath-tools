/// \file
/// \brief The ssap_and_prc class definitions

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

#include "ssap_and_prc.hpp"

#include <boost/lexical_cast.hpp>

#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/file/prc_scores_file/prc_scores_entry.hpp"
#include "cath/file/ssap_scores_file/ssap_scores_entry.hpp"
#include "cath/score/score_classification/rbf_model.hpp"

#include <cmath>
#include <string>

using namespace ::cath::common;
using namespace ::cath::file;
using namespace ::cath::homcheck;
using namespace ::cath::score;

using ::boost::optional;
using ::std::ostream;
using ::std::ostringstream;
using ::std::string;

/// \brief Ctor from a ssap_scores_entry and a prc_scores_entry
///
/// \pre prm_ssap and prm_prc must have matching name_1 and name_2 values,
///      else an invalid_argument_exception will be thrown
ssap_and_prc::ssap_and_prc(ssap_scores_entry prm_ssap, ///< The SSAP scores from which this should be constructed
                           prc_scores_entry  prm_prc   ///< The PRC scores from which this should be constructed
                           ) : the_ssap             { std::move( prm_ssap               ) },
                               the_prc              { std::move( prm_prc                ) },
                               magic_function_score { magic_function( the_ssap, the_prc ) } {
	if (
		( the_ssap.get_name_1() != the_prc.get_name_1() )
		||
		( the_ssap.get_name_2() != the_prc.get_name_2() )
		) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot construct a ssap_and_prc from a ssap_scores_entry and prc_scores_entry with mismatching names"));
	}
}

/// \brief Calculate and store the SVM score using the specified SVM RBF model
void ssap_and_prc::calculate_svm_score(const rbf_model &prm_svm ///< The SVM RBF model with which to calculate the score
                                       ) {
	svm_score = get_score( prm_svm, *this );
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

/// \brief Getter for the magic function score (calculated during construction)
///
/// \seealso magic_function()
const double & ssap_and_prc::get_magic_function_score() const {
	return magic_function_score;
}

/// \brief Getter for the possibly-set SVM score
const optional<double> & ssap_and_prc::get_svm_score_opt() const {
	return svm_score;
}

/// \brief Calculate the magic function associated with the specified SSAP and PRC results
///        (working on the assumption that their results reflect the same comparison)
///
/// At present, this does not check that the two results have matching IDs
///
/// The magic function is defined as: `ssap_score - log10( prc_evalue )`
double cath::homcheck::magic_function(const ssap_scores_entry &prm_ssap, ///< The SSAP result with which the magic function should be computed
                                      const prc_scores_entry  &prm_prc   ///< The PRC  result with which the magic function should be computed
                                      ) {
	return prm_ssap.get_ssap_score() - std::log10( prm_prc.get_evalue() );
}

/// \brief Return whether this ssap_and_prc result has an SVM score set
///
/// \relates ssap_and_prc
bool cath::homcheck::has_svm_score(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc to query
                                   ) {
	return static_cast<bool>( prm_ssap_and_prc.get_svm_score_opt() );
}

/// \brief Return the SVM score set in the specified ssap_and_prc
///
/// \pre `has_svm_score( prm_ssap_and_prc )` else an invalid_argument_exception is thrown
///
/// \relates ssap_and_prc
const double & cath::homcheck::get_svm_score(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc to query
                                             ) {
	if ( ! has_svm_score( prm_ssap_and_prc ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("Cannot retrieve a ssap_and_prc result's SVM score if it hasn't been set"));
	}
	return *prm_ssap_and_prc.get_svm_score_opt();
}

/// \brief Getter for the SSAP length_1 of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_ssap_length_1(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc result to query
                                                 ) {
	return prm_ssap_and_prc.get_ssap().get_length_1();
}

/// \brief Getter for the SSAP length_2 of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_ssap_length_2(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc result to query
                                                 ) {
	return prm_ssap_and_prc.get_ssap().get_length_2();
}

/// \brief Getter for the SSAP score of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const double & cath::homcheck::get_ssap_score(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc result to query
                                              ) {
	return prm_ssap_and_prc.get_ssap().get_ssap_score();
}

/// \brief Getter for the SSAP num_equivs of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_ssap_num_equivs(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc result to query
                                                   ) {
	return prm_ssap_and_prc.get_ssap().get_num_equivs();
}

/// \brief Getter for the SSAP overlap_pc of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const double & cath::homcheck::get_ssap_overlap_pc(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc result to query
                                                   ) {
	return prm_ssap_and_prc.get_ssap().get_overlap_pc();
}

/// \brief Getter for the SSAP seq_id_pc of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const double & cath::homcheck::get_ssap_seq_id_pc(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc result to query
                                                  ) {
	return prm_ssap_and_prc.get_ssap().get_seq_id_pc();
}

/// \brief Getter for the SSAP RMSD of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const double & cath::homcheck::get_ssap_rmsd(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc result to query
                                             ) {
	return prm_ssap_and_prc.get_ssap().get_rmsd();
}

/// \brief Getter for the PRC start_1 of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_prc_start_1(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc result to query
                                               ) {
	return prm_ssap_and_prc.get_prc().get_start_1();
}

/// \brief Getter for the PRC end_1 of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_prc_end_1(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc result to query
                                             ) {
	return prm_ssap_and_prc.get_prc().get_end_1();
}

/// \brief Getter for the PRC length_1 of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_prc_length_1(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc result to query
                                                ) {
	return prm_ssap_and_prc.get_prc().get_length_1();
}

/// \brief Getter for the PRC hit_num of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_prc_hit_num(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc result to query
                                               ) {
	return prm_ssap_and_prc.get_prc().get_hit_num();
}

/// \brief Getter for the PRC start_2 of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_prc_start_2(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc result to query
                                               ) {
	return prm_ssap_and_prc.get_prc().get_start_2();
}

/// \brief Getter for the PRC end_2 of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_prc_end_2(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc result to query
                                             ) {
	return prm_ssap_and_prc.get_prc().get_end_2();
}

/// \brief Getter for the PRC length_2 of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const size_t & cath::homcheck::get_prc_length_2(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc result to query
                                                ) {
	return prm_ssap_and_prc.get_prc().get_length_2();
}

/// \brief Getter for the PRC simple of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const double & cath::homcheck::get_prc_simple(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc result to query
                                              ) {
	return prm_ssap_and_prc.get_prc().get_simple();
}

/// \brief Getter for the PRC reverse of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const double & cath::homcheck::get_prc_reverse(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc result to query
                                               ) {
	return prm_ssap_and_prc.get_prc().get_reverse();
}

/// \brief Getter for the PRC evalue of the specified ssap_and_prc object
///
/// \relates ssap_and_prc
const double & cath::homcheck::get_prc_evalue(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc result to query
                                              ) {
	return prm_ssap_and_prc.get_prc().get_evalue();
}

/// \brief Simple to_string() overload for ssap_and_prc
///
/// \relates ssap_and_prc
string cath::homcheck::to_string(const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc to be output as a string
                                 ) {
	ostringstream evalue_ss;
	evalue_ss << get_prc_evalue( prm_ssap_and_prc );
	return "ssap_and_prc[query_id:"
		+ prm_ssap_and_prc.get_query_id()
		+ "; match_id:"
		+ prm_ssap_and_prc.get_match_id()
		+ "; ssap_score:"
		+ std::to_string( get_ssap_score( prm_ssap_and_prc ) )
		+ "; prc_evalue:"
		+ evalue_ss.str()
		+ ( has_svm_score( prm_ssap_and_prc ) ? ( "; SVM:" + ::std::to_string( get_svm_score( prm_ssap_and_prc ) ) ) : "" )
		+ "; magic_function:"
		+ std::to_string( prm_ssap_and_prc.get_magic_function_score() )
		+ "]";
}

/// \brief Simple insertion operator for ssap_and_prc
///
/// \relates ssap_and_prc
ostream & cath::homcheck::operator<<(ostream            &prm_ostream,     ///< The ostream to which the ssap_and_prc should be output
                                     const ssap_and_prc &prm_ssap_and_prc ///< The ssap_and_prc to output
                                     ) {
	prm_ostream << to_string( prm_ssap_and_prc );
	return prm_ostream;
}

