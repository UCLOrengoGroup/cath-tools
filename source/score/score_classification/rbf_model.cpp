/// \file
/// \brief The rbf_model class definitions

#include "rbf_model.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/optional.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/numeric.hpp>

#include "common/algorithm/transform_build.hpp"
#include "common/boost_addenda/string_algorithm/split_build.hpp"
#include "common/file/open_fstream.hpp"
#include "common/size_t_literal.hpp"
#include "exception/not_implemented_exception.hpp"
#include "exception/runtime_error_exception.hpp"
#include "file/prc_scores_file/prc_scores_entry.hpp"
#include "file/ssap_scores_file/ssap_scores_entry.hpp"
#include "score/homcheck_tools/ssap_and_prc.hpp"
#include "score/score_classification/value_list_scaling.hpp"

#include <cmath>
#include <fstream>

using namespace cath;
using namespace cath::common;
using namespace cath::file;
using namespace cath::homcheck;
using namespace cath::score;
using namespace std;

using boost::accumulate;
using boost::algorithm::is_any_of;
using boost::algorithm::is_space;
using boost::algorithm::token_compress_on;
using boost::contains;
using boost::filesystem::path;
using boost::icontains;
using boost::irange;
using boost::numeric_cast;
using boost::optional;

constexpr value_list_scaling rbf_model::PRC_EVALUE_SCALING;
constexpr value_list_scaling rbf_model::PRC_REVERSE_SCALING;
constexpr value_list_scaling rbf_model::PRC_SIMPLE_SCALING;
constexpr value_list_scaling rbf_model::SSAP_NUM_EQUIVS_SCALING;
constexpr value_list_scaling rbf_model::SSAP_OVERLAP_PC_SCALING;
constexpr value_list_scaling rbf_model::SSAP_RMSD_SCALING;
constexpr value_list_scaling rbf_model::SSAP_SEQ_ID_PC_SCALING;
constexpr value_list_scaling rbf_model::SSAP_SCORE_SCALING;

// Temporarily suppress a false-positive GCC compiler warning
//
// It seems GCC has a problem that's currently making it
// falsely warn about parsed_gamma and parsed_b in parse_rbf_model()
// (two boost::optional<double> variables) being passed to this ctor
// when it can't confirm their values are initialised.
//
// That code does an explicit check they're both initialised and throws before
// calling this ctor if not.
//
// This comment appears to record the fix:
//
//     https://gcc.gnu.org/bugzilla/show_bug.cgi?id=47679#c32
//
// Since that comment's dated 2015-09-20, the fix should have
// made it into GCC 5.3 (2015-12-04) but not GCC 5.2 (2015-07-16).
#ifdef  __GNUC__
#ifndef __clang__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#endif

/// \brief TODOCUMENT
rbf_model::rbf_model(const double                                &arg_gamma,       ///< The SVM RBF gamma parameter
                     const double                                &arg_b,           ///< The SVM RBF b parameter
                     const vector<pair<double, double_octuple> > &arg_model_values ///< The SVM RBF support vectors
                     ) : gamma           ( arg_gamma        ),
                         b               ( arg_b            ),
                         support_vectors ( arg_model_values ) {
}

// Stop suppressing the GCC warning again
#ifdef  __GNUC__
#ifndef __clang__
#pragma GCC diagnostic pop
#endif
#endif

/// \brief Make the standard list of scores for use with an rbf_model from the specified ssap_and_prc
double_octuple rbf_model::make_standard_scores(const ssap_and_prc &arg_ssap_and_prc ///< The SSAP and PRC scores from which to make the scores to be used with an rbf_model
                                               ) {
	return make_tuple(
		scale_value_copy( PRC_EVALUE_SCALING,                     log10( get_prc_evalue     ( arg_ssap_and_prc ) ) ),
		scale_value_copy( PRC_REVERSE_SCALING,                           get_prc_reverse    ( arg_ssap_and_prc )   ),
		scale_value_copy( PRC_SIMPLE_SCALING,                            get_prc_simple     ( arg_ssap_and_prc )   ),
		scale_value_copy( SSAP_NUM_EQUIVS_SCALING, numeric_cast<double>( get_ssap_num_equivs( arg_ssap_and_prc ) ) ),
		scale_value_copy( SSAP_OVERLAP_PC_SCALING,                       get_ssap_overlap_pc( arg_ssap_and_prc )   ),
		scale_value_copy( SSAP_RMSD_SCALING,                             get_ssap_rmsd      ( arg_ssap_and_prc )   ),
		scale_value_copy( SSAP_SEQ_ID_PC_SCALING,                        get_ssap_seq_id_pc ( arg_ssap_and_prc )   ),
		scale_value_copy( SSAP_SCORE_SCALING,                            get_ssap_score     ( arg_ssap_and_prc )   )
	);
}

/// \brief Make the standard list of scores for use with an rbf_model from the specified SSAP scores and PRC scores
double_octuple rbf_model::make_standard_scores(const prc_scores_entry  &arg_prc, ///< The PRC scores from which to make the scores to be used with an rbf_model
                                               const ssap_scores_entry &arg_ssap ///< The SSAP scores from which to make the scores to be used with an rbf_model
                                               ) {
	return make_tuple(
		scale_value_copy( PRC_EVALUE_SCALING,                     log10(  arg_prc.get_evalue()     ) ),
		scale_value_copy( PRC_REVERSE_SCALING,                            arg_prc.get_reverse()      ),
		scale_value_copy( PRC_SIMPLE_SCALING,                             arg_prc.get_simple()       ),
		scale_value_copy( SSAP_NUM_EQUIVS_SCALING, numeric_cast<double>( arg_ssap.get_num_equivs() ) ),
		scale_value_copy( SSAP_OVERLAP_PC_SCALING,                       arg_ssap.get_overlap_pc()   ),
		scale_value_copy( SSAP_RMSD_SCALING,                             arg_ssap.get_rmsd()         ),
		scale_value_copy( SSAP_SEQ_ID_PC_SCALING,                        arg_ssap.get_seq_id_pc()    ),
		scale_value_copy( SSAP_SCORE_SCALING,                            arg_ssap.get_ssap_score()   )
	);
}

/// \brief Get the SVM the this RBF model assigns to the input vector
double rbf_model::get_score(const double_octuple &arg_to_be_scored ///< The vector to be scored
                            ) const {
	return accumulate(
		support_vectors,
		0.0,
		[&] (const double &x, const tuple<double, double_octuple> &support_vector_entry) {
			const double_octuple &sup_vec = get<1>( support_vector_entry );
			return x + (
					get<0>( support_vector_entry ) * ::std::exp(
					0.0 - gamma * (
						   ( ( get<0>( arg_to_be_scored ) - get<0>( sup_vec ) ) * ( get<0>( arg_to_be_scored ) - get<0>( sup_vec ) ) )
						 + ( ( get<1>( arg_to_be_scored ) - get<1>( sup_vec ) ) * ( get<1>( arg_to_be_scored ) - get<1>( sup_vec ) ) )
						 + ( ( get<2>( arg_to_be_scored ) - get<2>( sup_vec ) ) * ( get<2>( arg_to_be_scored ) - get<2>( sup_vec ) ) )
						 + ( ( get<3>( arg_to_be_scored ) - get<3>( sup_vec ) ) * ( get<3>( arg_to_be_scored ) - get<3>( sup_vec ) ) )
						 + ( ( get<4>( arg_to_be_scored ) - get<4>( sup_vec ) ) * ( get<4>( arg_to_be_scored ) - get<4>( sup_vec ) ) )
						 + ( ( get<5>( arg_to_be_scored ) - get<5>( sup_vec ) ) * ( get<5>( arg_to_be_scored ) - get<5>( sup_vec ) ) )
						 + ( ( get<6>( arg_to_be_scored ) - get<6>( sup_vec ) ) * ( get<6>( arg_to_be_scored ) - get<6>( sup_vec ) ) )
						 + ( ( get<7>( arg_to_be_scored ) - get<7>( sup_vec ) ) * ( get<7>( arg_to_be_scored ) - get<7>( sup_vec ) ) )
					)
				)
			);
		}
	) - b;
}

/// \brief Parse an rbf_model from an SVM-light RBF model file
///
/// \relates rbf_model
rbf_model cath::score::parse_rbf_model(istream &arg_model_istream ///< The istream containing the SVM-light RBF model file data
                                       ) {
	// Prepare some constants that are useful in the parsing
	constexpr size_t NUM_VECTOR_COMPONENTS       = 8;
	constexpr size_t NUM_COMP_PARTS              = 2;
	const     string svmlight_header_string      = "SVM-light Version V";
	const     string gamma_comment_string        = "kernel parameter -g";
	const     string b_comment_string            = "threshold b";
	const     string support_vecs_comment_string = "each following line is a SV";

	// Check that the first line is the expected header, else throw an exception
	string line_string;
	getline( arg_model_istream, line_string );
	if ( ! icontains( line_string, svmlight_header_string ) ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to parse data as SVM-light model because it doesn't begin with an SVM-light header"));
	}

	// Attempt to parse the rest of the header lines and extract the gamma and b parameters
	optional<double> parsed_gamma;
	optional<double> parsed_b;
	while ( getline( arg_model_istream, line_string ) ) {
		if ( contains( line_string, gamma_comment_string ) ) {
			parsed_gamma = stod( split_build<str_vec>( line_string, is_space(), token_compress_on ).front() );
		}
		if ( contains( line_string, b_comment_string ) ) {
			parsed_b     = stod( split_build<str_vec>( line_string, is_space(), token_compress_on ).front() );
		}
		if ( contains( line_string, support_vecs_comment_string ) ) {
			break;
		}
	}

	// If the gamma or b parameters weren't successfully parsed, then throw an exception
	if ( ! parsed_gamma || ! parsed_b ) {
		BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to parse data as SVM-light model because unable to find gamma and/or b parameters"));
	}

	// Parse the remaining lines as support vectors
	vector<pair<double, double_octuple> > support_vectors;
	while ( getline( arg_model_istream, line_string ) ) {

		// Split the line on whitespace and if there aren't at least NUM_VECTOR_COMPONENTS + 1 parts, throw an exception
		const auto line_parts = split_build<str_vec>( line_string, is_space(), token_compress_on );
		if ( line_parts.size() < NUM_VECTOR_COMPONENTS + 1 ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to parse data as SVM-light model because cannot find enough components in support vector"));
		}

		// Extract the weight from the first part
		const auto weight  = stod( line_parts.front() );

		// Parse the next NUM_VECTOR_COMPONENTS parts of the line into vectors,
		// assuming they're preceded with 1:, 2:, 3:, ..., :8 respectively
		const auto sup_vec = transform_build<doub_vec>(
			irange( 1_z, NUM_VECTOR_COMPONENTS + 1 ),
			[&] (const size_t &x) {
				const auto line_part = line_parts[ x ];
				const auto comp_pair = split_build<str_vec>( line_part, is_any_of( ":" ) );
				if ( comp_pair.size() != NUM_COMP_PARTS ) {
					BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to parse data as SVM-light model because support vector component doesn't contain two (colon-separated) parts"));
				}
				if ( stoul( comp_pair.front() ) != x ) {
					BOOST_THROW_EXCEPTION(runtime_error_exception("Unable to parse data as SVM-light model because component index is not as expected"));
				}
				return stod( comp_pair.back() );
			}
		);

		// Add the parsed data from this line to the support vectors
		support_vectors.emplace_back(
			weight,
			make_tuple(
				sup_vec[ 0 ],
				sup_vec[ 1 ],
				sup_vec[ 2 ],
				sup_vec[ 3 ],
				sup_vec[ 4 ],
				sup_vec[ 5 ],
				sup_vec[ 6 ],
				sup_vec[ 7 ]
			)
		);
	}

	// Return an rbf_model constructed from the parsed data
	return {
		*parsed_gamma,
		*parsed_b,
		support_vectors
	};
}

/// \brief Parse an rbf_model from an SVM-light RBF model file
///
/// \relates rbf_model
rbf_model cath::score::parse_rbf_model(const path &arg_svmlight_model_file ///< The SVM-light RBF model file
                                       ) {
	ifstream model_stream;
	open_ifstream( model_stream, arg_svmlight_model_file );
	const auto the_model = parse_rbf_model( model_stream );
	model_stream.close();
	return the_model;
}

/// \brief Get the score that the specified SVM assigns to the specified ssap_and_prc result
///
/// \relates rbf_model
double cath::score::get_score(const rbf_model    &arg_rbf_model,    ///< The SVM RBF model with which to score the SSAP/PRC data
                              const ssap_and_prc &arg_ssaps_and_prc ///< The SSAP/PRC data to score
                              ) {
	return arg_rbf_model.get_score( rbf_model::make_standard_scores( arg_ssaps_and_prc ) );
}

/// \brief Get the score that the specified SVM assigns to the specified SSAP scores and PRC scores
///
/// \relates rbf_model
double cath::score::get_score(const rbf_model         &arg_rbf_model, ///< The SVM RBF model with which to score the SSAP/PRC data
                              const prc_scores_entry  &arg_prc,       ///< The PRC data to score
                              const ssap_scores_entry &arg_ssap       ///< The SSAP data to score
                              ) {
	return arg_rbf_model.get_score( rbf_model::make_standard_scores( arg_prc, arg_ssap ) );
}

