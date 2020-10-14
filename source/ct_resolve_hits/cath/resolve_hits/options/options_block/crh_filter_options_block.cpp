/// \file
/// \brief The crh_filter_options_block class definitions

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

#include "crh_filter_options_block.hpp"

#include <boost/format.hpp>

#include "cath/common/clone/make_uptr_clone.hpp"
#include "cath/common/program_options/prog_opt_num_range.hpp"

#include <iostream>

using namespace cath;
using namespace cath::common;
using namespace cath::opts;
using namespace cath::rslv;

using boost::format;
using boost::none;
using boost::program_options::options_description;
using boost::program_options::validation_error;
using boost::program_options::value;
using boost::program_options::variables_map;
using std::string;
using std::unique_ptr;

/// \brief The option name for the worst permissible evalue before a hit is ignored
const string crh_filter_options_block::PO_WORST_PERMISSIBLE_EVALUE   { "worst-permissible-evalue"   };

/// \brief The option name for the worst permissible bitscore before a hit is ignored
const string crh_filter_options_block::PO_WORST_PERMISSIBLE_BITSCORE { "worst-permissible-bitscore" };

/// \brief The option name for the worst permissible cath-resolve-hits score before a hit is ignored
const string crh_filter_options_block::PO_WORST_PERMISSIBLE_SCORE    { "worst-permissible-score"    };

/// \brief The option name for the query IDs on which to filter the input, if any are present
const string crh_filter_options_block::PO_FILTER_QUERY_ID            { "filter-query-id"            };

/// \brief The option name for the maximum number of query IDs to process
const string crh_filter_options_block::PO_LIMIT_QUERIES              { "limit-queries"              };

/// \brief The option name for the (optional) minimum coverage fraction of an HMM for a hit to be considered
const string crh_filter_options_block::PO_MIN_HMM_COVERAGE           { "min-hmm-coverage"           };

/// \brief The option name for the (optional) minimum coverage fraction of an HMM for a discontinuous hit (/^\dc_{32}$/) to be considered
const string crh_filter_options_block::PO_MIN_DC_HMM_COVERAGE        { "min-dc-hmm-coverage"        };

/// \brief A standard do_clone method
unique_ptr<options_block> crh_filter_options_block::do_clone() const {
	return { make_uptr_clone( *this ) };
}

/// \brief Define this block's name (used as a header for the block in the usage)
string crh_filter_options_block::do_get_block_name() const {
	return "Hit filtering";
}

/// \brief Add this block's options to the provided options_description
void crh_filter_options_block::do_add_visible_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                     const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                     ) {
	const string bitscore_varname { "<bitscore>" };
	const string evalue_varname   { "<evalue>"   };
	const string id_varname       { "<id>"       };
	const string num_varname      { "<num>"      };
	const string score_varname    { "<score>"    };

	const auto worst_permissible_evalue_notifier   = [&] (const resscr_t &x) { the_spec.set_worst_permissible_evalue  ( x ); };
	const auto worst_permissible_bitscore_notifier = [&] (const resscr_t &x) { the_spec.set_worst_permissible_bitscore( x ); };
	const auto worst_permissible_score_notifier    = [&] (const resscr_t &x) { the_spec.set_worst_permissible_score   ( x ); };
	const auto filter_query_ids_notifier           = [&] (const str_vec  &x) { the_spec.set_filter_query_ids          ( x ); };
	const auto limit_queries_notifier              = [&] (const size_t   &x) { the_spec.set_limit_queries             ( x ); };

	prm_desc.add_options()
		(
			( PO_WORST_PERMISSIBLE_EVALUE ).c_str(),
			value< prog_opt_num_range<resscr_t, 0, 1'000'000'000'000'000 > >()
				->value_name    ( evalue_varname                                      )
				->notifier      ( worst_permissible_evalue_notifier                   )
				->default_value (
					crh_filter_spec::DEFAULT_WORST_PERMISSIBLE_EVALUE,
					( format( "%.2g" ) % crh_filter_spec::DEFAULT_WORST_PERMISSIBLE_EVALUE ).str()
				),
			( "Ignore any hits with an evalue worse than " + evalue_varname ).c_str()
		)
		(
			( PO_WORST_PERMISSIBLE_BITSCORE ).c_str(),
			value< prog_opt_num_range<resscr_t, 0, 1'000'000'000'000'000 > >()
				->value_name    ( bitscore_varname                                    )
				->notifier      ( worst_permissible_bitscore_notifier                 )
				->default_value ( crh_filter_spec::DEFAULT_WORST_PERMISSIBLE_BITSCORE ),
			( "Ignore any hits with a bitscore worse than " + bitscore_varname ).c_str()
		)
		(
			( PO_WORST_PERMISSIBLE_SCORE ).c_str(),
			value< prog_opt_num_range<resscr_t, 0, 1'000'000'000'000'000 > >()
				->value_name    ( score_varname                                       )
				->notifier      ( worst_permissible_score_notifier                    ),
			( "Ignore any hits with a score worse than " + score_varname ).c_str()
		)
		(
			( PO_FILTER_QUERY_ID ).c_str(),
			value<str_vec>()
				->value_name    ( id_varname                                          )
				->notifier      ( filter_query_ids_notifier                           ),
			( "Ignore all input data except that for query protein(s) " + id_varname
			  + "\n(may be specified multiple times for multiple query proteins)" ).c_str()
		)
		(
			( PO_LIMIT_QUERIES ).c_str(),
			value<size_t>()
				->value_name    ( num_varname                                         )
				->notifier      ( limit_queries_notifier                              )
				->implicit_value( 1                                                   ),
			( "Only process the first " + num_varname
			  + " query protein(s) encountered in the input data" ).c_str()
		);
}

/// \brief Add a hidden option to the options_description for the hmm coverage options
void crh_filter_options_block::do_add_hidden_options_to_description(options_description &prm_desc,           ///< The options_description to which the options are added
                                                                    const size_t        &/*prm_line_length*/ ///< The line length to be used when outputting the description (not very clearly documented in Boost)
                                                                    ) {
	const string percent_varname{ "<percent>" };

	const auto check_pc_fn = [] (const double &x, const string &y) {
		using std::to_string;
		if ( ! ( x >= 0.0 && x <= 100.0 ) ) {
			throw validation_error{
				validation_error::invalid_option_value,
				y,
				to_string( x )
			};
		}
	};

	const auto min_hmm_cvg_ntfr    = [&] (const double &x) { check_pc_fn( x, PO_MIN_HMM_COVERAGE    ); the_spec.set_min_hmm_coverage_frac   ( x / 100.0 ); };
	const auto min_dc_hmm_cvg_ntfr = [&] (const double &x) { check_pc_fn( x, PO_MIN_DC_HMM_COVERAGE ); the_spec.set_min_dc_hmm_coverage_frac( x / 100.0 ); };

	prm_desc.add_options()
		(
			PO_MIN_HMM_COVERAGE.c_str(),
			value<double>()
				->value_name    ( percent_varname     )
				->notifier      ( min_hmm_cvg_ntfr    )
				->implicit_value( 50.0                ),
			( "[IN PUBLIC GENE3D COMMAND] In hmmsearch_out input, ignore any hits for which 100.0 * ( hmm_to + 1 - hmm_from ) / hmm_length  <  "
				+ percent_varname ).c_str()
		)
		(
			PO_MIN_DC_HMM_COVERAGE.c_str(),
			value<double>()
				->value_name    ( percent_varname     )
				->notifier      ( min_dc_hmm_cvg_ntfr )
				->implicit_value( 80.0                ),
			( R"([IN PUBLIC GENE3D COMMAND] In hmmsearch_out input, ignore any /^\dc_{32}$/ hits for which 100.0 * ( hmm_to + 1 - hmm_from ) / hmm_length  <  )"
				+ percent_varname + "\n"
				R"((overriding any --worst-hmm-coverage value for those /^\dc_{32}$/ hits))" ).c_str()
		);
}

/// \brief Generate a description of any problem that makes the specified crh_filter_options_block invalid
///        or none otherwise
str_opt crh_filter_options_block::do_invalid_string(const variables_map &/*prm_variables_map*/ ///< The variables map, which options_blocks can use to determine which options were specified, defaulted etc
                                                    ) const {
	if (the_spec.get_limit_queries() && ! the_spec.get_filter_query_ids().empty() ) {
		return "Cannot specify both --" + PO_FILTER_QUERY_ID + " and --" + PO_LIMIT_QUERIES;
	}
	return none;
}

/// \brief Return all options names for this block
str_vec crh_filter_options_block::do_get_all_options_names() const {
	return {
		crh_filter_options_block::PO_WORST_PERMISSIBLE_EVALUE,
		crh_filter_options_block::PO_WORST_PERMISSIBLE_BITSCORE,
		crh_filter_options_block::PO_WORST_PERMISSIBLE_SCORE,
		crh_filter_options_block::PO_FILTER_QUERY_ID,
		crh_filter_options_block::PO_LIMIT_QUERIES,
	};
}

/// \brief Getter for the crh_filter_spec that the crh_filter_options_block configures
const crh_filter_spec & crh_filter_options_block::get_crh_filter_spec() const {
	return the_spec;
}
