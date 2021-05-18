/// \file
/// \brief The crh_filter_spec class definitions

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

#include "crh_filter_spec.hpp"

#include "cath/common/exception/invalid_argument_exception.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::rslv;

/// \brief Getter for the worst permissible evalue before a hit is ignored
const resscr_t & crh_filter_spec::get_worst_permissible_evalue() const {
	return worst_permissible_evalue;
}

/// \brief Getter for the worst permissible cath-resolve-hits bitscore before a hit is ignored
const resscr_t & crh_filter_spec::get_worst_permissible_bitscore() const {
	return worst_permissible_bitscore;
}

/// \brief Getter for the worst permissible cath-resolve-hits score before a hit is ignored
const resscr_opt & crh_filter_spec::get_worst_permissible_score() const {
	return worst_permissible_score;
}

/// \brief Getter for the query IDs on which to filter the input, if any are present
const str_vec & crh_filter_spec::get_filter_query_ids() const {
	return filter_query_ids;
}

/// \brief Getter for the (optional) maximum number of queries to process
const size_opt & crh_filter_spec::get_limit_queries() const {
	return limit_queries;
}

/// \brief Getter for the default value for the worst permissible evalue before a hit is ignored
const doub_opt & crh_filter_spec::get_min_hmm_coverage_frac() const {
	return min_hmm_coverage_frac;
}

/// \brief Getter for the default value for the worst permissible bitscore before a hit is ignored
const doub_opt & crh_filter_spec::get_min_dc_hmm_coverage_frac() const {
	return min_dc_hmm_coverage_frac;
}

/// \brief Setter for the worst permissible evalue before a hit is ignored
crh_filter_spec & crh_filter_spec::set_worst_permissible_evalue(const resscr_t &prm_worst_permissible_evalue ///< The worst permissible evalue before a hit is ignored
                                                                ) {
	worst_permissible_evalue = prm_worst_permissible_evalue;
	return *this;
}

/// \brief Setter for the worst permissible cath-resolve-hits score before a hit is ignored
crh_filter_spec & crh_filter_spec::set_worst_permissible_bitscore(const resscr_t &prm_worst_permissible_bitscore ///< The worst permissible cath-resolve-hits score before a hit is ignored
                                                                  ) {
	worst_permissible_bitscore = prm_worst_permissible_bitscore;
	return *this;
}

/// \brief Setter for the worst permissible cath-resolve-hits score before a hit is ignored
crh_filter_spec & crh_filter_spec::set_worst_permissible_score(const resscr_opt &prm_worst_permissible_score ///< The worst permissible cath-resolve-hits score before a hit is ignored
                                                               ) {
	worst_permissible_score = prm_worst_permissible_score;
	return *this;
}

/// \brief Setter for the query IDs on which to filter the input, if any are present
crh_filter_spec & crh_filter_spec::set_filter_query_ids(const str_vec &prm_filter_query_ids ///< The query IDs on which to filter the input, if any are present
                                                        ) {
	filter_query_ids = prm_filter_query_ids;
	return *this;
}

/// \brief Setter for the (optional) maximum number of queries to process
crh_filter_spec & crh_filter_spec::set_limit_queries(const size_opt &prm_limit_queries ///< The (optional) maximum number of queries to process
                                                     ) {
	limit_queries = prm_limit_queries;
	return *this;
}

/// \brief Setter for the default value for the worst permissible evalue before a hit is ignored
crh_filter_spec & crh_filter_spec::set_min_hmm_coverage_frac(const doub_opt &prm_min_hmm_coverage_frac ///< The default value for the worst permissible evalue before a hit is ignored
                                                             ) {
	if ( ! ( prm_min_hmm_coverage_frac >= 0.0 && prm_min_hmm_coverage_frac <= 1.0 ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("min_hmm_coverage_frac must be between 0 and 1 (inclusive)"));
	}
	min_hmm_coverage_frac = prm_min_hmm_coverage_frac;
	return *this;
}

/// \brief Setter for the default value for the worst permissible bitscore before a hit is ignored
crh_filter_spec & crh_filter_spec::set_min_dc_hmm_coverage_frac(const doub_opt &prm_min_dc_hmm_coverage_frac ///< The default value for the worst permissible bitscore before a hit is ignored
                                                                ) {
	if ( ! ( prm_min_dc_hmm_coverage_frac >= 0.0 && prm_min_dc_hmm_coverage_frac <= 1.0 ) ) {
		BOOST_THROW_EXCEPTION(invalid_argument_exception("min_dc_hmm_coverage_frac must be between 0 and 1 (inclusive)"));
	}
	min_dc_hmm_coverage_frac = prm_min_dc_hmm_coverage_frac;
	return *this;
}

/// \brief Make a crh_filter_spec that accepts all valid input hits
///
/// \relates crh_filter_spec
crh_filter_spec cath::rslv::make_accept_all_filter_spec() {
	return crh_filter_spec{}
		.set_worst_permissible_evalue  ( 1'000'000'000'000'000.0f )
		.set_worst_permissible_bitscore( 0.0                      );
}