/// \file
/// \brief The hmmsearch_parser class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_FILE_DETAIL_HMMSEARCH_PARSER_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_FILE_DETAIL_HMMSEARCH_PARSER_H

#include <boost/algorithm/string/predicate.hpp>

#include "common/string/string_parse_tools.hpp"
#include "exception/runtime_error_exception.hpp"
#include "resolve_hits/file/cath_id_score_category.hpp"
#include "resolve_hits/file/detail/hmmsearch_aln.hpp"
#include "resolve_hits/read_and_process_hits/read_and_process_mgr.hpp"
#include "resolve_hits/resolve_hits_type_aliases.hpp"

#include <string>
#include <regex>

namespace cath {
	namespace rslv {
		namespace detail {

			/// \brief Contain data to represent an hmmsearch output summary line for a hit
			struct hmmsearch_summary final {
				/// \brief The hit's bitscore
				double   bitscore;

				/// \brief The hit's "env" from residue index
				residx_t env_from;

				/// \brief The hit's "env" to residue index
				residx_t env_to;

				/// \brief The hit's "ali" from residue index
				residx_t ali_from;

				/// \brief The hit's "ali" to residue index
				residx_t ali_to;

				/// \brief Whether the conditional and independent evalues are "suspicious" under the CATH-Gene3D rules
				///        (see `cath-resolve-hits --cath-rules-help`)
				bool     evalues_are_susp;
			};

			/// Type alias for a vector of hmmsearch_summary entries
			using hmmsearch_summary_vec = std::vector<hmmsearch_summary>;

			/// \brief Parse hmmsearch output data
			class hmmsearch_parser final {
			private:
				/// \brief The current line (allows memory to be reused between calls to getline)
				std::string line;

				/// \brief (A reference to) the istream from which the data should be read
				std::reference_wrapper<std::istream> the_istream;

				/// \brief The current query ID
				std::string query_id;

				/// \brief A regular expression for detecting summary lines
				std::regex is_summary_line_regex{ R"(^\s*\d+\s+[\!\?])" };

				/// \brief The result of parsing any summaries
				hmmsearch_summary_vec summaries;

				/// \brief A counter for working through summary lines whilst parsing the following alignments
				size_t summary_ctr = 0;

				/// \brief The current alignment being parsed
				hmmsearch_aln the_aln;

				/// \brief Whether a hit has already been skipped for having a negative bitscore
				bool skipped_for_negtv_bitscore = false;

			public:
				explicit hmmsearch_parser(std::istream &);
				hmmsearch_parser(const std::istream &&) = delete;

				bool end_of_istream() const;

				bool line_is_empty() const;
				bool line_is_at_block() const;
				bool line_is_at_aln() const;
				bool line_is_at_pipeline_stats() const;
				bool line_is_summary_header() const;
				bool line_is_summary() const;

				bool alignment_is_empty() const;

				void advance_line();
				void advance_line_to_nonempty();
				void advance_line_to_block();
				void advance_line_until_next_aln();

				void parse_summary_from_header(const bool &);
				void parse_alignment_section();
				void finish_alignment(read_and_process_mgr &,
				                      const bool &,
				                      const residx_t &,
				                      const bool &);

				std::string & get_line();
				const std::string & get_line() const;

				std::string & get_query_id();
				const std::string & get_query_id() const;

				hmmsearch_summary_vec & get_summaries();
				const hmmsearch_summary_vec & get_summaries() const;

				static constexpr size_t LINE_BITSCORE_OFFSET    =  2;
				static constexpr size_t LINE_COND_EVALUE_OFFSET =  4;
				static constexpr size_t LINE_INDP_EVALUE_OFFSET =  5;
				static constexpr size_t LINE_ALI_FROM_OFFSET    =  9;
				static constexpr size_t LINE_ALI_TO_OFFSET      = 10;
				static constexpr size_t LINE_ENV_FROM_OFFSET    = 12;
				static constexpr size_t LINE_ENV_TO_OFFSET      = 13;
			};

			/// \brief Ctor from the istream from which the data should be read
			inline hmmsearch_parser::hmmsearch_parser(std::istream &arg_istream ///< The istream from which the data should be read
			                                          ) : the_istream( arg_istream ) {
			}

			/// \brief Whether the end of the input data has been reached
			inline bool hmmsearch_parser::end_of_istream() const {
				return ! the_istream.get();
			}

			/// \brief Whether the current line is empty
			inline bool hmmsearch_parser::line_is_empty() const {
				return line.empty();
			}

			/// \brief Whether the current line is at the start of a block
			inline bool hmmsearch_parser::line_is_at_block() const {
				return boost::algorithm::starts_with( line, ">> " );
			}

			/// \brief Whether the current line is at an alignment
			inline bool hmmsearch_parser::line_is_at_aln() const {
				return boost::algorithm::starts_with( line, "  == domain" );
			}

			/// \brief Whether the current line is the start of pipeline statistics
			inline bool hmmsearch_parser::line_is_at_pipeline_stats() const {
				return boost::algorithm::starts_with( line, "Internal pipeline statistics summary" );
			}

			/// \brief Whether the current line is a summary header
			inline bool hmmsearch_parser::line_is_summary_header() const {
				return boost::algorithm::contains( line, "core  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc" );
			}

			/// \brief Whether the current line is a summary line
			inline bool hmmsearch_parser::line_is_summary() const {
				return regex_search( line, is_summary_line_regex );
			}

			/// \brief Whether the parsed alignment is empty
			inline bool hmmsearch_parser::alignment_is_empty() const {
				return the_aln.empty();
			}

			/// \brief Advance one line
			inline void hmmsearch_parser::advance_line() {
				if ( ! getline( the_istream.get(), line ) ) {
					BOOST_THROW_EXCEPTION(common::runtime_error_exception("Unexpectedly hit end of hmmsearch output file"));
				}
			}

			/// \brief Advance lines until an non-empty line is reached
			inline void hmmsearch_parser::advance_line_to_nonempty() {
				if ( line_is_empty() ) {
					while ( getline( the_istream.get(), line ) && line_is_empty() ) {}
				}
			}

			/// \brief Advance lines until the start of a new block
			inline void hmmsearch_parser::advance_line_to_block() {
				if ( ! line_is_at_block() ) {
					while ( getline( the_istream.get(), line ) && ! line_is_at_block() ) {}
				}

				if ( ! end_of_istream() ) {
					const auto query_id_begin_itr = common::find_itr_before_first_non_space( next( common::cbegin( line ), 3 ), common::cend( line ) );
					const auto query_id_end_itr   = common::find_itr_before_first_space    ( query_id_begin_itr,                common::cend( line ) );

					query_id.assign( query_id_begin_itr, query_id_end_itr );
				}
			}

			/// \brief Advance lines until the start of the next alignment
			inline void hmmsearch_parser::advance_line_until_next_aln() {
				if ( ! line_is_at_aln() ) {
					while ( getline( the_istream.get(), line ) && ! line_is_at_aln() ) {}
				}
			}

			/// \brief Parse a summary from the current summary header line
			inline void hmmsearch_parser::parse_summary_from_header(const bool &arg_apply_cath_policies /// Whether to apply CATH-Gene3D policies (see `cath-resolve-hits --cath-rules-help`)
			                                                        ) {
				summaries.clear();
				summary_ctr = 0;
				advance_line();
				advance_line();

				while ( line_is_summary() ) {
					const auto bitscore_itrs    = common::find_field_itrs( line, LINE_BITSCORE_OFFSET                                                          );
					const auto cond_evalue_itrs = common::find_field_itrs( line, LINE_COND_EVALUE_OFFSET, 1 + LINE_BITSCORE_OFFSET,    bitscore_itrs.second    );
					const auto indp_evalue_itrs = common::find_field_itrs( line, LINE_INDP_EVALUE_OFFSET, 1 + LINE_COND_EVALUE_OFFSET, cond_evalue_itrs.second );
					const auto ali_from_itrs    = common::find_field_itrs( line, LINE_ALI_FROM_OFFSET,    1 + LINE_INDP_EVALUE_OFFSET, indp_evalue_itrs.second );
					const auto ali_to_itrs      = common::find_field_itrs( line, LINE_ALI_TO_OFFSET,      1 + LINE_ALI_FROM_OFFSET,    ali_from_itrs.second    );
					const auto env_from_itrs    = common::find_field_itrs( line, LINE_ENV_FROM_OFFSET,    1 + LINE_ALI_TO_OFFSET,      ali_to_itrs.second      );
					const auto env_to_itrs      = common::find_field_itrs( line, LINE_ENV_TO_OFFSET,      1 + LINE_ENV_FROM_OFFSET,    env_from_itrs.second    );

					summaries.push_back( hmmsearch_summary{
						common::parse_double_from_field( bitscore_itrs.first, bitscore_itrs.second ),
						common::parse_uint_from_field  ( env_from_itrs.first, env_from_itrs.second ),
						common::parse_uint_from_field  ( env_to_itrs.first,   env_to_itrs.second   ),
						arg_apply_cath_policies
							? common::parse_uint_from_field ( ali_from_itrs.first, ali_from_itrs.second )
							: 0,
						arg_apply_cath_policies
							? common::parse_uint_from_field ( ali_to_itrs.first,   ali_to_itrs.second   )
							: 0,
						hmmer_evalues_are_suspicious(
							common::parse_double_from_field( cond_evalue_itrs.first, cond_evalue_itrs.second ),
							common::parse_double_from_field( indp_evalue_itrs.first, indp_evalue_itrs.second )
						)
					} );

					advance_line();
				}
			}

			/// \brief Parse the current alignment section
			inline void hmmsearch_parser::parse_alignment_section() {
				advance_line();
				the_aln.add_a( line );
				advance_line();
				advance_line();
				the_aln.add_b( line );
				advance_line();
				advance_line();
				advance_line_to_nonempty();
			}

			/// \brief Finishe the current alignment
			inline void hmmsearch_parser::finish_alignment(read_and_process_mgr &arg_read_and_process_mgr, ///< The read_and_process_mgr to which complete hits should be added
			                                               const bool           &arg_apply_cath_policies,  ///< Whether to apply CATH-Gene3D policies (see `cath-resolve-hits --cath-rules-help`)
			                                               const residx_t       &arg_min_gap_length,       ///< The minimum length for a gap to be considered a gap
			                                               const bool           &arg_parse_hmmsearch_aln   ///< Whether to parse the hmmsearch alignment information for outputting later
			                                               ) {
				auto            aln_results   = the_aln.process_aln( arg_min_gap_length, arg_parse_hmmsearch_aln );
				std::string    &id_a          = std::get<0>( aln_results );
				hit_seg_vec    &segs          = std::get<1>( aln_results );
				auto            alnd_rngs_opt = boost::make_optional(
					arg_parse_hmmsearch_aln,
					std::get<2>( aln_results )
				);

				const auto     &summ          = summaries[ summary_ctr ];

				if ( summ.bitscore <= 0 ) {
					if ( ! skipped_for_negtv_bitscore ) {
						BOOST_LOG_TRIVIAL( warning ) << "Skipping at least one hit (eg between \""
							<< query_id
							<< "\" and \""
							<< id_a
							<< "\" with bitscore "
							<< summ.bitscore
							<< ") for having a negative bitscore, which cannot currently be handled."
							<< " It's typically not a problem to exclude such weak hits.";
						skipped_for_negtv_bitscore = true;
					}
				}
				else {
					const auto      id_score_cat  = cath_score_category_of_id( id_a, arg_apply_cath_policies );
					const bool      apply_dc_cat  = ( id_score_cat == cath_id_score_category::DC_TYPE );
					const residx_t &start         = apply_dc_cat ? summ.ali_from : summ.env_from;
					const residx_t &stop          = apply_dc_cat ? summ.ali_to   : summ.env_to;

					// If applying the CATH discontinuous policy, erase any segments after the first
					if ( apply_dc_cat ) {
						segs.erase(
							std::next( common::cbegin( segs ) ),
							           common::cend  ( segs )
						);
					}

					// Overwrite the final start/stop with those parsed from the summary
					// (env_from/env_to normally; ali_from/ali_to when applying the CATH discontinuous policy)
					segs.front().set_start_arrow( arrow_before_res( start ) );
					segs.back ().set_stop_arrow ( arrow_after_res ( stop  ) );

					arg_read_and_process_mgr.add_hit(
						query_id,
						segs,
						std::move( id_a ),
						summ.bitscore / bitscore_divisor( arg_apply_cath_policies, id_score_cat, summ.evalues_are_susp ),
						hit_score_type::BITSCORE,
						std::move( alnd_rngs_opt )
					);
				}

				++summary_ctr;
				the_aln.reset();
			}

			/// \brief Non-const getter for the current line
			inline std::string & hmmsearch_parser::get_line() {
				return line;
			}

			/// \brief Const getter for the current line
			inline const std::string & hmmsearch_parser::get_line() const {
				return line;
			}

			/// \brief Non-const getter for the query ID
			inline std::string & hmmsearch_parser::get_query_id() {
				return query_id;
			}

			/// \brief Const getter for query ID
			inline const std::string & hmmsearch_parser::get_query_id() const {
				return query_id;
			}

			/// \brief Non-const getter for the parsed summaries
			inline hmmsearch_summary_vec & hmmsearch_parser::get_summaries() {
				return summaries;
			}

			/// \brief Const getter for the parsed summaries
			inline const hmmsearch_summary_vec & hmmsearch_parser::get_summaries() const {
				return summaries;
			}

		} // namespace detail
	} // namespace rslv
} // namespace cath

#endif
