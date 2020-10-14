/// \file
/// \brief The hmmer_parser class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_FILE_DETAIL_HMMER_PARSER_HPP
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_FILE_DETAIL_HMMER_PARSER_HPP

#include <boost/algorithm/string/predicate.hpp>
#include <boost/optional/optional_io.hpp>

#include "cath/common/boost_addenda/make_string_ref.hpp"
#include "cath/common/exception/out_of_range_exception.hpp"
#include "cath/common/exception/runtime_error_exception.hpp"
#include "cath/common/string/string_parse_tools.hpp"
#include "cath/resolve_hits/file/cath_id_score_category.hpp"
#include "cath/resolve_hits/file/detail/hmmer_aln.hpp"
#include "cath/resolve_hits/file/hmmer_format.hpp"
#include "cath/resolve_hits/read_and_process_hits/read_and_process_mgr.hpp"
#include "cath/resolve_hits/resolve_hits_type_aliases.hpp"

#include <regex>
#include <string>

namespace cath {
	namespace rslv {
		namespace detail {

			/// \brief Whether the hmm coverage is OK or has failed a required threshold
			///
			/// The hmm threshold is the fraction of the HMM involved in the hmmsearch hit
			/// (`( hmm_to + 1 - hmm_from ) / hmm_length`)
			enum class hmm_coverage : bool {
				OK,     ///< The hmm coverage is OK
				TOO_LOW ///< The hmm coverage has failed a required threshold
			};

			/// \brief Contain data to represent a HMMER output summary line for a hit
			struct hmmer_summary final {
				/// \brief The hit's bitscore
				double        bitscore;

				/// \brief The hit's "env" from residue index
				seq::residx_t env_from;

				/// \brief The hit's "env" to residue index
				seq::residx_t env_to;

				/// \brief The hit's "ali" from residue index
				seq::residx_t ali_from;

				/// \brief The hit's "ali" to residue index
				seq::residx_t ali_to;

				/// \brief Whether the conditional and independent evalues are "suspicious" under the CATH-Gene3D rules
				///        (see `cath-resolve-hits --cath-rules-help`)
				bool          evalues_are_susp;

				/// \brief The hit's conditional evalue
				double        conditional_evalue;

				/// \brief The hit's independent evalue
				double        independent_evalue;

				/// \brief Whether the hmm coverage for the hit is OK
				hmm_coverage  hmm_coverage_is_ok;
			};

			/// Type alias for a vector of hmmer_summary entries
			using hmmer_summary_vec = std::vector<hmmer_summary>;

			/// \brief Parse HMMER output data
			class hmmer_parser final {
			private:
				/// \brief The HMMER format to parse
				hmmer_format format;

				/// \brief The current line (allows memory to be reused between calls to getline)
				std::string line;

				/// \brief (A reference to) the istream from which the data should be read
				std::reference_wrapper<std::istream> the_istream;

				/// \brief The current query ID
				str_opt query_id;

				/// \brief A regular expression for detecting summary lines
				std::regex is_summary_line_regex{ R"(^\s*\d+\s+[\!\?])" };

				/// \brief The result of parsing any summaries
				hmmer_summary_vec summaries;

				/// \brief A counter for working through summary lines whilst parsing the following alignments
				size_t summary_ctr = 0;

				/// \brief The current alignment being parsed
				hmmer_aln the_aln;

				/// \brief Whether a hit has already been skipped for having a negative bitscore
				bool skipped_for_negtv_bitscore = false;

				/// \brief The match ID parsed from the last prefix line (/^Query: /) if any
				str_opt  prefix_match_id;

				/// \brief The hmm length parsed from the last prefix line (/^Query: /) if any
				size_opt prefix_hmm_length;

				template <typename T>
				inline static auto get_query_id_impl(T &prm_hmmer_parser) -> decltype( *prm_hmmer_parser.query_id );

				void advance_to_block_or_prefix();
				bool line_is_at_prefix() const;
				bool line_is_empty() const;
				bool line_is_summary() const;

				str_opt & get_prefix_id();
				str_opt & get_block_id();

			public:
				hmmer_parser(const hmmer_format &,
				             std::istream &);
				hmmer_parser(const hmmer_format &,
				             const std::istream &&) = delete;

				bool end_of_istream() const;

				bool line_is_at_block() const;
				bool line_is_at_aln() const;
				bool line_is_at_pipeline_stats() const;
				bool line_is_summary_header() const;

				bool alignment_is_empty() const;

				void advance_line();
				void advance_line_to_nonempty();
				void advance_line_to_block();
				void advance_line_until_next_aln();

				void parse_summary_from_header(const bool &,
				                               const crh_filter_spec &);
				void parse_alignment_section();
				void finish_alignment(read_and_process_mgr &,
				                      const bool &,
				                      const seq::residx_t &,
				                      const bool &);

				std::string & get_line();
				const std::string & get_line() const;

				std::string & get_query_id();
				const std::string & get_query_id() const;

				hmmer_summary_vec & get_summaries();
				const hmmer_summary_vec & get_summaries() const;

				static constexpr size_t LINE_BITSCORE_OFFSET    =  2;
				static constexpr size_t LINE_COND_EVALUE_OFFSET =  4;
				static constexpr size_t LINE_INDP_EVALUE_OFFSET =  5;
				static constexpr size_t LINE_HMM_FROM_OFFSET    =  6;
				static constexpr size_t LINE_HMM_TO_OFFSET      =  7;
				static constexpr size_t LINE_ALI_FROM_OFFSET    =  9;
				static constexpr size_t LINE_ALI_TO_OFFSET      = 10;
				static constexpr size_t LINE_ENV_FROM_OFFSET    = 12;
				static constexpr size_t LINE_ENV_TO_OFFSET      = 13;
			};

			/// \brief const-agnostic implementation of get_query_id()
			///
			/// See GSL rule: Pro.Type.3: Don't use const_cast to cast away const (i.e., at all)
			/// (https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md#Pro-type-constcast)
			template <typename T>
			inline auto hmmer_parser::get_query_id_impl(T &prm_hmmer_parser ///< The (const / non-const) hmmer_parser on which the get_query_id() is being called
			                                            ) -> decltype( *prm_hmmer_parser.query_id ) {
				/* the complex logic around getting a possibly-const reference to my_bar */
				if ( ! prm_hmmer_parser.query_id ) {
					BOOST_THROW_EXCEPTION(common::runtime_error_exception(
						"Whilst parsing "
						+ ( prm_hmmer_parser.format == hmmer_format::HMMSEARCH ? "hmmsearch"s : "hmmscan"s )
						+ " input, was unable to retrieve a query ID when required."
						+ (
							prm_hmmer_parser.prefix_match_id
							?
								" Another ID ("
								+ *prm_hmmer_parser.prefix_match_id
								+ ") was present - perhaps this is actually "
								+ ( prm_hmmer_parser.format != hmmer_format::HMMSEARCH ? "hmmsearch"s : "hmmscan"s )
								+ " output?"
							: ""s
						)
					));
				}
				return *prm_hmmer_parser.query_id;
			}


			/// \brief Advance the istream to the next block or prefix
			inline void hmmer_parser::advance_to_block_or_prefix() {
				if ( ! line_is_at_prefix() && ! line_is_at_block() ) {
					while (
						getline( the_istream.get(), line )
						&&
						! line_is_at_prefix()
						&&
						! line_is_at_block()
						&&
						! end_of_istream()
					) {}
				}
			}

			/// \brief Whether the current line is at a prefix
			inline bool hmmer_parser::line_is_at_prefix() const {
				return boost::algorithm::starts_with( line, "Query: " );
			}

			/// \brief Whether the current line is empty
			inline bool hmmer_parser::line_is_empty() const {
				return line.empty();
			}

			/// \brief Whether the current line is a summary line
			inline bool hmmer_parser::line_is_summary() const {
				return regex_search( line, is_summary_line_regex );
			}

			/// \brief Get whichever ID is associated with the prefix
			///
			/// This is different when parsing hmmsearch / hmmscan output
			inline str_opt & hmmer_parser::get_prefix_id() {
				return ( format == hmmer_format::HMMSEARCH ) ? prefix_match_id : query_id;
			}

			/// \brief Get whichever ID is associated with the block
			///
			/// This is different when parsing hmmsearch / hmmscan output
			inline str_opt & hmmer_parser::get_block_id() {
				return ( format == hmmer_format::HMMSEARCH ) ? query_id        : prefix_match_id;
			}

			/// \brief Ctor from the istream from which the data should be read
			inline hmmer_parser::hmmer_parser(const hmmer_format &prm_format, ///< The HMMER format to parse
			                                  std::istream       &prm_istream ///< The istream from which the data should be read
			                                  ) : format     { prm_format  },
			                                      the_istream{ prm_istream } {
			}

			/// \brief Whether the end of the input data has been reached
			inline bool hmmer_parser::end_of_istream() const {
				return ! the_istream.get();
			}

			/// \brief Whether the current line is at the start of a block
			inline bool hmmer_parser::line_is_at_block() const {
				return boost::algorithm::starts_with( line, ">> " );
			}

			/// \brief Whether the current line is at an alignment
			inline bool hmmer_parser::line_is_at_aln() const {
				return boost::algorithm::starts_with( line, "  == domain" );
			}

			/// \brief Whether the current line is the start of pipeline statistics
			inline bool hmmer_parser::line_is_at_pipeline_stats() const {
				return boost::algorithm::starts_with( line, "Internal pipeline statistics summary" );
			}

			/// \brief Whether the current line is a summary header
			inline bool hmmer_parser::line_is_summary_header() const {
				return boost::algorithm::contains( line, "core  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc" );
			}

			/// \brief Whether the parsed alignment is empty
			inline bool hmmer_parser::alignment_is_empty() const {
				return the_aln.empty();
			}

			/// \brief Advance one line
			inline void hmmer_parser::advance_line() {
				if ( ! getline( the_istream.get(), line ) ) {
					BOOST_THROW_EXCEPTION(common::runtime_error_exception("Unexpectedly hit end of hmmsearch/hmmscan output file"));
				}
			}

			/// \brief Advance lines until an non-empty line is reached
			inline void hmmer_parser::advance_line_to_nonempty() {
				if ( line_is_empty() ) {
					while ( getline( the_istream.get(), line ) && line_is_empty() ) {}
				}
			}

			/// \brief Advance lines until the start of a new block
			inline void hmmer_parser::advance_line_to_block() {
				prefix_hmm_length = boost::none;
				prefix_match_id   = boost::none;

				while ( ! line_is_at_block() && ! end_of_istream() ) {
					advance_to_block_or_prefix();

					if ( line_is_at_prefix() && ! end_of_istream() ) {
						const auto line_end = common::cend( line );
						const auto match_id_begin_itr = common::find_itr_before_first_non_space( next( common::cbegin( line ), 7 ), line_end );
						const auto match_id_end_itr   = common::find_itr_before_first_space    ( match_id_begin_itr,                line_end );
						const auto length_pre_itr     = common::find_itr_before_first_non_space( match_id_end_itr,                  line_end );

						if (
							distance( length_pre_itr, line_end ) < 5
							||
							*std::prev( line_end ) != ']'
							||
							(
								! boost::algorithm::starts_with( common::make_string_ref( length_pre_itr, line_end ), "[M=" )
								&&
								! boost::algorithm::starts_with( common::make_string_ref( length_pre_itr, line_end ), "[L=" )
							)

						) {
							BOOST_THROW_EXCEPTION(common::out_of_range_exception(
								"Cannot parse HMM length out of hmmsearch/hmmscan line \""
								+ line
								+ "\""
							));
						}

						get_prefix_id().emplace( match_id_begin_itr, match_id_end_itr );
						prefix_hmm_length = common::parse_uint_from_field( next( length_pre_itr, 3 ), prev( line_end ) );

						getline( the_istream.get(), line );
					}
				}

				if ( ! end_of_istream() ) {
					const auto query_id_begin_itr = common::find_itr_before_first_non_space( next( common::cbegin( line ), 3 ), common::cend( line ) );
					const auto query_id_end_itr   = common::find_itr_before_first_space    ( query_id_begin_itr,                common::cend( line ) );

					get_block_id().emplace( query_id_begin_itr, query_id_end_itr );
				}
			}

			/// \brief Advance lines until the start of the next alignment
			inline void hmmer_parser::advance_line_until_next_aln() {
				if ( ! line_is_at_aln() ) {
					while ( getline( the_istream.get(), line ) && ! line_is_at_aln() ) {}
				}
			}

			/// \brief Parse a summary from the current summary header line
			inline void hmmer_parser::parse_summary_from_header(const bool            &prm_apply_cath_policies, ///< Whether to apply CATH-Gene3D policies (see `cath-resolve-hits --cath-rules-help`)
			                                                    const crh_filter_spec &prm_filter_spec          ///< The filter spec, used to determine which hits to exclude on inadequate hmm coverage
			                                                    ) {
				// Get the min hmm coverage for this specific match ID
				const doub_opt min_hmm_coverage = static_cast<bool>( prefix_match_id )
					? hmm_coverage_for_match( prm_filter_spec, *prefix_match_id )
					: boost::none;

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

					const double conditional_evalue = common::parse_double_from_field( cond_evalue_itrs.first, cond_evalue_itrs.second );
					const double independent_evalue = common::parse_double_from_field( indp_evalue_itrs.first, indp_evalue_itrs.second );

					// Calculate the hmm coverage status
					const hmm_coverage hmm_coverage_status = [&] {
						if ( min_hmm_coverage && prefix_hmm_length ) {
							const auto hmm_from_itrs = common::find_field_itrs( line, LINE_HMM_FROM_OFFSET, 1 + LINE_INDP_EVALUE_OFFSET, indp_evalue_itrs.second );
							const auto hmm_to_itrs   = common::find_field_itrs( line, LINE_HMM_TO_OFFSET,   1 + LINE_HMM_FROM_OFFSET,    hmm_from_itrs.second    );
							const auto hmm_from      = common::parse_double_from_field( hmm_from_itrs.first, hmm_from_itrs.second );
							const auto hmm_to        = common::parse_double_from_field( hmm_to_itrs.first,   hmm_to_itrs.second   );
							const auto hmm_coverage  = (
								( hmm_to + 1.0 - hmm_from ) / debug_numeric_cast<double>( *prefix_hmm_length )
							);
							if ( hmm_coverage < *min_hmm_coverage ) {
								return detail::hmm_coverage::TOO_LOW;
							}
						}
						return detail::hmm_coverage::OK;
					} ();

					summaries.push_back( hmmer_summary{
						common::parse_double_from_field( bitscore_itrs.first, bitscore_itrs.second ),
						common::parse_uint_from_field  ( env_from_itrs.first, env_from_itrs.second ),
						common::parse_uint_from_field  ( env_to_itrs.first,   env_to_itrs.second   ),
						prm_apply_cath_policies
							? common::parse_uint_from_field ( ali_from_itrs.first, ali_from_itrs.second )
							: 0,
						prm_apply_cath_policies
							? common::parse_uint_from_field ( ali_to_itrs.first,   ali_to_itrs.second   )
							: 0,
						hmmer_evalues_are_suspicious(
							conditional_evalue,
							independent_evalue
						),
						conditional_evalue,
						independent_evalue,
						hmm_coverage_status
					} );

					advance_line();
				}
			}

			/// \brief Parse the current alignment section
			inline void hmmer_parser::parse_alignment_section() {
				if ( boost::algorithm::ends_with( line, " RF" ) ) {
					advance_line();
				}
				the_aln.add_a( line );
				advance_line();
				advance_line();
				the_aln.add_b( line );
				advance_line();
				advance_line();
				advance_line_to_nonempty();
			}

			/// \brief Finish the current alignment
			inline void hmmer_parser::finish_alignment(read_and_process_mgr &prm_read_and_process_mgr, ///< The read_and_process_mgr to which complete hits should be added
			                                           const bool           &prm_apply_cath_policies,  ///< Whether to apply CATH-Gene3D policies (see `cath-resolve-hits --cath-rules-help`)
			                                           const seq::residx_t  &prm_min_gap_length,       ///< The minimum length for a gap to be considered a gap
			                                           const bool           &prm_parse_hmmer_aln       ///< Whether to parse the HMMER alignment information for outputting later
			                                           ) {
				auto              aln_results   = the_aln.process_aln( prm_min_gap_length, prm_parse_hmmer_aln );
				std::string      &id_a          = std::get<0>( aln_results );
				seq::seq_seg_vec &segs          = std::get<1>( aln_results );
				auto              alnd_rngs_opt = boost::make_optional(
					prm_parse_hmmer_aln,
					std::get<2>( aln_results )
				);

				if ( prefix_match_id && *prefix_match_id != id_a ) {
					BOOST_THROW_EXCEPTION(common::runtime_error_exception(
						"Whilst parsing "
						+ ( format == hmmer_format::HMMSEARCH ? "hmmsearch"s : "hmmscan"s )
						+ " output, found ID \""
						+ id_a
						+ "\" that mismatches previously recorded ID \""
						+ *prefix_match_id
						+ "\"."
						+ (
							( query_id && *query_id == id_a )
							?
								+ " But the other ID matches so perhaps this is actually "
								+ ( format != hmmer_format::HMMSEARCH ? "hmmsearch"s : "hmmscan"s )
								+ " output?"
							: ""
						)
					));
				}

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
				else if ( summ.hmm_coverage_is_ok == hmm_coverage::OK ) {
					const auto           id_score_cat  = cath_score_category_of_id( id_a, prm_apply_cath_policies );
					const bool           apply_dc_cat  = ( id_score_cat == cath_id_score_category::DC_TYPE );
					const seq::residx_t &start         = apply_dc_cat ? summ.ali_from : summ.env_from;
					const seq::residx_t &stop          = apply_dc_cat ? summ.ali_to   : summ.env_to;

					// If applying the CATH discontinuous policy, erase any segments after the first
					if ( apply_dc_cat ) {
						segs.erase(
							std::next( common::cbegin( segs ) ),
							           common::cend  ( segs )
						);
					}

					// Overwrite the final start/stop with those parsed from the summary
					// (env_from/env_to normally; ali_from/ali_to when applying the CATH discontinuous policy)
					segs.front().set_start_arrow( seq::arrow_before_res( start ) );
					segs.back ().set_stop_arrow ( seq::arrow_after_res ( stop  ) );

					hit_extras_store extras;
					if ( prm_parse_hmmer_aln ) {
						extras.push_back< hit_extra_cat::ALND_RGNS >( to_string( std::get<2>( aln_results ) ) );
					}
					extras.push_back< hit_extra_cat::COND_EVAL >( summ.conditional_evalue );
					extras.push_back< hit_extra_cat::INDP_EVAL >( summ.independent_evalue );

					prm_read_and_process_mgr.add_hit(
						*query_id,
						std::move( segs ),
						std::move( id_a ),
						summ.bitscore / bitscore_divisor( prm_apply_cath_policies, summ.evalues_are_susp ),
						hit_score_type::BITSCORE,
						std::move( extras )
					);
				}

				++summary_ctr;
				the_aln.reset();
			}

			/// \brief Non-const getter for the current line
			inline std::string & hmmer_parser::get_line() {
				return line;
			}

			/// \brief Const getter for the current line
			inline const std::string & hmmer_parser::get_line() const {
				return line;
			}

			/// \brief Non-const getter for the query ID
			inline std::string & hmmer_parser::get_query_id() {
				return get_query_id_impl( *this );
			}

			/// \brief Const getter for query ID
			inline const std::string & hmmer_parser::get_query_id() const {
				return get_query_id_impl( *this );
			}

			/// \brief Non-const getter for the parsed summaries
			inline hmmer_summary_vec & hmmer_parser::get_summaries() {
				return summaries;
			}

			/// \brief Const getter for the parsed summaries
			inline const hmmer_summary_vec & hmmer_parser::get_summaries() const {
				return summaries;
			}

		} // namespace detail
	} // namespace rslv
} // namespace cath

#endif
