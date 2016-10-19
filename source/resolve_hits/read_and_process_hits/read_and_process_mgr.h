/// \file
/// \brief The read_and_process_mgr class header

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

#ifndef _CATH_TOOLS_SOURCE_RESOLVE_HITS_READ_AND_PROCESS_HITS_READ_AND_PROCESS_MGR_H
#define _CATH_TOOLS_SOURCE_RESOLVE_HITS_READ_AND_PROCESS_HITS_READ_AND_PROCESS_MGR_H

#include <boost/optional.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/utility/string_ref.hpp>

#include "common/algorithm/sort_uniq_build.h"
#include "common/type_aliases.h"
#include "exception/runtime_error_exception.h"
#include "resolve_hits/calc_hit.h"
#include "resolve_hits/calc_hit_list.h"
#include "resolve_hits/options/spec/crh_filter_spec.h"
#include "resolve_hits/read_and_process_hits/hits_processor/hits_processor.h"

#include <future>
#include <unordered_map>

namespace cath { namespace rslv { class crh_input_spec; } }
namespace cath { namespace rslv { class crh_spec; } }

namespace cath {
	namespace rslv {

		// Things different managers might or might not want to do:
		//  * read in extra data (evalue(s)) / don't
		//  * resolve on way / at end /not at all
		//  * throw away data immediately after resolve / don't
		//  * put new hit in correct place immediately / don't
		//  * emplace new hit in new thread / same thread
		//  * remove strictly worse hits during emplace / don't

		/// \brief Manage the hits being read from the input file
		///
		/// This stores the hits, indexed by query_id in a unordered_map.
		///
		/// Client code should call add_hit() for each new hit and then process_all_outstanding()
		/// at the end of the input
		///
		/// If `input_hits_are_grouped` then this can trigger asynchronous processing
		/// of a block of hits for a given query_id when it reaches the end. To allow
		/// this, it stores the previous hit's query ID and when a hit with a different
		/// query ID is encountered, it triggers processing of the hits associated with
		/// the previous query ID.
		///
		/// If `! input_hits_are_grouped`, then the results must all be processed synchronously
		/// at the end.
		///
		///
		/// This class separates out the reading code from the code that processes
		/// the results as they come in.
		///
		/// Any asynchronous processing is done using in the static member function process_query_id(),
		/// which has access to:
		///  * (non-const) the relevant hits data
		///  * (non-const) the (reference-wrapped) ostream
		///
		/// (This is OK because std::unordered_map guarantees to not invalidate references/pointers to an element's key/data
		///  (unless the element is erased), even on operations that do invalidate iterators to the element)
		///
		/// The call to async takes its own copy of the query_id string. The worker thread should have access to
		/// no other data members.
		class read_and_process_mgr final {
		private:
			/// \brief (A unique_ptr to) the hits_processor that will process the hits
			///
			/// Could easily be a clone_ptr if there's any reason to want read_and_process_mgr to be copyable
			std::unique_ptr<detail::hits_processor> processor_ptr;

			/// \brief The filter spec to define how to filter the hits
			crh_filter_spec the_filter_spec;

			/// \brief A type alias for an unordered_map from the query_id string to the corresponding hits data
			using str_hit_list_umap = std::unordered_map<std::string, full_hit_list>;

			/// \brief An unordered_map from the query_id string to the corresponding hits data
			str_hit_list_umap hit_list_by_query_id;

			/// \brief Whether or not (the user guaranteees that) the input hits data is presorted
			///        to ensure all hits for a query_id are consecutive
			///
			/// When this is true, it's possible to process and erase a block of hits for a query_id when the
			/// end of the block is detected, which can be faster and more memory efficient.
			bool input_hits_are_grouped = DEFAULT_INPUT_HITS_ARE_GROUPED;

			/// \brief A string that can be reused for holding a local copy of the query ID
			///        for hashing to avoid having to reallocate on every call to add_hit()
			std::string temp_hashable_query_id;

			/// \brief A type alias for a pair of string and hit data (non-const) reference
			using str_hits_ref_pair = std::pair<std::string, full_hit_list &>;

			/// \brief The query_id of the previous result and a reference to the corresponding data (if any)
			///
			/// This is used when `input_hits_are_grouped` to indicate that a calc_hit_list under this query_id
			/// can be processed as soon as a result is found with a different query id
			///
			/// The hit data reference avoids many of the repeated calls to the hashing function when input_hits_are_grouped
			boost::optional<str_hits_ref_pair> prev_query_id_and_hits_ref;

			/// \brief A record of the query_id associated with the latest asynchronous processing job
			///        so that the thread can detect if it finds itself trying to add another hit to
			///        data that it assumed was finished and which it has passed off to an async worker
			///        thread.
			///
			/// This is only used by the main thread, never by the async worker thread.
			str_opt to_be_erased_query_id;

			/// \brief The std::future for the std::async() processing.
			///
			/// At present, this isn't used to pass back results, just to allow the main thread to wait
			/// until an async processing job is complete.
			std::future<void> resolve_future;

			static void process_query_id(detail::hits_processor &,
			                             const std::string &,
			                             const crh_filter_spec &,
			                             full_hit_list &);

			void trigger_async_process_query_id(const std::string &);

			void wait_for_any_active_work() const;

		public:
			/// \brief The default value for whether input hits can be assumed to be pre-sorted
			///
			/// This is false - better not to assume this guarantee
			static constexpr bool DEFAULT_INPUT_HITS_ARE_GROUPED = false;

			explicit read_and_process_mgr(const detail::hits_processor &,
			                              const crh_filter_spec &,
			                              const bool & = DEFAULT_INPUT_HITS_ARE_GROUPED);

			void add_hit(const boost::string_ref &,
			             const hit_seg_vec &,
			             std::string &&,
			             const double &,
			             const hit_score_type &);

			void process_all_outstanding();

			const crh_filter_spec & get_filter_spec() const;
		};

		const str_vec & get_filter_query_ids(const read_and_process_mgr &);

		bool should_skip_query_id(const read_and_process_mgr &,
		                          const std::string &);

		bool should_skip_query_id(const read_and_process_mgr &,
		                          const boost::string_ref &);

		read_and_process_mgr make_read_and_process_mgr(const detail::hits_processor &,
		                                               const crh_spec &);

		read_and_process_mgr make_read_and_process_mgr(const detail::hits_processor &,
		                                               const crh_filter_spec &,
		                                               const crh_input_spec &);

		read_and_process_mgr make_read_and_process_mgr(std::ostream &,
		                                               const crh_spec &);

		/// \brief Process the specified data
		///
		/// This is called directly in process_all_outstanding() and through async in trigger_async_process_query_id()
		inline void read_and_process_mgr::process_query_id(detail::hits_processor &arg_processor,   ///< The hits_processor to use to process the hits
		                                                   const std::string      &arg_query_id,    ///< The query ID
		                                                   const crh_filter_spec  &arg_filter_spec, ///< The filter spec to define how to filter the hits
		                                                   full_hit_list          &arg_full_hits    ///< The hits to process
		                                                   ) {
			arg_processor.process_hits_for_query( arg_query_id, arg_filter_spec, arg_full_hits );
		}

		/// \brief Trigger asynchronous processing of the data corresponding to the specified protein_query_id
		inline void read_and_process_mgr::trigger_async_process_query_id(const std::string &arg_query_id ///< The protein_query_id
		                                                                 ) {
			// Wait until any previously triggered asynchronous processing is complete
			wait_for_any_active_work();

			// Mark to_be_erased_query_id with this query ID so it's possible to detect if there's any
			// attempt to modify the data for this query ID while it's being processed
			to_be_erased_query_id = arg_query_id;

			// Start asynchronous processing of the work to process these hits and output the results
			resolve_future = async(
				std::launch::async,
				process_query_id,
				std::ref( *processor_ptr ),
				arg_query_id,
				the_filter_spec,
				std::ref( hit_list_by_query_id[ arg_query_id ] )
			);
		}

		/// \brief Check that any previously triggered async computation is complete
		inline void read_and_process_mgr::wait_for_any_active_work() const {
			if ( resolve_future.valid() ) {
				resolve_future.wait();
			}
		}

		/// \brief Ctor from the ostream to which the results should be written
		inline read_and_process_mgr::read_and_process_mgr(const detail::hits_processor &arg_hits_processor,        ///< The hits_processor to use to process the hits
		                                                  const crh_filter_spec        &arg_filter_spec,           ///< The filter spec to define how to filter the hits
		                                                  const bool                   &arg_input_hits_are_grouped ///< Whether or not the input hits are guaranteed to be presorted
		                                                  ) : processor_ptr          { arg_hits_processor.clone() },
		                                                      the_filter_spec        { arg_filter_spec            },
		                                                      input_hits_are_grouped { arg_input_hits_are_grouped } {
		}

		/// \brief Add a new hit for the current query_id
		///
		/// \pre `is_active()` else an invalid_argument_exception will be thrown
		inline void read_and_process_mgr::add_hit(const boost::string_ref &arg_query_id,  ///< A string_ref of the query_id
		                                          const hit_seg_vec       &arg_segments,  ///< Any fragments of the new hit
		                                          std::string            &&arg_label,     ///< The label associated with the new hit
		                                          const double            &arg_score,     ///< The score associated with the new hit
		                                          const hit_score_type    &arg_score_type ///< The type of the score
		                                          ) {
			// If this hit's score doesn't meet the filter and such hits don't need to be kept, then skip it
			if ( ! processor_ptr->parse_hits_that_fail_score_filter() ) {
				if ( ! score_passes_filter( the_filter_spec, arg_score, arg_score_type ) ) {
					return;
				}
			}

			// Store the query_id in a local string temp_hashable_query_id, which can then be
			// used for hashing
			temp_hashable_query_id.assign( arg_query_id.data(), arg_query_id.length() );

			// If input_hits_are_grouped and previous hits had a different query_id, trigger an async worker thread to
			// process the hits associated with that previous query_id
			if ( input_hits_are_grouped && prev_query_id_and_hits_ref && prev_query_id_and_hits_ref->first != temp_hashable_query_id ) {
				trigger_async_process_query_id( prev_query_id_and_hits_ref->first );
				prev_query_id_and_hits_ref = boost::none;
			}

			// If this query_id matches to_be_erased_query_id then something's wrong
			// (probably that the input data violates a given guarantee that it'd be grouped by query_id)
			if ( to_be_erased_query_id && *to_be_erased_query_id == temp_hashable_query_id ) {

				// ** Should probably catch this exception and exit with error code somwhere **
				BOOST_THROW_EXCEPTION(common::runtime_error_exception(
					"Attempt to add a hit for a query_id "
					+ temp_hashable_query_id
					+ " that previously appeared to be finished - suggests the input data breaks the guarantee to be grouped by query_id"
				));
			}

			// If this is the same query_id as the previous query, then just re-use the cached reference to that query's hits data
			// otherwise, look it up in hit_list_by_query_id
			auto &the_hits = prev_query_id_and_hits_ref ? prev_query_id_and_hits_ref->second
			                                            : hit_list_by_query_id[ temp_hashable_query_id ];

			// Add the new hit to the query's hits data
			the_hits.emplace_back(
				arg_segments,
				move( arg_label ),
				arg_score,
				arg_score_type
			);

			// If the input hits are presorted then ensure prev_query_id_and_hits_ref is up-to-date
			// with this hit
			if ( input_hits_are_grouped && ! prev_query_id_and_hits_ref ) {
				prev_query_id_and_hits_ref = str_hits_ref_pair{ temp_hashable_query_id, the_hits };
			}
		}

		/// \brief Process all outstanding data
		inline void read_and_process_mgr::process_all_outstanding() {
			// If input_hits_are_grouped, there may be asynchronous processing going on, so wait until that's
			// finished before processing this data (else the results may get interleaved)
			if ( input_hits_are_grouped ) {
				wait_for_any_active_work();
			}


			// Get a sorted list of all the query IDs
			const auto sorted_query_ids = common::sort_build<str_vec>(
				hit_list_by_query_id | boost::adaptors::map_keys
			);

			// Loop over the sorted query IDs
			for (const auto &query_id : sorted_query_ids) {
				// Grab the data for this query_id and process it
				auto &query_id_full_hits = hit_list_by_query_id[ query_id ];
				if ( ! query_id_full_hits.empty() ) {
					process_query_id(
						*processor_ptr,
						query_id,
						the_filter_spec,
						query_id_full_hits
					);
				}
			}

			// Clear all data in hit_list_by_query_id
			// and wipe prev_query_id_and_hits_ref and to_be_erased_query_id
			hit_list_by_query_id.clear();
			prev_query_id_and_hits_ref = boost::none;
			to_be_erased_query_id      = boost::none;

			// Signal the hits_processor to finish all work
			processor_ptr->finish_work();
		}

		/// \brief Getter for the filter spec to define how to filter the hits 
		inline const crh_filter_spec & read_and_process_mgr::get_filter_spec() const {
			return the_filter_spec;
		}

		/// \brief Get the filter query IDs from the crh_filter_spec in the specified read_and_process_mgr
		///
		/// \relates read_and_process_mgr
		inline const str_vec & get_filter_query_ids(const read_and_process_mgr &arg_read_and_process_mgr ///< The read_and_process_mgr to query
		                                            ) {
			return arg_read_and_process_mgr.get_filter_spec().get_filter_query_ids();
		}

		/// \brief Whether the specified query ID should be skipped according to the crh_filter_spec in the specified read_and_process_mgr
		///
		/// \relates read_and_process_mgr
		inline bool should_skip_query_id(const read_and_process_mgr &arg_read_and_process_mgr, ///< The read_and_process_mgr to query
		                                 const std::string          &arg_query_id              ///< The query ID to test
		                                 ) {
			return should_skip_query_id( arg_read_and_process_mgr.get_filter_spec(), arg_query_id );
		}

		/// \brief Whether the specified query ID should be skipped according to the crh_filter_spec in the specified read_and_process_mgr
		///
		/// \relates read_and_process_mgr
		inline bool should_skip_query_id(const read_and_process_mgr &arg_read_and_process_mgr, ///< The read_and_process_mgr to query
		                                 const boost::string_ref    &arg_query_id              ///< The query ID to test
		                                 ) {
			return should_skip_query_id( arg_read_and_process_mgr.get_filter_spec(), arg_query_id );
		}

	} // namespace rslv
} // namespace cath

#endif
