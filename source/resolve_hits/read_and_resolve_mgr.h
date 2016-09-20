/// \file
/// \brief The read_and_resolve_mgr class header

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

#ifndef READ_AND_RESOLVE_MGR_H_INCLUDED
#define READ_AND_RESOLVE_MGR_H_INCLUDED

#include <boost/optional.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/utility/string_ref.hpp>

#include "common/algorithm/sort_uniq_build.h"
#include "common/type_aliases.h"
#include "exception/runtime_error_exception.h"
#include "resolve_hits/hit.h"
#include "resolve_hits/hit_list.h"

#include <future>
#include <unordered_map>

namespace cath {
	namespace rslv {

		/// \brief Manage the hits being read from the input file
		///
		/// This stores the hits, indexed by query_id in a unordered_map.
		///
		/// Client code should call add_hit() for each new hit and then process_all_outstanding()
		/// at the end of the input
		///
		/// If `input_hits_are_presorted` then this can trigger asynchronous processing
		/// of a block of hits for a given query_id when it reaches the end. To allow
		/// this, it stores the previous hit's query ID and when a hit with a different
		/// query ID is encountered, it triggers processing of the hits associated with
		/// the previous query ID.
		///
		/// If `! input_hits_are_presorted`, then the results must all be processed synchronously
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
		class read_and_resolve_mgr final {
		private:
			/// \brief A type alias for an unordered_map from the query_id string to the corresponding hits data
			using str_hit_list_umap = std::unordered_map<std::string, std::pair<hit_vec, str_vec> >;

			/// \brief An unordered_map from the query_id string to the corresponding hits data
			str_hit_list_umap hit_list_by_query_id;

			/// \brief Whether or not (the user guaranteees that) the input hits data is presorted
			///        to ensure all hits for a query_id are consecutive
			///
			/// When this is true, it's possible to process and erase a block of hits for a query_id when the
			/// end of the block is detected, which can be faster and more memory efficient.
			bool input_hits_are_presorted = DEFAULT_INPUT_HITS_ARE_PRESORTED;

			/// \brief A string that can be reused for holding a local copy of the query ID
			///        for hashing to avoid having to reallocate on every call to add_hit()
			std::string temp_hashable_query_id;

			/// \brief A type alias for a pair of string and hit data (non-const) reference
			using str_hits_ref_pair = std::pair<std::string, std::pair<hit_vec, str_vec> &>;

			/// \brief The query_id of the previous result and a reference to the corresponding data (if any)
			///
			/// This is used when `input_hits_are_presorted` to indicate that a hit_list under this query_id
			/// can be processed as soon as a result is found with a different query id
			///
			/// The hit data reference avoids many of the repeated calls to the hashing function when input_hits_are_presorted
			boost::optional<str_hits_ref_pair> prev_query_id_and_hits_ref;

			/// \brief A record of the query_id associated with the latest asynchronous processing job
			///        so that the thread can detect if it finds itself trying to add another hit to
			///        data that it assumed was finished and which it has passed off to an async worker
			///        thread.
			///
			/// This is only used by the main thread, never by the async worker thread.
			opt_str to_be_erased_query_id;

			/// \brief (A reference_wrapper to) the output stream to which the results should be written
			std::reference_wrapper<std::ostream> output_stream;


			/// \brief The std::future for the std::async() processing.
			///
			/// At present, this isn't used to pass back results, just to allow the main thread to wait
			/// until an async processing job is complete.
			std::future<void> resolve_future;

			static void process_query_id(const std::string &arg_query_id,
			                             hit_vec &,
			                             str_vec &,
			                             std::ostream &);

			void trigger_async_process_query_id(const std::string &);

			void wait_for_any_active_work() const;

		public:

			/// \brief The default value for whether input hits can be assumed to be pre-sorted
			///
			/// This is false - better not to assume this guarantee
			static constexpr bool DEFAULT_INPUT_HITS_ARE_PRESORTED = false;

			explicit read_and_resolve_mgr(std::ostream &,
			                              const bool & = DEFAULT_INPUT_HITS_ARE_PRESORTED);

			void add_hit(const boost::string_ref &,
			             const res_arrow &,
			             const res_arrow &,
			             const hit_seg_vec &,
			             const resscr_t &,
			             std::string &&);

			void process_all_outstanding();
		};

		

		/// \brief Check that any previously triggered async computation is complete
		inline void read_and_resolve_mgr::wait_for_any_active_work() const {
			if ( resolve_future.valid() ) {
				resolve_future.wait();
			}
		}

		/// \brief Ctor from the ostream to which the results should be written
		inline read_and_resolve_mgr::read_and_resolve_mgr(std::ostream &arg_output_stream,           ///< The ostream to which the results should be written
		                                                  const bool   &arg_input_hits_are_presorted ///< Whether or not the input hits are guaranteed to be presorted
		                                                  ) : input_hits_are_presorted{ arg_input_hits_are_presorted },
		                                                      output_stream           { arg_output_stream            } {
		}

		/// \brief Add a new hit for the current query_id
		///
		/// \pre `is_active()` else an invalid_argument_exception will be thrown
		inline void read_and_resolve_mgr::add_hit(const boost::string_ref &arg_query_id,    ///< A string_ref of the query_id
		                                          const res_arrow         &arg_start_arrow, ///< The start boundary of the new hit
		                                          const res_arrow         &arg_stop_arrow,  ///< The stop  boundary of the new hit
		                                          const hit_seg_vec       &arg_fragments,   ///< Any fragments of the new hit
		                                          const resscr_t          &arg_score,       ///< The score associated with the new hit
		                                          std::string            &&arg_label        ///< The label associated with the new hit
		                                          ) {
			// Store the query_id in a local string temp_hashable_query_id, which can then be
			// used for hashing
			temp_hashable_query_id.assign( arg_query_id.data(), arg_query_id.length() );

			// If input_hits_are_presorted and previous hits had a different query_id, trigger an async worker thread to
			// process the hits associated with that previous query_id
			if ( input_hits_are_presorted && prev_query_id_and_hits_ref && prev_query_id_and_hits_ref->first != temp_hashable_query_id ) {
				trigger_async_process_query_id( prev_query_id_and_hits_ref->first );
				prev_query_id_and_hits_ref = boost::none;
			}

			// If this query_id matches to_be_erased_query_id then something's wrong
			// (probably that the input data violates a given guarantee that it'd be)
			if ( to_be_erased_query_id && *to_be_erased_query_id == temp_hashable_query_id ) {

				// ** Should probably catch this exception and exit with error code somwhere **
				BOOST_THROW_EXCEPTION(common::runtime_error_exception(
					"Attempt to add a hit for a query_id "
					+ temp_hashable_query_id
					+ " that previously appeared to be finished - suggests the input data breaks the guarantee to be sorted by query_id"
				));
			}

			// If this is the same query_id as the previous query, then just re-use the cached reference to that query's hits data
			// otherwise, look it up in hit_list_by_query_id
			auto &the_hits = prev_query_id_and_hits_ref ? prev_query_id_and_hits_ref->second
			                                            : hit_list_by_query_id[ temp_hashable_query_id ];

			// Add the new hit to the query's hits data
			the_hits.first.emplace_back(
				arg_start_arrow,
				arg_stop_arrow,
				arg_fragments,
				arg_score,
				the_hits.second.size()
			);
			the_hits.second.emplace_back( move( arg_label ) );

			// If the input hits are presorted then ensure prev_query_id_and_hits_ref is up-to-date
			// with this hit
			if ( input_hits_are_presorted && ! prev_query_id_and_hits_ref ) {
				prev_query_id_and_hits_ref = str_hits_ref_pair{ temp_hashable_query_id, the_hits };
			}
		}

		/// \brief Process all outstanding data
		inline void read_and_resolve_mgr::process_all_outstanding() {
			// If input_hits_are_presorted, there may be asynchronous processing going on, so wait until that's
			// finished before processing this data (else the results may get interleaved)
			if ( input_hits_are_presorted ) {
				wait_for_any_active_work();
			}


			// Get a sorted list of all the query IDs
			const auto sorted_query_ids = common::sort_build<str_vec>(
				hit_list_by_query_id | boost::adaptors::map_keys
			);

			// Loop over the sorted query IDs
			for (const auto &query_id : sorted_query_ids) {
				// Grab the data for this query_id and process it
				auto &query_id_hit_list_pair = hit_list_by_query_id[ query_id ];
				if ( ! query_id_hit_list_pair.first.empty() ) {
					process_query_id(
						query_id,
						query_id_hit_list_pair.first,
						query_id_hit_list_pair.second,
						output_stream.get()
					);
				}
			}

			// Clear all data in hit_list_by_query_id
			// and wipe prev_query_id_and_hits_ref and to_be_erased_query_id
			hit_list_by_query_id.clear();
			prev_query_id_and_hits_ref = boost::none;
			to_be_erased_query_id      = boost::none;
		}

	}
}

#endif
