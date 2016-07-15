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

#include "common/type_aliases.h"
#include "resolve_hits/hit.h"
#include "resolve_hits/hit_list.h"

#include <future>

namespace cath {
	namespace rslv {

		/// \brief Manage the hits being read from the input file
		///
		/// To make the code efficient, we want to be able to do stuff to the
		/// results that are being read from the input file before the whole input
		/// file has been parsed (particularly because the file may contain input
		/// data for many query_ids and there's no reason they should all have to be
		/// simultaneously held in memory).
		///
		/// This class separates out the reading code from the code that processes
		/// the results as they come in.
		///
		/// This provides separation between the reading and the processing
		class read_and_resolve_mgr final {
		private:
			/// \brief Whether this read_and_resolve_mgr is currently active (and has a query_id)
			bool active = false;

			/// \brief The query_id of the results that are currently coming in
			std::string query_id;

			/// \brief The hits that have been read in so far
			hit_vec the_hits;

			/// \brief The labels that have been read in so far
			str_vec labels;


			/// \brief (A reference_wrapper to) the output stream to which the results should be written
			std::reference_wrapper<std::ostream> output_stream;


			/// \brief The future for the std::async() computations to sort and resolve the hits
			///        for a query once they've all been read in
			std::future<void> resolve_future;

			/// \brief The query_id of any data that's currently be sorted & resolved
			std::string computation_query_id;

			/// \brief Any hits that are currently being sorted & resolved
			hit_vec computation_hits;

			/// \brief The labels of any hits that are currently being sorted & resolved
			str_vec computation_labels;

		public:
			explicit read_and_resolve_mgr(std::ostream &);

			void set_query_id(std::string &&);

			const std::string & get_query_id() const;

			void add_hit(const res_arrow &,
			             const res_arrow &,
			             const hit_seg_vec &,
			             const resscr_t &,
			             std::string &&);

			void complete();

			bool is_active() const;

			void final_wait() const;
		};

		/// \brief Ctor from the ostream to which the results should be written
		inline read_and_resolve_mgr::read_and_resolve_mgr(std::ostream &arg_output_stream ///< The ostream to which the results should be written
		                                                  ) : output_stream( arg_output_stream ) {
		}

		/// \brief Set the query_id for a new set of data
		///
		/// \pre `! is_active()` else an invalid_argument_exception will be thrown
		inline void read_and_resolve_mgr::set_query_id(std::string &&arg_query_id ///< The query_id for the new set of data
		                                               ) {
			if ( is_active() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot set query_id from a read_and_resolve_mgr that is already active because a query_id has already been set"));
			}
			query_id = std::move( arg_query_id );
			active = true;
		}

		/// \brief Getter for the current query_id
		///
		/// \pre `is_active()` else an invalid_argument_exception will be thrown
		inline const std::string & read_and_resolve_mgr::get_query_id() const {
			if ( ! is_active() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot get query_id from a read_and_resolve_mgr that isn't active because no query_id has been set"));
			}
			return query_id;
		}

		/// \brief Add a new hit for the current query_id
		///
		/// \pre `is_active()` else an invalid_argument_exception will be thrown
		inline void read_and_resolve_mgr::add_hit(const res_arrow    &arg_start_arrow, ///< The start boundary of the new hit
		                                          const res_arrow    &arg_stop_arrow,  ///< The stop  boundary of the new hit
		                                          const hit_seg_vec  &arg_fragments,   ///< Any fragments of the new hit
		                                          const resscr_t     &arg_score,       ///< The score associated with the new hit
		                                          std::string       &&arg_label        ///< The label associated with the new hit
		                                          ) {
			if ( ! is_active() ) {
				BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot add hit to read_and_resolve_mgr that is not active because no query_id has been set"));
			}
			the_hits.emplace_back(
				arg_start_arrow,
				arg_stop_arrow,
				arg_fragments,
				arg_score,
				labels.size()
			);
			labels.emplace_back( move( arg_label ) );
		}

		/// \brief Getter for whether this is currently active (accepting data for a previously specified query_id)
		inline bool read_and_resolve_mgr::is_active() const {
			return active;
		}

		/// \brief Call complete on the specified read_and_resolve_mgr if it's currently active
		///
		/// \relates read_and_resolve_mgr
		inline void complete_if_active(read_and_resolve_mgr &arg_read_and_resolve_mgr ///< The read_and_resolve_mgr to complete if it's active
		                               ) {
			if ( arg_read_and_resolve_mgr.is_active() ) {
				arg_read_and_resolve_mgr.complete();
			}
		}

	}
}

#endif
