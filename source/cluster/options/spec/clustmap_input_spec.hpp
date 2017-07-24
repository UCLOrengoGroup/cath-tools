/// \file
/// \brief The clustmap_input_spec class header

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

#ifndef _CATH_TOOLS_SOURCE_CLUSTER_OPTIONS_SPEC_CLUSTMAP_INPUT_SPEC_H
#define _CATH_TOOLS_SOURCE_CLUSTER_OPTIONS_SPEC_CLUSTMAP_INPUT_SPEC_H

#include <boost/optional.hpp>

#include "common/path_type_aliases.hpp"

namespace cath {
	namespace clust {

		/// \brief Specify the inputs for cath-map-clusters
		class clustmap_input_spec final {
		private:
			/// \brief The cluster-membership file for the working clusters
			boost::filesystem::path working_clustmemb_file;

			/// \brief An optional file specify a cluster-membership file for map-from clusters
			path_opt                map_from_clustmemb_file;

			/// \brief Whether to read batches from working_clustmemb_file (rather than cluster membership directly)
			bool                    read_batches_from_input = DEFAULT_READ_BATCHES_FROM_INPUT;

		public:
			/// \brief Default value for whether to read batches from working_clustmemb_file (rather than cluster membership directly)
			static constexpr bool DEFAULT_READ_BATCHES_FROM_INPUT = false;

			/// \brief Default ctor
			clustmap_input_spec() = default;

			const boost::filesystem::path & get_working_clustmemb_file() const;
			const path_opt & get_map_from_clustmemb_file() const;
			const bool & get_read_batches_from_input() const;

			clustmap_input_spec & set_working_clustmemb_file(const boost::filesystem::path &);
			clustmap_input_spec & set_map_from_clustmemb_file(const path_opt &);
			clustmap_input_spec & set_read_batches_from_input(const bool &);
		};

		/// \brief Getter for the cluster-membership file for the working clusters
		inline const boost::filesystem::path & clustmap_input_spec::get_working_clustmemb_file() const {
			return working_clustmemb_file;
		}

		/// \brief Getter for an optional file specify a cluster-membership file for map-from clusters
		inline const path_opt & clustmap_input_spec::get_map_from_clustmemb_file() const {
			return map_from_clustmemb_file;
		}

		/// \brief Getter for whether to read batches from working_clustmemb_file (rather than cluster membership directly)
		inline const bool & clustmap_input_spec::get_read_batches_from_input() const {
			return read_batches_from_input;
		}

		/// \brief Setter for the cluster-membership file for the working clusters
		inline clustmap_input_spec & clustmap_input_spec::set_working_clustmemb_file(const boost::filesystem::path &arg_working_clustmemb_file ///< The cluster-membership file for the working clusters
		                                                                             ) {
			working_clustmemb_file = arg_working_clustmemb_file;
			return *this;
		}

		/// \brief Setter for an optional file specify a cluster-membership file for map-from clusters
		inline clustmap_input_spec & clustmap_input_spec::set_map_from_clustmemb_file(const path_opt &arg_map_from_clustmemb_file ///< An optional file specify a cluster-membership file for map-from clusters
		                                                                              ) {
			map_from_clustmemb_file = arg_map_from_clustmemb_file;
			return *this;
		}

		/// \brief Setter for whether to read batches from working_clustmemb_file (rather than cluster membership directly)
		inline clustmap_input_spec & clustmap_input_spec::set_read_batches_from_input(const bool &arg_read_batches_from_input ///< Whether to read batches from working_clustmemb_file (rather than cluster membership directly)
		                                                                              ) {
			read_batches_from_input = arg_read_batches_from_input;
			return *this;
		}

	} // namespace clust
} // namespace cath

#endif
