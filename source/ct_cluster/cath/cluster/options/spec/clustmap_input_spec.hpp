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

#ifndef _CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_OPTIONS_SPEC_CLUSTMAP_INPUT_SPEC_HPP
#define _CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_OPTIONS_SPEC_CLUSTMAP_INPUT_SPEC_HPP

#include <filesystem>

#include "cath/common/path_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"

namespace cath {
	namespace clust {

		/// \brief Specify the inputs for cath-map-clusters
		class clustmap_input_spec final {
		private:
			/// \brief The cluster-membership file for the working clusters
			::std::filesystem::path working_clustmemb_file;

			/// \brief An optional file specify a cluster-membership file for map-from clusters
			path_opt                map_from_clustmemb_file;

			/// \brief Whether to read batches from working_clustmemb_file (rather than cluster membership directly)
			bool                    read_batches_from_input = DEFAULT_READ_BATCHES_FROM_INPUT;

		public:
			/// \brief Default value for whether to read batches from working_clustmemb_file (rather than cluster membership directly)
			static constexpr bool DEFAULT_READ_BATCHES_FROM_INPUT = false;

			/// \brief Default ctor
			clustmap_input_spec() = default;

			[[nodiscard]] const ::std::filesystem::path &get_working_clustmemb_file() const;
			[[nodiscard]] const path_opt &               get_map_from_clustmemb_file() const;
			[[nodiscard]] const bool &                   get_read_batches_from_input() const;

			clustmap_input_spec & set_working_clustmemb_file(const ::std::filesystem::path &);
			clustmap_input_spec & set_map_from_clustmemb_file(const path_opt &);
			clustmap_input_spec & set_read_batches_from_input(const bool &);
		};

		str_opt get_invalid_description(const clustmap_input_spec &);

		/// \brief Getter for the cluster-membership file for the working clusters
		inline const ::std::filesystem::path & clustmap_input_spec::get_working_clustmemb_file() const {
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
		inline clustmap_input_spec & clustmap_input_spec::set_working_clustmemb_file(const ::std::filesystem::path &prm_working_clustmemb_file ///< The cluster-membership file for the working clusters
		                                                                             ) {
			working_clustmemb_file = prm_working_clustmemb_file;
			return *this;
		}

		/// \brief Setter for an optional file specify a cluster-membership file for map-from clusters
		inline clustmap_input_spec & clustmap_input_spec::set_map_from_clustmemb_file(const path_opt &prm_map_from_clustmemb_file ///< An optional file specify a cluster-membership file for map-from clusters
		                                                                              ) {
			map_from_clustmemb_file = prm_map_from_clustmemb_file;
			return *this;
		}

		/// \brief Setter for whether to read batches from working_clustmemb_file (rather than cluster membership directly)
		inline clustmap_input_spec & clustmap_input_spec::set_read_batches_from_input(const bool &prm_read_batches_from_input ///< Whether to read batches from working_clustmemb_file (rather than cluster membership directly)
		                                                                              ) {
			read_batches_from_input = prm_read_batches_from_input;
			return *this;
		}

	} // namespace clust
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_OPTIONS_SPEC_CLUSTMAP_INPUT_SPEC_HPP
