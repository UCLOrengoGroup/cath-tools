/// \file
/// \brief The clustmap_output_spec class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_OPTIONS_SPEC_CLUSTMAP_OUTPUT_SPEC_HPP
#define _CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_OPTIONS_SPEC_CLUSTMAP_OUTPUT_SPEC_HPP

#include "cath/common/path_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"

namespace cath {
	namespace clust {

		/// \brief Specify the outputs for cath-map-clusters
		class clustmap_output_spec final {
		private:
			/// \brief A batch ID to optionally append as a last column to the output
			str_opt  append_batch_id;

			/// \brief An optional file to which output should be redirected
			path_opt output_to_file;

			/// \brief An optional file to which a Markdown summary should be written
			path_opt summarise_to_file;

			/// \brief Whether to print out per-domain mapping info
			bool print_domain_mapping = DEFAULT_PRINT_DOMAIN_MAPPING;

		public:
			/// \brief Default value for whether to print out per-domain mapping info
			static constexpr bool DEFAULT_PRINT_DOMAIN_MAPPING = false;

			/// \brief Default ctor
			clustmap_output_spec() = default;

			[[nodiscard]] const str_opt & get_append_batch_id() const;
			[[nodiscard]] const path_opt &get_output_to_file() const;
			[[nodiscard]] const path_opt &get_summarise_to_file() const;
			[[nodiscard]] const bool &    get_print_domain_mapping() const;

			clustmap_output_spec & set_append_batch_id(const str_opt &);
			clustmap_output_spec & set_output_to_file(const path_opt &);
			clustmap_output_spec & set_summarise_to_file(const path_opt &);
			clustmap_output_spec & set_print_domain_mapping(const bool &);
		};

		/// \brief Getter for a batch ID to optionally append as a last column to the output
		inline const str_opt & clustmap_output_spec::get_append_batch_id() const {
			return append_batch_id;
		}

		/// \brief Getter for an optional file to which output should be redirected
		inline const path_opt & clustmap_output_spec::get_output_to_file() const {
			return output_to_file;
		}

		/// \brief Getter for an optional file to which a Markdown summary should be written
		inline const path_opt & clustmap_output_spec::get_summarise_to_file() const {
			return summarise_to_file;
		}


		/// \brief Getter for whether to print out per-domain mapping info
		inline const bool & clustmap_output_spec::get_print_domain_mapping() const {
			return print_domain_mapping;
		}

		/// \brief Setter for a batch ID to optionally append as a last column to the output
		inline clustmap_output_spec & clustmap_output_spec::set_append_batch_id(const str_opt &prm_append_batch_id ///< A batch ID to optionally append as a last column to the output
		                                                                        ) {
			append_batch_id = prm_append_batch_id;
			return *this;
		}

		/// \brief Setter for an optional file to which output should be redirected
		inline clustmap_output_spec & clustmap_output_spec::set_output_to_file(const path_opt &prm_output_to_file ///< An optional file to which output should be redirected
		                                                                       ) {
			output_to_file = prm_output_to_file;
			return *this;
		}

		/// \brief Setter for an optional file to which a Markdown summary should be written
		inline clustmap_output_spec & clustmap_output_spec::set_summarise_to_file(const path_opt &prm_summarise_to_file ///< An optional file to which a Markdown summary should be written
		                                                                          ) {
			summarise_to_file = prm_summarise_to_file;
			return *this;
		}

		/// \brief Setter for whether to print out per-domain mapping info
		inline clustmap_output_spec & clustmap_output_spec::set_print_domain_mapping(const bool &prm_print_domain_mapping ///< Whether to print out per-domain mapping info
		                                                                             ) {
			print_domain_mapping = prm_print_domain_mapping;
			return *this;
		}

	} // namespace clust
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_CLUSTER_CATH_CLUSTER_OPTIONS_SPEC_CLUSTMAP_OUTPUT_SPEC_HPP
