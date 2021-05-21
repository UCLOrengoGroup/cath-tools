/// \file
/// \brief The superfamily_of_domain class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_HOMCHECK_TOOLS_SUPERFAMILY_OF_DOMAIN_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_HOMCHECK_TOOLS_SUPERFAMILY_OF_DOMAIN_HPP

#include <filesystem>
#include <regex>
#include <string>
#include <unordered_map>

#include "cath/common/type_aliases.hpp"

namespace cath {
	namespace homcheck {

		namespace detail {
			/// \brief Simple predicate function object to return whether a string is a valid superfamily_id
			class is_valid_superfamily_id final {
			private:
				/// \brief The regular expression used to determine whether a string is a valid CATH superfamily ID
				::std::regex SUPERFAMILY_ID_REGEX{ R"(^\d+\.\d+\.\d+\.\d+$)" };

			public:
				/// \brief Simple predicate function operator to return whether a string is a valid superfamily_id
				inline bool operator()(const std::string &prm_superfamily_id_string ///< The string to check
				                       ) const {
				    return regex_search( prm_superfamily_id_string, SUPERFAMILY_ID_REGEX );
				}
			};

			class is_valid_cath_node_id final {
			private:
				/// \brief The regular expression used to determine whether a string is a valid CATH superfamily ID
				::std::regex NODE_ID_REGEX{ R"(^\d+(\.\d+){0,3}$)" };

			public:
				/// \brief Simple predicate function operator to return whether a string is a valid cath_node_id
				inline bool operator()(const std::string &prm_cath_node_id_string ///< The string to check
									   ) const {
					return regex_search( prm_cath_node_id_string, NODE_ID_REGEX );
				}
			};

			std::string fold_of_superfamily_id(const std::string &);

		} // namespace detail

		/// \brief A lookup from domain_id to the superfamily in which that domain is currently classified
		///        (or will be classified after actions suggested by this code)
		class superfamily_of_domain final {
		private:
			/// \brief Hash from domain_id to superfamily ID
			///
			/// Pre-existing superfamily_ids satisfy meet is_valid_superfamily_id()
			/// New superfamily IDs are built from pre-existing ones and look like: 2.60.40.new_sf_in_fold_of_1cukA01
			std::unordered_map<std::string, std::string> sf_of_dom;

			static bool is_created_sf(const std::string &);

		public:
			superfamily_of_domain() = default;
			explicit superfamily_of_domain(const str_str_pair_vec &);

			size_t size() const;

			bool is_in_new_superfamily(const std::string &) const;
			bool has_superfamily_of_domain(const std::string &) const;
			const std::string & get_superfamily_of_domain(const std::string &) const;
			bool is_in_created_sf(const std::string &) const;

			void add_domain_in_new_sf_in_fold_of_domain(const std::string &,
			                                            const std::string &);

			// See lines ~600-700 in /usr/local/svn/source/cathcgi/trunk/update/HomCheck.pl
//			$oNewHClassification = $oDB->addChildNodeToClassificationNoCommit($oNewTClassification, $strNewHName, $strNewHComment);
//			$aDomainMessages     = $oDB->assignDomainNoCommit($oDomain, $oRelatedDomain,  $domainHistoryAssignmentEventId,  'DOMAIN_ASSIGNED', $strComment, $oNewHClassification, undef, $xmlDataStructure, my $preserveSolid = undef);
//			void add_(const std::string &);
		};

		superfamily_of_domain parse_superfamily_of_domain(std::istream &);
		superfamily_of_domain parse_superfamily_of_domain(const ::std::filesystem::path &);
		superfamily_of_domain parse_superfamily_of_domain(const std::string &);

	} // namespace homcheck
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SCORE_HOMCHECK_TOOLS_SUPERFAMILY_OF_DOMAIN_HPP
