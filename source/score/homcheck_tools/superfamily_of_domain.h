/// \file
/// \brief The superfamily_of_domain class header

#ifndef SUPERFAMILY_OF_DOMAIN_H_INCLUDED
#define SUPERFAMILY_OF_DOMAIN_H_INCLUDED

#include <boost/filesystem/path.hpp>

#include "common/type_aliases.h"

#include <regex>
#include <string>
#include <unordered_map>


namespace cath {
	namespace homcheck {

		namespace detail {
			/// \brief Simple predicate function object to return whether a string is a valid superfamily_id
			class is_valid_superfamily_id final {
			private:
				static const std::regex SUPERFAMILY_ID_REGEX;
			public:
				/// \brief Simple predicate function operator to return whether a string is a valid superfamily_id
				inline bool operator()(const std::string &arg_superfamily_id_string ///< The string to check
				                       ) const {
				    return regex_search( arg_superfamily_id_string, SUPERFAMILY_ID_REGEX );
				}
			};

			class is_valid_cath_node_id final {
			private:
				static const std::regex NODE_ID_REGEX;
			public:
				/// \brief Simple predicate function operator to return whether a string is a valid cath_node_id
				inline bool operator()(const std::string &arg_cath_node_id_string ///< The string to check
									   ) const {
					return regex_search( arg_cath_node_id_string, NODE_ID_REGEX );
				}
			};
			std::string fold_of_superfamily_id(const std::string &);
		}

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
			static const std::string NEW_SF_CORE_STRING;

			superfamily_of_domain() = default;
			superfamily_of_domain(const str_str_pair_vec &);

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
		superfamily_of_domain parse_superfamily_of_domain(const boost::filesystem::path &);
		superfamily_of_domain parse_superfamily_of_domain(const std::string &);

	}
}

#endif
