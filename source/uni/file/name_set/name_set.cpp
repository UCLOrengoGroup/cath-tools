/// \file
/// \brief The name_set class definitions

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

#include "name_set.hpp"

#include "common/hash/hash_value_combine.hpp"

#include <iostream>

using namespace cath::common;

using std::string;
using std::ostream;

/// \brief If the specified name_set has a domain_name_from_regions, return that,
///        else if it has a specified_id, return that, else return the name_from_acq
///
/// \relates name_set
string cath::file::get_domain_or_specified_or_name_from_acq(const name_set &prm_name_set ///< The name_set to query
                                                            ) {
	return
		prm_name_set.get_domain_name_from_regions() ? *prm_name_set.get_domain_name_from_regions() :
		prm_name_set.get_specified_id()             ? *prm_name_set.get_specified_id()             :
		                                               prm_name_set.get_name_from_acq();
}

/// \brief Generate a string describing the specified name_set
///
/// \relates name_set
string cath::file::to_string(const name_set &prm_name_set ///< The name_set to describe
                             ) {
	return
		  R"(name_set[name_from_acq:")"
		+ prm_name_set.get_name_from_acq()
		+ R"(")"
		+ (
			prm_name_set.get_primary_source_file()
			? ( R"(, prim_src_file:")" + prm_name_set.get_primary_source_file()->string() + R"(")" )
			: ""
		)
		+ (
			prm_name_set.get_specified_id()
			? ( R"(, spec_id:")" + *prm_name_set.get_specified_id() + R"(")" )
			: ""
		)
		+ (
			prm_name_set.get_domain_name_from_regions()
			? ( R"(, dom_id:")" + *prm_name_set.get_domain_name_from_regions() + R"(")" )
			: ""
		)
		+ "]";
}

/// \brief Insert a description of the specified name_set into the specified ostream
///
/// \relates name_set
ostream & cath::file::operator<<(ostream        &prm_os,      ///< The ostream into which the description should be inserted
                                 const name_set &prm_name_set ///< The name_set to describe
                                 ) {
	prm_os << to_string( prm_name_set );
	return prm_os;
}


/// \brief Hash the details of the specified name_set into the specified seed value
///
/// \relates name_set
void cath::file::non_crypto_hash(size_t         &prm_init_hash_value, ///< The initial hash seed
                                 const name_set &prm_name_set         ///< The name_set to include in the hash
                                 ) {
	// hash_value_combine( prm_init_hash_value, std::hash<size_t>()( index ) );

	// hash_value_combine( prm_init_hash_value, std::hash<std::string>()( prm_name_set.get_name_from_acq() ) );
	// hash_value_combine( prm_init_hash_value, std::hash<std::string>()( prm_name_set.get_primary_source_file()      ) );
	// hash_value_combine( prm_init_hash_value, std::hash<std::string>()( prm_name_set.get_specified_id()             ) );
	// hash_value_combine( prm_init_hash_value, std::hash<std::string>()( prm_name_set.get_domain_name_from_regions() ) );

	const string   &name_from_acq                = prm_name_set.get_name_from_acq();
	const path_opt &primary_source_file_opt      = prm_name_set.get_primary_source_file();
	const str_opt  &specified_id_opt             = prm_name_set.get_specified_id();
	const str_opt  &domain_name_from_regions_opt = prm_name_set.get_domain_name_from_regions();

	// Follow libstdc++'s lead of attempting to make the value used for none/nullopt an "unusual" value
	static constexpr size_t NULL_HASH_VAL = static_cast<size_t>( -3333 );

	hash_value_combine( prm_init_hash_value,                                std::hash<std::string>()( name_from_acq                     )                 );
	hash_value_combine( prm_init_hash_value, primary_source_file_opt      ? std::hash<std::string>()( primary_source_file_opt->string() ) : NULL_HASH_VAL );
	hash_value_combine( prm_init_hash_value, specified_id_opt             ? std::hash<std::string>()( *specified_id_opt                 ) : NULL_HASH_VAL );
	hash_value_combine( prm_init_hash_value, domain_name_from_regions_opt ? std::hash<std::string>()( *domain_name_from_regions_opt     ) : NULL_HASH_VAL );

	// prm_init_hash_value;
	// prm_name_set;
}

/// \brief Hash the details of the specified name_set into the specified seed value
///
/// \relates name_set
size_t cath::file::non_crypto_hash_copy(size_t          prm_init_hash_value, ///< The initial hash seed
                                        const name_set &prm_name_set         ///< The name_set to include in the hash
                                        ) {
	non_crypto_hash( prm_init_hash_value, prm_name_set );
	return prm_init_hash_value;
}

