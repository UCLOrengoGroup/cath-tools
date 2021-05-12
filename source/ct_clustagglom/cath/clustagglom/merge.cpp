/// \file
/// \brief The merge class definitions

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

#include "merge.hpp"

#include <filesystem>
#include <fstream>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/format.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "cath/common/boost_addenda/string_algorithm/split_build.hpp"
#include "cath/common/exception/runtime_error_exception.hpp"
#include "cath/common/file/open_fstream.hpp"
#include "cath/common/type_aliases.hpp"

using namespace ::cath;
using namespace ::cath::common;
using namespace ::cath::clust;

using ::boost::adaptors::transformed;
using ::boost::algorithm::join;
using ::boost::format;
using ::boost::is_space;
using ::boost::token_compress_on;
using ::boost::trim_copy;
using ::std::filesystem::path;
using ::std::ifstream;
using ::std::istream;
using ::std::ofstream;
using ::std::ostream;
using ::std::string;
using ::std::to_string;

/// \brief Generate a string describing the specified merge
///
/// \relates merge
string cath::clust::to_string(const merge &prm_merge ///< The merge to describe
                              ) {
	using ::std::to_string;
	return
		  "merge["
		+ ( format( R"(%5d)" ) % prm_merge.node_a     ).str()
		+ " + "
		+ ( format( R"(%5d)" ) % prm_merge.node_b     ).str()
		+ " -> "
		+ ( format( R"(%5d)" ) % prm_merge.merge_node ).str()
		+ " ("
		+ to_string( prm_merge.dissim )
		+ ")]";
}

/// \brief Generate a string describing the specified merge_vec
///
/// \relates merge
string cath::clust::to_string(const merge_vec &prm_merges ///< The merge_vec to describe
                              ) {
	return join(
		prm_merges
			| transformed( [] (const merge &x) { return to_string( x ); } ),
		"\n"
	);
}

/// \brief Insert a description of the specified merge_list into the specified ostream
void cath::clust::write_merge_list(ostream         &prm_os,    ///< The ostream into which the description should be inserted
                                   const merge_vec &prm_merges ///< The merge_list to describe
                                   ) {
	using ::std::to_string;
	prm_os << join(
		prm_merges
			| transformed( [] (const merge &x) {
				return
					  ( format( R"(%5d)" ) % x.node_a     ).str()
					+ "\t"
					+ ( format( R"(%5d)" ) % x.node_b     ).str()
					+ "\t"
					+ ( format( R"(%5d)" ) % x.merge_node ).str()
					+ "\t"
					+ to_string( x.dissim )
					+ "\n";
			} ),
		""
	);
}

/// \brief Write a description of the specified merge_list into the specified file
void cath::clust::write_merge_list(const path      &prm_out_file, ///< The file to which the description should be written
                                   const merge_vec &prm_merges    ///< The merge_list to describe
                                   ) {
	ofstream out_stream = open_ofstream( prm_out_file );
	write_merge_list( out_stream, prm_merges );
	out_stream.close();
}

/// \brief Read a merge_list from the specified istream
///
/// \relates merge
merge_vec cath::clust::read_merge_list(istream &prm_merges_istream ///< The istream to read the merge_list from
                                       ) {
	string line_string;
	merge_vec results;
	while ( getline( prm_merges_istream, line_string ) ) {
		const str_vec line_parts = split_build<str_vec>( trim_copy( line_string ), is_space(), token_compress_on );

		if ( line_parts.size() != 4 ) {
			BOOST_THROW_EXCEPTION(runtime_error_exception(
				"Whilst reading from merge_list file, unable to find exactly four parts in line:\n"
				+ line_string
			));
		}

		results.emplace_back(
			stoul( line_parts[ 0 ] ),
			stoul( line_parts[ 1 ] ),
			stoul( line_parts[ 2 ] ),
			stod ( line_parts[ 3 ] )
		);
	}

	return results;
}

/// \brief Read a merge_list from the specified file
///
/// \relates merge
merge_vec cath::clust::read_merge_list(const path &prm_merges_file ///< The file to read the merge_list from
                                       ) {
	ifstream merges_istream = open_ifstream( prm_merges_file );
	const auto merges = read_merge_list( merges_istream );
	merges_istream.close();
	return merges;
}


/// \brief Insert a description of the specified merge into the specified ostream
///
/// \relates merge
ostream & cath::clust::operator<<(ostream     &prm_os,   ///< The ostream into which the description should be inserted
                                  const merge &prm_merge ///< The merge to describe
                                  ) {
	prm_os << to_string( prm_merge );
	return prm_os;
}
