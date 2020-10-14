/// \file
/// \brief The superposition_outputter_list class definitions

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

#include "superposition_outputter_list.hpp"

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "cath/common/boost_addenda/ptr_container/unique_ptr_functions.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/outputter/superposition_outputter/superposition_outputter.hpp"

using namespace ::cath::opts;
using namespace ::cath::sup;

using ::boost::adaptors::transformed;
using ::boost::algorithm::any_of;
using ::boost::algorithm::join;
using ::boost::string_ref;
using ::std::ostream;
using ::std::string;

/// \brief TODOCUMENT
void superposition_outputter_list::push_back(const superposition_outputter &prm_outputter ///< TODOCUMENT
                                             ) {
	cath::common::push_back( outputters, prm_outputter.clone() );
}

/// \brief Return the number of outputters currently in the list
size_t superposition_outputter_list::size() const {
	return outputters.size();
}

/// \brief TODOCUMENT
bool superposition_outputter_list::empty() const {
	return outputters.empty();
}

/// \brief TODOCUMENT
superposition_outputter_list::const_iterator superposition_outputter_list::begin() const {
	return common::cbegin( outputters );
}

/// \brief TODOCUMENT
superposition_outputter_list::const_iterator superposition_outputter_list::end() const {
	return common::cend( outputters );
}

/// \brief Generate a string describing the specified superposition_outputter_list
///
/// \relates superposition_outputter_list
string cath::opts::to_string(const superposition_outputter_list &prm_superposition_outputters ///< Generate a string describing the specified superposition_outputter_list
                             ) {
	return "superposition_outputter_list[ "
		+ join(
			prm_superposition_outputters
				| transformed( [] (const superposition_outputter &x) {
					return x.get_name();
				} ),
			", "
		)
		+ " ]";
}

/// \brief Insert a description of the specified superposition_outputter_list into the specified ostream
///
/// \relates superposition_outputter_list
ostream & cath::opts::operator<<(ostream                            &prm_os,                      ///< The ostream into which the description should be inserted
                                 const superposition_outputter_list &prm_superposition_outputters ///< Generate a string describing the specified superposition_outputter_list
                                 ) {
	prm_os << to_string( prm_superposition_outputters );
	return prm_os;
}

/// \brief TODOCUMENT
///
/// \relates superposition_outputter_list
void cath::opts::use_all_superposition_outputters(const superposition_outputter_list &prm_superposition_outputters, ///< TODOCUMENT
                                                  const superposition_context        &prm_superposition_context,    ///< TODOCUMENT
                                                  ostream                            &prm_stdout,                   ///< TODOCUMENT
                                                  ostream                            &/*prm_stderr*/                ///< TODOCUMENT
                                                  ) {
	// For each of the superposition_outputters specified by the cath_superpose_options, output the superposition
	for (const superposition_outputter &outputter : prm_superposition_outputters) {
		outputter.output_superposition( prm_superposition_context, prm_stdout );
	}
}

/// \brief TODOCUMENT
///
/// \relates superposition_outputter_list
bool cath::opts::any_superposition_outputters_involve_display_spec(const superposition_outputter_list &prm_superposition_outputters ///< TODOCUMENT
                                                                   ) {
	return any_of( prm_superposition_outputters, [] (const superposition_outputter &x) { return x.involves_display_spec(); } );
}
