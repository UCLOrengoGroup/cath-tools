/// \file
/// \brief The alignment_outputter_list class definitions

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

#include "alignment_outputter_list.hpp"

#include "cath/common/boost_addenda/ptr_container/unique_ptr_functions.hpp"
#include "cath/common/cpp14/cbegin_cend.hpp"
#include "cath/outputter/alignment_outputter/alignment_outputter.hpp"

using namespace ::cath::align;
using namespace ::cath::common;
using namespace ::cath::opts;
using namespace ::std;

/// \brief TODOCUMENT
void alignment_outputter_list::push_back(const alignment_outputter &prm_outputter ///< TODOCUMENT
                                         ) {
	cath::common::push_back( outputters, prm_outputter.clone() );
}

/// \brief TODOCUMENT
bool alignment_outputter_list::empty() const {
	return outputters.empty();
}

/// \brief TODOCUMENT
alignment_outputter_list::const_iterator alignment_outputter_list::begin() const {
	return common::cbegin( outputters );
}

/// \brief TODOCUMENT
alignment_outputter_list::const_iterator alignment_outputter_list::end() const {
	return common::cend( outputters );
}

/// \brief TODOCUMENT
///
/// \relates alignment_outputter_list
void cath::opts::use_all_alignment_outputters(const alignment_outputter_list &prm_alignment_outputters, ///< TODOCUMENT
                                              const alignment_context        &prm_alignment_context,    ///< TODOCUMENT
                                              ostream                        &prm_stdout,               ///< TODOCUMENT
                                              ostream                        &/*prm_stderr*/            ///< TODOCUMENT
                                              ) {
	// For each of the alignment_outputters specified by the cath_superpose_options, output the alignment
	for (const alignment_outputter &outputter : prm_alignment_outputters) {
		outputter.output_alignment(prm_alignment_context, prm_stdout);
	}
}

/// \brief TODOCUMENT
///
/// \relates alignment_outputter_list
bool cath::opts::any_alignment_outputters_involve_display_spec(const alignment_outputter_list &prm_alignment_outputters ///< TODOCUMENT
                                                               ) {
	// For each of the superposition_outputters specified by the cath_superpose_options, output the superposition
	for (const alignment_outputter &outputter : prm_alignment_outputters) {
		if ( outputter.involves_display_spec() ) {
			return true;
		}
	}
	return false;
}
