/// \file
/// \brief The dyn_prog_aligner class definitions

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

#include "dyn_prog_aligner.hpp"

#include "cath/alignment/alignment.hpp"
#include "cath/common/clone/check_uptr_clone_against_this.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::align::gap;
using namespace ::cath::common;
using namespace ::std;

/// \brief Standard approach to achieving a virtual copy-ctor
unique_ptr<dyn_prog_aligner> dyn_prog_aligner::clone() const {
	return check_uptr_clone_against_this( do_clone(), *this );
}

/// \brief NVI pass-through to virtual do_align() method
score_alignment_pair dyn_prog_aligner::align(const dyn_prog_score_source &prm_scorer,      ///< TODOCUMENT
                                             const gap_penalty           &prm_gap_penalty, ///< The gap penalty to be applied for each gap step (ie for opening OR extending a gap)
                                             const size_type             &prm_window_width ///< TODOCUMENT
                                             ) const {
	return do_align(prm_scorer, prm_gap_penalty, prm_window_width);
}

