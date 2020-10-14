/// \file
/// \brief The dyn_prog_score_source class definitions

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

#include "dyn_prog_score_source.hpp"

#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/common/temp_check_offset_1.hpp"

using namespace ::cath;
using namespace ::cath::align;

/// \brief An NVI pass-through method to get the number of elements in the first sequence
size_t dyn_prog_score_source::get_length_a() const {
	return do_get_length_a();
}

/// \brief An NVI pass-through method to get the number of elements in the second sequence
size_t dyn_prog_score_source::get_length_b() const {
	return do_get_length_b();
}
