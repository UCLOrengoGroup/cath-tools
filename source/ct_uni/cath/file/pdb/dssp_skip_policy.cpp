/// \file
/// \brief The dssp_skip_policy class definitions

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

#include "dssp_skip_policy.hpp"

using namespace ::cath::file;

static_assert( angle_skipping_of_dssp_skip_policy( dssp_skip_policy::SKIP__BREAK_ANGLES           ) == dssp_skip_angle_skipping::BREAK_ANGLES      );
static_assert( angle_skipping_of_dssp_skip_policy( dssp_skip_policy::DONT_SKIP__BREAK_ANGLES      ) == dssp_skip_angle_skipping::BREAK_ANGLES      );
static_assert( angle_skipping_of_dssp_skip_policy( dssp_skip_policy::DONT_SKIP__DONT_BREAK_ANGLES ) == dssp_skip_angle_skipping::DONT_BREAK_ANGLES );

static_assert( res_skipping_of_dssp_skip_policy  ( dssp_skip_policy::SKIP__BREAK_ANGLES           ) == dssp_skip_res_skipping::SKIP                );
static_assert( res_skipping_of_dssp_skip_policy  ( dssp_skip_policy::DONT_SKIP__BREAK_ANGLES      ) == dssp_skip_res_skipping::DONT_SKIP           );
static_assert( res_skipping_of_dssp_skip_policy  ( dssp_skip_policy::DONT_SKIP__DONT_BREAK_ANGLES ) == dssp_skip_res_skipping::DONT_SKIP           );
