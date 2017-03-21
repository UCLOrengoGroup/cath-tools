/// \file
/// \brief The hit_boundary_output class definitions

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

#include "hit_boundary_output.hpp"

using namespace cath::rslv;

static_assert(   hit_boundary_output_of_output_trimmed_hits( true                         ) == hit_boundary_output::TRIMMED, "" );
static_assert(   hit_boundary_output_of_output_trimmed_hits( false                        ) == hit_boundary_output::ORIG,    "" );
static_assert(   means_output_trimmed_hits                 ( hit_boundary_output::TRIMMED ),                                 "" );
static_assert( ! means_output_trimmed_hits                 ( hit_boundary_output::ORIG    ),                                 "" );
