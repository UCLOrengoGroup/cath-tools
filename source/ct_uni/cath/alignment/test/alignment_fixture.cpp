/// \file
/// \brief The alignment_fixture class definitions

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

#include "alignment_fixture.hpp"

#include "cath/alignment/align_type_aliases.hpp"
#include "cath/common/boost_addenda/range/to_vector.hpp"
#include "cath/common/size_t_literal.hpp"
#include "cath/structure/entry_querier/residue_querier.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::common;
using namespace ::std;

using ::std::nullopt;

alignment alignment_fixture::aln_a_a() {
	return alignment{ { aln_list_a | to_vector, aln_list_a | to_vector } };
}

alignment alignment_fixture::aln_a_b() {
	return alignment{ { aln_list_a | to_vector, aln_list_b | to_vector } };
}

alignment alignment_fixture::aln_b_a() {
	return alignment{ { aln_list_b | to_vector, aln_list_a | to_vector } };
}

alignment alignment_fixture::aln_long_long() {
	return alignment{ { aln_list_long | to_vector, aln_list_long | to_vector } };
}
