/// \file
/// \brief The dyn_prog_score_source_fixture class definitions

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

#include "dyn_prog_score_source_fixture.hpp"

#include "cath/structure/entry_querier/residue_querier.hpp"

using namespace ::cath;
using namespace ::cath::align;
using namespace ::cath::common;
using namespace ::std;

/// \brief An example score matrix for making an example old_matrix_dyn_prog_score_source
static score_vec_of_vec example_old_score_matrix() {
	return score_vec_of_vec{
		{ 0, 0, 0, 0, 0 },
		{ 0, 2, 6, 7, 1 },
		{ 0, 8, 4, 3, 9 },
	};
}

/// \brief An example score matrix for making an example old_matrix_dyn_prog_score_source
static float_score_vec_vec example_new_score_matrix() {
	return {
		{ 2.0, 6.0, 7.0, 1.0 },
		{ 8.0, 4.0, 3.0, 9.0 },
	};
}

/// \brief An example mask matrix for making an example mask_dyn_prog_score_source
static bool_vec_of_vec example_mask_matrix() {
	return bool_vec_of_vec{
		{ false, false, false, false, false },
		{ false,  true, false,  true, false },
		{ false, false,  true, false,  true },
	};
}

///// \brief Make an example entry_querier_dyn_prog_score_source for testing
//const entry_querier_dyn_prog_score_source dyn_prog_score_source_fixture::make_example_entry_querier_dyn_prog_score_source() {
//	return entry_querier_dyn_prog_score_source(
//		EXAMPLE_RESIDUE_QUERIER(),
//		example_protein_a,
//		example_protein_b,
//		0,
//		0
//	);
//}

/// \brief Make an example mask_dyn_prog_score_source for testing
const mask_dyn_prog_score_source & dyn_prog_score_source_fixture::make_example_mask_dyn_prog_score_source() {
	static const bool_vec_of_vec static_example_mask_matrix = example_mask_matrix();
	static const mask_dyn_prog_score_source the_mask_dyn_prog_score_source(
		static_example_mask_matrix,
		make_example_old_matrix_dyn_prog_score_source()
	);
	return the_mask_dyn_prog_score_source;
}

/// \brief Make an example old_matrix_dyn_prog_score_source for testing
const old_matrix_dyn_prog_score_source & dyn_prog_score_source_fixture::make_example_old_matrix_dyn_prog_score_source() {
	static const score_vec_of_vec static_example_old_score_matrix = example_old_score_matrix();
	static const old_matrix_dyn_prog_score_source the_old_matrix_dyn_prog_score_source(
		static_example_old_score_matrix,
		2,
		4,
		4
	);
	return the_old_matrix_dyn_prog_score_source;
}

/// \brief Make an example new_matrix_dyn_prog_score_source for testing
const new_matrix_dyn_prog_score_source & dyn_prog_score_source_fixture::make_example_new_matrix_dyn_prog_score_source() {
	static const new_matrix_dyn_prog_score_source the_new_matrix_dyn_prog_score_source(
		example_new_score_matrix(),
		2,
		4
	);
	return the_new_matrix_dyn_prog_score_source;
}

/// \brief Make an example sequence_string_dyn_prog_score_source for testing
const sequence_string_dyn_prog_score_source & dyn_prog_score_source_fixture::make_example_sequence_string_dyn_prog_score_source() {
	static const sequence_string_dyn_prog_score_source the_sequence_string_dyn_prog_score_source(
		string{ SEQUENCE_STRING_A },
		string{ SEQUENCE_STRING_B }
	);
	return the_sequence_string_dyn_prog_score_source;
}


///// \brief Template specialisation to make an example entry_querier_dyn_prog_score_source
//template <>
//entry_querier_dyn_prog_score_source dyn_prog_score_source_fixture::make_example_dyn_prog_score_source<entry_querier_dyn_prog_score_source>() {
//	return make_example_entry_querier_dyn_prog_score_source();
//}

/// \brief Template specialisation to make an example mask_dyn_prog_score_source
template <>
mask_dyn_prog_score_source dyn_prog_score_source_fixture::make_example_dyn_prog_score_source<mask_dyn_prog_score_source>() {
	return make_example_mask_dyn_prog_score_source();
}

/// \brief Template specialisation to make an example old_matrix_dyn_prog_score_source
template <>
old_matrix_dyn_prog_score_source dyn_prog_score_source_fixture::make_example_dyn_prog_score_source<old_matrix_dyn_prog_score_source>() {
	return make_example_old_matrix_dyn_prog_score_source();
}

/// \brief Template specialisation to make an example new_matrix_dyn_prog_score_source
template <>
new_matrix_dyn_prog_score_source dyn_prog_score_source_fixture::make_example_dyn_prog_score_source<new_matrix_dyn_prog_score_source>() {
	return make_example_new_matrix_dyn_prog_score_source();
}

/// \brief Template specialisation to make an example sequence_string_dyn_prog_score_source
template <>
sequence_string_dyn_prog_score_source dyn_prog_score_source_fixture::make_example_dyn_prog_score_source<sequence_string_dyn_prog_score_source>() {
	return make_example_sequence_string_dyn_prog_score_source();
}
