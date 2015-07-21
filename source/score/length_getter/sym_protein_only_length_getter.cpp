/// \file
/// \brief The sym_protein_only_length_getter class definitions

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#include "sym_protein_only_length_getter.h"

#include <boost/assign/ptr_list_inserter.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include "score/length_getter/geometric_mean_length_getter.h"
#include "score/length_getter/length_of_longer_getter.h"
#include "score/length_getter/length_of_shorter_getter.h"
#include "score/length_getter/mean_length_getter.h"

using namespace boost::logic;
using namespace cath::align;
using namespace cath::score;
using namespace std;

using boost::assign::ptr_push_back;
using boost::ptr_vector;

/// \brief TODOCUMENT
unique_ptr<protein_only_length_getter> sym_protein_only_length_getter::do_protein_only_clone() const {
	return { sym_protein_only_clone() };
}

/// \brief TODOCUMENT
unique_ptr<sym_protein_only_length_getter> sym_protein_only_length_getter::sym_protein_only_clone() const {
	return do_sym_protein_only_clone();
}

/// \brief TODOCUMENT
ptr_vector<sym_protein_only_length_getter> cath::score::get_all_sym_protein_only_length_getters() {
	ptr_vector<sym_protein_only_length_getter> all_sym_protein_only_length_getters;
	ptr_push_back< length_of_longer_getter      >( all_sym_protein_only_length_getters )( );
	ptr_push_back< length_of_shorter_getter     >( all_sym_protein_only_length_getters )( );
	ptr_push_back< mean_length_getter           >( all_sym_protein_only_length_getters )( );
	ptr_push_back< geometric_mean_length_getter >( all_sym_protein_only_length_getters )( );
	return all_sym_protein_only_length_getters;
}

