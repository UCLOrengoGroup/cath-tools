/// \file
/// \brief The alnd_rgn definitions

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

#include "alnd_rgn.hpp"

#include <boost/range/adaptor/transformed.hpp>
#include <boost/algorithm/string/join.hpp>

#include "cath/seq/seq_seg.hpp"

using namespace ::cath::rslv;
using namespace ::cath::seq;

using ::boost::algorithm::join;
using ::boost::adaptors::transformed;
using ::std::string;

/// \brief Ctor from start on each sequence and length
alnd_rgn::alnd_rgn(seq_arrow       prm_start_res_a, ///< The start of the aligned region in the first  sequence
                   seq_arrow       prm_start_res_b, ///< The start of the aligned region in the second sequence
                   const residx_t &prm_length       ///< The length of the aligned region
                   ) noexcept : start_res_a( std::move( prm_start_res_a ) ),
                                start_res_b( std::move( prm_start_res_b ) ),
                                length     ( prm_length                   ) {
}

/// \brief Getter for the start of the aligned region in the first  sequence
const seq_arrow & alnd_rgn::get_start_res_a() const {
	return start_res_a;
}

/// \brief Getter for the start of the aligned region in the second sequence
const seq_arrow & alnd_rgn::get_start_res_b() const {
	return start_res_b;
}

/// \brief Getter for the length of the aligned region
const residx_t & alnd_rgn::get_length() const {
	return length;
}

/// \brief Generate a string describing the specified alnd_rgn
///
/// \relates alnd_rgn
string cath::rslv::to_string(const alnd_rgn &prm_alnd_rgn ///< The alnd_rgn to describe
	                         ) {
	return to_simple_seg_string(
			prm_alnd_rgn.get_start_res_a(),
			prm_alnd_rgn.get_start_res_a() + prm_alnd_rgn.get_length()
		)
		+ ","
		+ to_simple_seg_string(
			prm_alnd_rgn.get_start_res_b(),
			prm_alnd_rgn.get_start_res_b() + prm_alnd_rgn.get_length()
		);
}

/// \brief Generate a string describing the specified alnd_rgns
///
/// \relates alnd_rgn
string cath::rslv::to_string(const alnd_rgn_vec &prm_alnd_rgns ///< The alnd_rgns to describe
	                         ) {
	return join(
		prm_alnd_rgns
			| transformed( [] (const alnd_rgn &x) { return to_string( x ); } ),
		";"
	);
}

