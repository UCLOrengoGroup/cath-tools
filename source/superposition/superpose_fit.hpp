/// \file
/// \brief The superpose_fit class header

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

#ifndef _CATH_TOOLS_SOURCE_SUPERPOSITION_SUPERPOSE_FIT_H
#define _CATH_TOOLS_SOURCE_SUPERPOSITION_SUPERPOSE_FIT_H

namespace cath { namespace geom { class coord_list; } }
namespace cath { namespace geom { class rotation; } }
namespace cath { namespace geom { namespace detail { class gsl_matrix_wrp; } } }

namespace cath {
	namespace sup {

		namespace detail {
			geom::detail::gsl_matrix_wrp cross_covariance_matrix(const geom::coord_list &,
			                                                     const geom::coord_list &);
		} // namespace detail

		geom::rotation superpose_fit_1st_to_2nd(const geom::coord_list &,
		                                        const geom::coord_list &);

		geom::rotation superpose_fit_2nd_to_1st(const geom::coord_list &,
		                                        const geom::coord_list &);

	} // namespace sup
} // namespace cath
#endif
