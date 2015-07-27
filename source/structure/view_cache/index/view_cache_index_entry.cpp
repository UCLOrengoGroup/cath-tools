/// \file
/// \brief The view_cache_index_entry class definitions

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

#include <boost/numeric/conversion/cast.hpp>

// #include <boost/geometry.hpp>
//#include <boost/geometry/geometries/point.hpp> // ***** TEMPORARY? *****
//#include <boost/geometry/geometries/geometries.hpp> // ***** TEMPORARY? *****

#include "structure/view_cache/index/view_cache_index_entry.h"

using namespace cath::index;
using namespace cath::index::detail;
using namespace cath::geom;
using namespace std;

using boost::numeric_cast;

//template <int T> class bytes_occupied;
//
//bytes_occupied< sizeof( index_type             ) > from_index_size;             // Currently  4,          4, even  2?
//bytes_occupied< sizeof( index_type             ) > to_index_size;               // Currently  4,          4, even  2?
//bytes_occupied< sizeof( view_type              ) > view_size;                   // Currently 24, move to 12
//bytes_occupied< sizeof( frame_quat_rot         ) > frame_size;                  // Currently 16,         16
//bytes_occupied< sizeof( angle_type             ) > from_phi_angle_size;         // Currently  8, move to  4
//bytes_occupied< sizeof( angle_type             ) > from_psi_angle_size;         // Currently  8, move to  4
//bytes_occupied< sizeof( angle_type             ) > to_phi_angle_size;           // Currently  8, move to  4
//bytes_occupied< sizeof( angle_type             ) > to_psi_angle_size;           // Currently  8, move to  4
//
//bytes_occupied< sizeof( view_cache_index_entry ) > view_cache_index_entry_size; // Currently 64, move to 52, even 48?
//
//bytes_occupied< sizeof( uint16_t               ) > uint16_t_size;               //


/// \brief Ctor that completely populates a view_cache_index_entry
view_cache_index_entry::view_cache_index_entry(const index_type  &arg_from_index, ///< The index of the from_residue in its protein
                                               const index_type  &arg_to_index,   ///< The index of the to_residue in its protein
                                               const view_type   &arg_view,       ///< The view of the to_residue from the from_residue
                                               const rotation    &arg_frame,      ///< The coordinate frame of the from_residue
                                               const angle_type  &arg_from_phi,   ///< The phi angle of the from_residue
                                               const angle_type  &arg_from_psi,   ///< The psi angle of the from_residue
                                               const angle_type  &arg_to_phi,     ///< The phi angle of the to_residue
                                               const angle_type  &arg_to_psi      ///< The psi angle of the to_residue
                                               ) : view           ( arg_view       ),
                                                   frame          (
                                                   	make_quat_rot_from_rotation<frame_quat_rot_type>(
                                                   		arg_frame
                                                   	)
                                                   ),
                                                   from_phi_angle ( shift_copy( arg_from_phi ) ),
                                                   from_psi_angle ( shift_copy( arg_from_psi ) ),
                                                   to_phi_angle   ( shift_copy( arg_to_phi   ) ),
                                                   to_psi_angle   ( shift_copy( arg_to_psi   ) ),
                                                   from_index     ( arg_from_index ),
                                                   to_index       ( arg_to_index   ) {
}

/// \brief TODOCUMENT
///
/// \relates view_cache_index_entry
ostream & cath::index::operator<<(ostream                      &arg_os,   ///< TODOCUMENT
                                  const view_cache_index_entry &arg_cache ///< TODOCUMENT
                                  ) {
  arg_os << "view_cache_index_entry[from-to-indices:";
  arg_os << arg_cache.get_from_index();
  arg_os << ",";
  arg_os << arg_cache.get_to_index();
  arg_os << ", view:[";
  arg_os << arg_cache.get_view().get<0>();
  arg_os << ", ";
  arg_os << arg_cache.get_view().get<1>();
  arg_os << ", ";
  arg_os << arg_cache.get_view().get<2>();
  arg_os << "], frame:";
  arg_os << arg_cache.get_frame();
  arg_os << ", from_phi:";
  arg_os << arg_cache.get_from_phi_angle();
  arg_os << ", from_psi:";
  arg_os << arg_cache.get_from_psi_angle();
  arg_os << ", to_phi:";
  arg_os << arg_cache.get_to_phi_angle();
  arg_os << ", to_psi:";
  arg_os << arg_cache.get_to_psi_angle();
  arg_os << "]";
  return arg_os;
}

/// \brief Factory function to build a view_cache_index_entry from a protein and the from/to indices
///
/// \relates view_cache_index_entry
view_cache_index_entry cath::index::detail::make_view_cache_index_entry(const protein &arg_protein,    ///< The protein containing the residues the view_cache_index_entry should represent
                                                                        const size_t  &arg_from_index, ///< The index of the from_residue the view_cache_index_entry should represent
                                                                        const size_t  &arg_to_index    ///< The index of the to_residue the view_cache_index_entry should represent
                                                                        ) {
	return view_cache_index_entry(
		numeric_cast<index_type>( arg_from_index ),
		numeric_cast<index_type>( arg_to_index   ),
		view_vector( arg_protein, arg_from_index, arg_to_index ),
		view_frame ( arg_protein, arg_from_index, arg_to_index ),
		convert_angle_type<angle_base_type>( arg_protein.get_residue_ref_of_index( arg_from_index ).get_phi_angle() ),
		convert_angle_type<angle_base_type>( arg_protein.get_residue_ref_of_index( arg_from_index ).get_psi_angle() ),
		convert_angle_type<angle_base_type>( arg_protein.get_residue_ref_of_index( arg_to_index   ).get_phi_angle() ),
		convert_angle_type<angle_base_type>( arg_protein.get_residue_ref_of_index( arg_to_index   ).get_psi_angle() )
	);
}
