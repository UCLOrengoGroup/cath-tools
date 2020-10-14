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

#include "structure/view_cache/index/view_cache_index_entry.hpp"

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
view_cache_index_entry::view_cache_index_entry(const index_type  &prm_from_index, ///< The index of the from_residue in its protein
                                               const index_type  &prm_to_index,   ///< The index of the to_residue in its protein
                                               view_type          prm_view,       ///< The view of the to_residue from the from_residue
                                               const rotation    &prm_frame,      ///< The coordinate frame of the from_residue
                                               const angle_type  &prm_from_phi,   ///< The phi angle of the from_residue
                                               const angle_type  &prm_from_psi,   ///< The psi angle of the from_residue
                                               const angle_type  &prm_to_phi,     ///< The phi angle of the to_residue
                                               const angle_type  &prm_to_psi      ///< The psi angle of the to_residue
                                               ) : view           ( std::move( prm_view      ) ),
                                                   frame          (
                                                   	make_quat_rot_from_rotation<frame_quat_rot_type>(
                                                   		prm_frame
                                                   	)
                                                   ),
                                                   from_phi_angle ( shift_copy( prm_from_phi ) ),
                                                   from_psi_angle ( shift_copy( prm_from_psi ) ),
                                                   to_phi_angle   ( shift_copy( prm_to_phi   ) ),
                                                   to_psi_angle   ( shift_copy( prm_to_psi   ) ),
                                                   from_index     ( prm_from_index ),
                                                   to_index       ( prm_to_index   ) {
}

/// \brief TODOCUMENT
///
/// \relates view_cache_index_entry
ostream & cath::index::operator<<(ostream                      &prm_os,   ///< TODOCUMENT
                                  const view_cache_index_entry &prm_cache ///< TODOCUMENT
                                  ) {
  prm_os << "view_cache_index_entry[from-to-indices:";
  prm_os << prm_cache.get_from_index();
  prm_os << ",";
  prm_os << prm_cache.get_to_index();
  prm_os << ", view:[";
  prm_os << prm_cache.get_view().get<0>();
  prm_os << ", ";
  prm_os << prm_cache.get_view().get<1>();
  prm_os << ", ";
  prm_os << prm_cache.get_view().get<2>();
  prm_os << "], frame:";
  prm_os << prm_cache.get_frame();
  prm_os << ", from_phi:";
  prm_os << prm_cache.get_from_phi_angle();
  prm_os << ", from_psi:";
  prm_os << prm_cache.get_from_psi_angle();
  prm_os << ", to_phi:";
  prm_os << prm_cache.get_to_phi_angle();
  prm_os << ", to_psi:";
  prm_os << prm_cache.get_to_psi_angle();
  prm_os << "]";
  return prm_os;
}

/// \brief Factory function to build a view_cache_index_entry from a protein and the from/to indices
///
/// \relates view_cache_index_entry
view_cache_index_entry cath::index::detail::make_view_cache_index_entry(const protein &prm_protein,    ///< The protein containing the residues the view_cache_index_entry should represent
                                                                        const size_t  &prm_from_index, ///< The index of the from_residue the view_cache_index_entry should represent
                                                                        const size_t  &prm_to_index    ///< The index of the to_residue the view_cache_index_entry should represent
                                                                        ) {
	return view_cache_index_entry(
		numeric_cast<index_type>( prm_from_index ),
		numeric_cast<index_type>( prm_to_index   ),
		view_vector( prm_protein, prm_from_index, prm_to_index ),
		view_frame ( prm_protein, prm_from_index, prm_to_index ),
		convert_angle_type<angle_base_type>( prm_protein.get_residue_ref_of_index( prm_from_index ).get_phi_angle() ),
		convert_angle_type<angle_base_type>( prm_protein.get_residue_ref_of_index( prm_from_index ).get_psi_angle() ),
		convert_angle_type<angle_base_type>( prm_protein.get_residue_ref_of_index( prm_to_index   ).get_phi_angle() ),
		convert_angle_type<angle_base_type>( prm_protein.get_residue_ref_of_index( prm_to_index   ).get_psi_angle() )
	);
}
