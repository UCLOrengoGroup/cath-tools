/// \file
/// \brief The residue class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_RESIDUE_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_RESIDUE_HPP

#include <cstddef>
#include <iosfwd>
#include <string_view>

#include <boost/operators.hpp>

#include "cath/biocore/residue_id.hpp"
#include "cath/common/exception/invalid_argument_exception.hpp"
#include "cath/file/pdb/dssp_skip_policy.hpp"
#include "cath/file/pdb/residue_makeup.hpp"
#include "cath/structure/geometry/angle.hpp"
#include "cath/structure/geometry/coord.hpp"
#include "cath/structure/geometry/rotation.hpp"
#include "cath/structure/protein/amino_acid.hpp"
#include "cath/structure/protein/sec_struc_type.hpp"

namespace cath {

	/// \brief TODOCUMENT
	///
	/// TODO: Make this constexpr
	///
	/// \todo This currently hard-codes the member+getter for two atoms coordinates: CA and CB
	///       These two should remain and remain mandatory but they should be complemented by
	///       N, C and O but these should be optional (with has_atom_coord method).
	///       This looks like a strong candidate for a template-based solution...
	///       Perhaps something like a Boost Fusion container that stores the coords based on enum values
	///       and template method accessors like `the_residue.has_atom_coord<ATOM_O>()` and
	///       `the_residue.get_atom_coord<ATOM_O>()`, and some static method of detecting which values
	///       are mandatory and which optional.
	///       The code that does this should be encapsulated in an implementation class because this class
	///       may be independently useful, eg for potentially applying common_atom_selection_policy to
	///       altering the distances used in SSAP scores.
	class residue final : private boost::equality_comparable<residue> {
	private:
		/// \brief TODOCUMENT
		residue_id     the_residue_id;

		/// \brief Single character representing the amino acid ('0' or an upper case letter)
		amino_acid     the_amino_acid;

		/// \brief Coordinates of the carbon alpha atom
		geom::coord    carbon_alpha_coord;

		/// \brief Coordinates of the carbon beta atom
		geom::coord    carbon_beta_coord;

		/// \brief TODOCUMENT
		size_t         sec_struc_number;

		/// \brief TODOCUMENT
		sec_struc_type the_sec_struc_type; ///< \todo Check: is it a problem that the sec files have values of H/S but the source compares this to H/E.
		                                   ///<       At least sec_struc.get_type() is used to check that one sec_struc's type matches another.

		/// \brief TODOCUMENT
		geom::rotation   frame;

		/// \brief The phi angle of this residue, in degrees within (0, 360.0] (possibly [0, 360.0] ? )
		geom::doub_angle phi_angle;

		/// \brief The psi angle of this residue, in degrees within (0, 360.0] (possibly [0, 360.0] ? )
		geom::doub_angle psi_angle;

		/// \brief TODOCUMENT
		size_t           access;

		constexpr static void check_phi_psi_angle(const geom::doub_angle &);

	  public:
		constexpr residue( const residue_id &,
		                   amino_acid,
		                   geom::coord,
		                   geom::coord,
		                   const size_t &,
		                   const sec_struc_type &,
		                   geom::rotation,
		                   geom::doub_angle,
		                   geom::doub_angle,
		                   const size_t & );

		void initialise();

		constexpr residue & set_amino_acid(const amino_acid &);
		constexpr residue & set_residue_sec_struc_number(const size_t &);
		constexpr residue & set_sec_struc_type(const sec_struc_type &);
		constexpr residue & set_access(const size_t &);

		[[nodiscard]] constexpr const residue_id &        get_pdb_residue_id() const;
		[[nodiscard]] constexpr amino_acid                get_amino_acid() const;
		[[nodiscard]] constexpr inline const geom::coord &get_carbon_alpha_coord() const;
		[[nodiscard]] constexpr inline const geom::coord &get_carbon_beta_coord() const;

		[[nodiscard]] constexpr size_t         get_sec_struc_number() const;
		[[nodiscard]] constexpr sec_struc_type get_sec_struc_type() const;

		[[nodiscard]] constexpr inline const geom::rotation &get_frame() const;
		[[nodiscard]] constexpr const geom::doub_angle &     get_phi_angle() const;
		[[nodiscard]] constexpr const geom::doub_angle &     get_psi_angle() const;
		[[nodiscard]] constexpr size_t                       get_access() const;

		/// \brief The default angles in degrees that's used for undetermined phi/psi angles
		///
		/// \todo Shouldn't the undetermined angles be handled more explicitly?
		///       The DSSP/WOLF files have the angles in the more natural range of [-180, 180]
		///       and the currently always get shifted into the (0, 360] range
		///       if they were left alone, that would leave 360.0 as a special "undetermined" value
		///       but this would require changes in residues_have_similar_area_angle_props()
		static constexpr geom::doub_angle DEFAULT_PHI_PSI = geom::ONE_REVOLUTION<double>;

		/// \brief Equality operator for residue
		///
		/// operator!= is provided by boost::equality_comparable
		///
		/// \relates residue
		///
		/// \param prm_lhs The first  residue to compare
		/// \param prm_rhs The second residue to compare
		friend constexpr bool operator==( const residue &prm_lhs, const residue &prm_rhs ) {
			// clang-format off
			return (
				( prm_lhs.get_pdb_residue_id()     == prm_rhs.get_pdb_residue_id()     )
				&&
				( prm_lhs.get_amino_acid()         == prm_rhs.get_amino_acid()         )
				&&
				( prm_lhs.get_carbon_alpha_coord() == prm_rhs.get_carbon_alpha_coord() )
				&&
				( prm_lhs.get_carbon_beta_coord()  == prm_rhs.get_carbon_beta_coord()  )
				&&
				( prm_lhs.get_sec_struc_number()   == prm_rhs.get_sec_struc_number()   )
				&&
				( prm_lhs.get_sec_struc_type()     == prm_rhs.get_sec_struc_type()     )
				&&
				( prm_lhs.get_frame()              == prm_rhs.get_frame()              )
				&&
				( prm_lhs.get_phi_angle()          == prm_rhs.get_phi_angle()          )
				&&
				( prm_lhs.get_psi_angle()          == prm_rhs.get_psi_angle()          )
				&&
				( prm_lhs.get_access()             == prm_rhs.get_access()             )
			);
			// clang-format off
		}

	};

	/// \brief Throw if angle is not in range [0, 360]
	constexpr void residue::check_phi_psi_angle( const geom::doub_angle &prm_angle ) {
		//	using ::boost::math::isfinite;
		if ( prm_angle < geom::ZERO_ANGLE<double> || prm_angle > geom::ONE_REVOLUTION<double> ) {
			BOOST_THROW_EXCEPTION( common::invalid_argument_exception(
			  "Phi/psi angle (" + ::boost::lexical_cast<::std::string>( prm_angle ) + ") is not in valid range [0, 360]" ) );
		}
	}

	/// \brief Ctor for residue
	///
	/// \param prm_residue_id         TODOCUMENT
	/// \param prm_amino_acid         TODOCUMENT
	/// \param prm_carbon_alpha_coord Coordinates of the carbon alpha atom
	/// \param prm_carbon_beta_coord  Coordinates of the carbon beta atom
	/// \param prm_sec_struc_number   TODOCUMENT
	/// \param prm_sec_struc_type     TODOCUMENT
	/// \param prm_frame              TODOCUMENT
	/// \param prm_phi                TODOCUMENT
	/// \param prm_psi                TODOCUMENT
	/// \param prm_access             TODOCUMENT
	constexpr residue::residue( const residue_id &    prm_residue_id,
	                            amino_acid            prm_amino_acid,
	                            geom::coord           prm_carbon_alpha_coord,
	                            geom::coord           prm_carbon_beta_coord,
	                            const size_t &        prm_sec_struc_number,
	                            const sec_struc_type &prm_sec_struc_type,
	                            geom::rotation        prm_frame,
	                            geom::doub_angle      prm_phi,
	                            geom::doub_angle      prm_psi,
	                            const size_t &        prm_access ) :
	        the_residue_id( prm_residue_id ),
	        the_amino_acid( std::move( prm_amino_acid ) ),
	        carbon_alpha_coord( std::move( prm_carbon_alpha_coord ) ),
	        carbon_beta_coord( std::move( prm_carbon_beta_coord ) ),
	        sec_struc_number( prm_sec_struc_number ),
	        the_sec_struc_type( prm_sec_struc_type ),
	        frame( std::move( prm_frame ) ),
	        phi_angle( std::move( prm_phi ) ),
	        psi_angle( std::move( prm_psi ) ),
	        access( prm_access ) {
		check_phi_psi_angle( get_phi_angle() );
		check_phi_psi_angle( get_psi_angle() );
	}

	/// \brief TODOCUMENT
	constexpr const geom::coord & residue::get_carbon_alpha_coord() const {
		return carbon_alpha_coord;
	}

	/// \brief TODOCUMENT
	constexpr const geom::coord & residue::get_carbon_beta_coord() const {
		return carbon_beta_coord;
	}

	/// \brief TODOCUMENT
	constexpr const geom::rotation & residue::get_frame() const {
		return frame;
	}

	/// \brief Setter for the amino acid
	constexpr residue & residue::set_amino_acid(const amino_acid &prm_amino_acid ///< The amino acid to set
	                                            ) {
		the_amino_acid = prm_amino_acid;
		return *this;
	}

	/// \brief Setter for the residue sec_struc number
	constexpr residue & residue::set_residue_sec_struc_number(const size_t &prm_sec_struc_number ///< The residue sec_struc number to set
	                                                          ) {
		sec_struc_number = prm_sec_struc_number;
		return *this;
	}

	/// \brief Setter for the sec_struc_type
	constexpr residue & residue::set_sec_struc_type(const sec_struc_type &prm_sec_struc_type ///< The sec_struc_type to set
	                                                ) {
		the_sec_struc_type = prm_sec_struc_type;
		return *this;
	}

	/// \brief Setter for the accessibility (calculated in a DSSP/wolf manner)
	constexpr residue & residue::set_access(const size_t &prm_accessibility ///< The accessibility (calculated in a DSSP/wolf manner) to set
	                                        ) {
		access = prm_accessibility;
		return *this;
	}

	/// \brief TODOCUMENT
	constexpr const residue_id & residue::get_pdb_residue_id() const {
		return the_residue_id;
	}

	/// \brief TODOCUMENT
	constexpr amino_acid residue::get_amino_acid() const {
		return the_amino_acid;
	}

	/// \brief TODOCUMENT
	constexpr size_t residue::get_sec_struc_number() const {
		return sec_struc_number;
	}

	/// \brief TODOCUMENT
	constexpr sec_struc_type residue::get_sec_struc_type() const {
		return the_sec_struc_type;
	}

	/// \brief TODOCUMENT
	constexpr const geom::doub_angle & residue::get_phi_angle() const {
		return phi_angle;
	}

	/// \brief TODOCUMENT
	constexpr const geom::doub_angle & residue::get_psi_angle() const {
		return psi_angle;
	}

	/// \brief TODOCUMENT
	constexpr size_t residue::get_access() const {
		return access;
	}

	/// \brief TODOCUMENT
	///
	/// \relates residue
	constexpr const chain_label & get_chain_label(const residue &prm_residue ///< TODOCUMENT
	                                              ) {
		return prm_residue.get_pdb_residue_id().get_chain_label();
	}

	/// \brief TODOCUMENT
	///
	/// \relates residue
	constexpr const residue_name & get_pdb_residue_name(const residue &prm_residue ///< TODOCUMENT
	                                                    ) {
		return prm_residue.get_pdb_residue_id().get_residue_name();
	}


	/// \brief TODOCUMENT
	///
	/// \relates residue
	constexpr const int & pdb_number(const residue &prm_residue ///< TODOCUMENT
	                                 ) {
		return get_pdb_residue_name( prm_residue ).residue_number();
	}

	/// \brief TODOCUMENT
	///
	/// \relates residue
	constexpr int pdb_number_or_value_if_null(const residue &prm_residue, ///< TODOCUMENT
	                                          const int     &prm_value
	                                          ) {
		return residue_number_or_value_if_null( get_pdb_residue_name( prm_residue ), prm_value );
	}

	/// \brief TODOCUMENT
	///
	/// \relates residue
	constexpr bool has_pdb_insert(const residue &prm_residue ///< TODOCUMENT
	                              ) {
		return           has_insert( get_pdb_residue_name( prm_residue ) );
	}

	/// \brief TODOCUMENT
	///
	/// \relates residue
	constexpr bool has_pdb_insert_or_value_if_null(const residue &prm_residue, ///< TODOCUMENT
	                                               const bool    &prm_value    ///< TODOCUMENT
	                                               ) {
		return has_insert_or_value_if_null( get_pdb_residue_name( prm_residue ), prm_value );
	}

	/// \brief TODOCUMENT
	///
	/// \relates residue
	constexpr const char & pdb_insert(const residue &prm_residue ///< TODOCUMENT
	                                  ) {
		if ( ! has_pdb_insert( prm_residue ) ) {
			BOOST_THROW_EXCEPTION(common::invalid_argument_exception("Cannot get PDB residue insert code of residue that has no insert code"));
		}
		return insert( get_pdb_residue_name( prm_residue ) );
	}

	/// \brief TODOCUMENT
	///
	/// \relates residue
	constexpr char pdb_insert_or_value_if_null(const residue &prm_residue, ///< TODOCUMENT
	                                           const char    &prm_value    ///< TODOCUMENT
	                                           ) {
		return insert_or_value_if_null( get_pdb_residue_name( prm_residue ), prm_value );
	}

	/// \brief TODOCUMENT
	///
	/// \relates residue
	constexpr char pdb_insert_or_value_if_null_or_absent(const residue &prm_residue, ///< TODOCUMENT
	                                                     const char    &prm_value    ///< TODOCUMENT
	                                                     ) {
		return insert_or_value_if_null_or_absent( get_pdb_residue_name( prm_residue ), prm_value );
	}

	/// \brief TODOCUMENT
	///
	/// \relates residue
	constexpr bool residue_matches_residue_id(const residue    &prm_residue,   ///< TODOCUMENT
	                                          const residue_id &prm_residue_id ///< TODOCUMENT
	                                          ) {
		return ( prm_residue.get_pdb_residue_id() == prm_residue_id );
	}

	/// \brief TODOCUMENT
	///
	/// \relates residue
	constexpr char get_amino_acid_letter_tolerantly(const residue &prm_residue ///< TODOCUMENT
	                                                ) {
		return prm_residue.get_amino_acid().get_letter_tolerantly();
	}

	/// \brief TODOCUMENT
	///
	/// \relates residue
	constexpr char_3_arr get_amino_acid_code(const residue &prm_residue ///< TODOCUMENT
	                                         ) {
		return prm_residue.get_amino_acid().get_code();
	}

	/// \brief TODOCUMENT
	///
	/// \relates residue
	constexpr ::std::string_view get_amino_acid_name(const residue &prm_residue ///< TODOCUMENT
	                                                 ) {
		return prm_residue.get_amino_acid().get_name();
	}

	/// \brief Wipe the residue's secondary structure number and set its label to sec_struc_type::COIL
	///
	/// \relates residue
	constexpr void wipe_secondary_structure(residue &prm_residue ///< The residue to be modified
	                                        ) {
		prm_residue.set_sec_struc_type          ( sec_struc_type::COIL );
		prm_residue.set_residue_sec_struc_number( 0                    );
	}

	/// \brief TODOCUMENT
	///
	/// \todo Sort out indexing and then get rid of this because a blank residue doesn't really make any sense
	inline constexpr residue NULL_RESIDUE( residue_id{},
	                                       amino_acid( 'X' ),
	                                       geom::ORIGIN_COORD,
	                                       geom::ORIGIN_COORD,
	                                       0, // Secondary structure number
	                                       sec_struc_type::COIL,
	                                       geom::IDENTITY_ROTATION,
	                                       residue::DEFAULT_PHI_PSI,
	                                       residue::DEFAULT_PHI_PSI,
	                                       0 // Accessibility
	);

	/// \brief Check whether the specified residue is a null residue
	///
	/// This currently does a full equality comparison with NULL_RESIDUE but
	/// could probably be made more efficient if there is call for that.
	///
	/// \relates residue
	constexpr bool is_null_residue(const residue &prm_residue ///< TODOCUMENT
	                               ) {
		return prm_residue == NULL_RESIDUE;
	}

	std::string ssap_legacy_alignment_left_side_string(const cath::residue &);

	std::string ssap_legacy_alignment_left_side_gap_string();

	std::string ssap_legacy_alignment_right_side_string(const cath::residue &);

	std::string ssap_legacy_alignment_right_side_gap_string();

	std::string get_pdb_residue_id_string(const residue &);

	int get_accessi_of_residue(const residue &);

	std::ostream & operator<<(std::ostream &,
	                          const residue &);

	geom::rotation construct_residue_frame(const geom::coord &,
	                                       const geom::coord &,
	                                       const geom::coord &);

	residue combine_residues_from_dssp_and_pdb(const residue &,
	                                           const residue &,
	                                           const file::dssp_skip_angle_skipping & = file::dssp_skip_angle_skipping::BREAK_ANGLES,
	                                           const file::residue_makeup           & = file::residue_makeup::ALL_PROPER_AMINO_ACIDS);

} // namespace cath
#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_STRUCTURE_PROTEIN_RESIDUE_HPP
