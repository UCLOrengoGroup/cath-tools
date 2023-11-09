/// \file
/// \brief The residue_id class header

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

#ifndef CATH_TOOLS_SOURCE_CT_BIOCORE_CATH_BIOCORE_RESIDUE_ID_HPP
#define CATH_TOOLS_SOURCE_CT_BIOCORE_CATH_BIOCORE_RESIDUE_ID_HPP

#include <istream>
#include <ostream>
#include <string>
#include <utility>

#include <boost/operators.hpp>

#include "cath/biocore/biocore_type_aliases.hpp"
#include "cath/biocore/chain_label.hpp"
#include "cath/biocore/residue_name.hpp"

namespace cath {

	/// \brief A PDB residue identifier, bundling the chain_label and residue_name
	///
	/// Converted to a string, this might look like "A:324C"
	class residue_id final : private boost::equality_comparable<residue_id> {
	  private:
		/// \brief The chain on which the residue belongs
		chain_label chain;

		/// \brief The residue_name in the PDB
		residue_name res_name;

	  public:
		constexpr residue_id() = default;
		constexpr residue_id( const chain_label &, residue_name );

		[[nodiscard]] constexpr const chain_label & get_chain_label() const;
		[[nodiscard]] constexpr const residue_name &get_residue_name() const;

		/// \brief Return whether the two specified residue_ids are identical
		///
		/// \relates residue_id
		///
		/// \param prm_lhs The first  residue_id to compare
		/// \param prm_rhs The second residue_id to compare
		friend constexpr bool operator==( const residue_id &prm_lhs, const residue_id &prm_rhs ) {
			// clang-format off
			return (
				prm_lhs.get_chain_label() == prm_rhs.get_chain_label()
				&&
				prm_lhs.get_residue_name() == prm_rhs.get_residue_name()
			);
			// clang-format on
		}
	};

	/// \brief Ctor for residue_id
	///
	/// \param prm_chain    The chain on which the residue belongs
	/// \param prm_res_name The residue_name in the PDB
	constexpr residue_id::residue_id( const chain_label &prm_chain, residue_name prm_res_name ) :
	        chain{ prm_chain }, res_name{ std::move( prm_res_name ) } {
	}

	/// \brief Getter for the chain on which the residue belongs
	constexpr const chain_label & residue_id::get_chain_label() const {
		return chain;
	}

	/// \brief Getter for The residue_name in the PDB
	constexpr const residue_name &residue_id::get_residue_name() const {
		return res_name;
	}

	/// \brief Whether this is a null residue_id
	///
	/// \param prm_residue_id The residue_id to query
	constexpr bool is_null( const residue_id &prm_residue_id ) {
		return prm_residue_id.get_residue_name().is_null();
	}

	/// \brief Make a residue_id from the chain code
	///
	/// \relates residue_id
	///
	/// \param prm_chain_char The chain code
	constexpr residue_id make_residue_id( const char &prm_chain_char ) {
		return { chain_label{ prm_chain_char }, residue_name{} };
	}

	/// \brief Make a residue_id from the chain code and residue number
	///
	/// \relates residue_id
	///
	/// \param prm_chain_char     The chain code
	/// \param prm_residue_number The residue number
	constexpr residue_id make_residue_id( const char &prm_chain_char, const int &prm_residue_number ) {
		return { chain_label{ prm_chain_char }, residue_name{ prm_residue_number } };
	}

	/// \brief Make a residue_id from the chain code, residue number and insert code
	///
	/// \relates residue_id
	///
	/// \param prm_chain_char     The chain code
	/// \param prm_residue_number The residue number
	/// \param prm_insert_code    The insert code
	constexpr residue_id make_residue_id( const char &prm_chain_char, const int &prm_residue_number, const char &prm_insert_code ) {
		return { chain_label{ prm_chain_char }, residue_name{ prm_residue_number, prm_insert_code } };
	}

	/// \brief Return whether the specified residue ID has a strictly negative residue number
	///
	/// \relates residue_id
	///
	/// \param prm_residue_id The residue_id to query
	constexpr bool has_strictly_negative_residue_number( const residue_id &prm_residue_id ) {
		return has_strictly_negative_residue_number( prm_residue_id.get_residue_name() );
	}

	std::string to_string(const residue_id &);

	std::ostream & operator<<(std::ostream &,
	                          const residue_id &);

	std::istream & operator>>(std::istream &,
	                          residue_id &);

	residue_id make_residue_id(const std::string &);

	chain_label_residue_id_vec_map get_residue_id_by_chain_label(const residue_id_vec &);

	chain_label_opt consistent_chain_label(const residue_id_vec &);
	bool have_consistent_chain_labels(const residue_id_vec &);
	bool has_any_strictly_negative_residue_numbers(const residue_id_vec &);

} // namespace cath

#endif // CATH_TOOLS_SOURCE_CT_BIOCORE_CATH_BIOCORE_RESIDUE_ID_HPP
