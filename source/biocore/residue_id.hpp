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

#ifndef _CATH_TOOLS_SOURCE_BIOCORE_RESIDUE_ID_H
#define _CATH_TOOLS_SOURCE_BIOCORE_RESIDUE_ID_H

#include "biocore/chain_label.hpp"
#include "biocore/residue_name.hpp"
#include "structure/structure_type_aliases.hpp"

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
		residue_id() = default;
		residue_id(const chain_label &,
		           residue_name);

		const chain_label & get_chain_label() const;
		const residue_name & get_residue_name() const;
	};

	/// \brief Ctor for residue_id
	inline residue_id::residue_id(const chain_label &arg_chain,   ///< The chain on which the residue belongs
	                              residue_name       arg_res_name ///< The residue_name in the PDB
	                              ) : chain   { arg_chain                 },
	                                  res_name{ std::move( arg_res_name ) } {
	}

	/// \brief Getter for the chain on which the residue belongs
	inline const chain_label & residue_id::get_chain_label() const {
		return chain;
	}

	/// \brief Getter for The residue_name in the PDB
	inline const residue_name & residue_id::get_residue_name() const {
		return res_name;
	}

	/// \brief Whether this is a null residue_id
	inline bool is_null(const residue_id &arg_residue_id ///< The residue_id to query
	                    ) {
		return arg_residue_id.get_residue_name().is_null();
	}

	/// \brief Return whether the two specified residue_ids are identical
	///
	/// \relates residue_id
	inline bool operator==(const residue_id &arg_residue_id_a, ///< The first  residue_id to compare
	                       const residue_id &arg_residue_id_b  ///< The second residue_id to compare
	                       ) {
		return (
			arg_residue_id_a.get_chain_label()  == arg_residue_id_b.get_chain_label()
			&&
			arg_residue_id_a.get_residue_name() == arg_residue_id_b.get_residue_name()
		);
	}

	/// \brief Make a residue_id from the chain code
	///
	/// \relates residue_id
	inline residue_id make_residue_id(const char &arg_chain_char ///< The chain code
	                                  ) {
		return { chain_label{ arg_chain_char }, residue_name{} };
	}

	/// \brief Make a residue_id from the chain code and residue number
	///
	/// \relates residue_id
	inline residue_id make_residue_id(const char &arg_chain_char,    ///< The chain code
	                                  const int  &arg_residue_number ///< The residue number
	                                  ) {
		return { chain_label{ arg_chain_char }, residue_name{ arg_residue_number } };
	}

	/// \brief Make a residue_id from the chain code, residue number and insert code
	///
	/// \relates residue_id
	inline residue_id make_residue_id(const char &arg_chain_char,     ///< The chain code
	                                  const int  &arg_residue_number, ///< The residue number
	                                  const char &arg_insert_code     ///< The insert code
	                                  ) {
		return { chain_label{ arg_chain_char }, residue_name{ arg_residue_number, arg_insert_code } };
	}

	std::string to_string(const residue_id &);

	std::ostream & operator<<(std::ostream &,
	                          const residue_id &);

	std::istream & operator>>(std::istream &,
	                          residue_id &);

	residue_id make_residue_id(const std::string &);

	bool has_strictly_negative_residue_number(const residue_id &);

	chain_label_residue_id_vec_map get_residue_id_by_chain_label(const residue_id_vec &);

	chain_label_opt consistent_chain_label(const residue_id_vec &);
	bool have_consistent_chain_labels(const residue_id_vec &);
	bool has_any_strictly_negative_residue_numbers(const residue_id_vec &);

} // namespace cath

#endif
