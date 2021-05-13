/// \file
/// \brief The residue_name class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_BIOCORE_CATH_BIOCORE_RESIDUE_NAME_HPP
#define _CATH_TOOLS_SOURCE_CT_BIOCORE_CATH_BIOCORE_RESIDUE_NAME_HPP

#include <optional>

#include <boost/operators.hpp>

namespace cath {

	/// \brief Represent a PDB residue name (eg 324A)
	///
	/// TODO: Come C++17, make this constexpr
	class residue_name final : private boost::equality_comparable<residue_name> {
	private:
		/// \brief The residue number
		int res_num = 0;

		/// \brief The (optional insert code)
		::std::optional<char> insert;

		/// \brief Whether this is a null residue
		bool is_null_residue_name = true;

		void sanity_check() const;
		void sanity_check_is_not_null_residue() const;

	public:
		residue_name();
		explicit residue_name(const int &);
		residue_name(const int &,
		             const char &);

		const bool & is_null() const;
		const int & residue_number() const;
		const ::std::optional<char> & opt_insert() const;
	};


	bool operator< (const residue_name &, const residue_name &) = delete; // Better to avoid residue_names (except for equality/inequality) because disallowing it avoids anyone mistakenly trying to do that to get the correct order of residues
	bool operator<=(const residue_name &, const residue_name &) = delete; // Better to avoid residue_names (except for equality/inequality) because disallowing it avoids anyone mistakenly trying to do that to get the correct order of residues
	bool operator> (const residue_name &, const residue_name &) = delete; // Better to avoid residue_names (except for equality/inequality) because disallowing it avoids anyone mistakenly trying to do that to get the correct order of residues
	bool operator>=(const residue_name &, const residue_name &) = delete; // Better to avoid residue_names (except for equality/inequality) because disallowing it avoids anyone mistakenly trying to do that to get the correct order of residues

	bool operator==(const residue_name &,
	                const residue_name &);

	std::string to_string(const residue_name &);
	std::ostream & operator<<(std::ostream &,
	                          const residue_name &);

	std::istream & operator>>(std::istream &,
	                          residue_name &);

	int residue_number_or_value_if_null(const residue_name &,
	                                    const int &);
	::std::optional<char> opt_insert_or_value_if_null(const residue_name &,
	                                                  const ::std::optional<char> &);
	bool has_insert(const residue_name &);
	bool has_insert_or_value_if_null(const residue_name &,
	                                 const bool &);
	const char & insert(const residue_name &);
	char insert_or_value_if_null(const residue_name &,
	                             const char &);
	char insert_or_value_if_null_or_absent(const residue_name &,
	                                       const char &);
	bool has_strictly_negative_residue_number(const residue_name &);

	std::string insert_string(const residue_name &);
	std::string make_residue_name_string_with_insert_or_space(const residue_name &);
	residue_name make_residue_name(const std::string &);
	residue_name make_residue_name_with_non_insert_char(const int &,
	                                                    const char &,
	                                                    const char & = ' ');


} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_BIOCORE_CATH_BIOCORE_RESIDUE_NAME_HPP
