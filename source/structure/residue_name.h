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

#ifndef RESIDUE_NAME_H_INCLUDED
#define RESIDUE_NAME_H_INCLUDED

#include <boost/operators.hpp>
#include <boost/optional.hpp>

namespace cath {

	/// \brief TODOCUMENT
	class residue_name final : private boost::equality_comparable<residue_name> {
	private:
		/// \brief TODOCUMENT
		bool is_null_residue_name = true;

		/// \brief TODOCUMENT
		int residue_number = 0;

		/// \brief TODOCUMENT
		boost::optional<char> insert_code;

		void sanity_check() const;
		void sanity_check_is_not_null_residue() const;

	public:
		residue_name();
		explicit residue_name(const int &);
		residue_name(const int &,
		             const char &);

		const bool & get_is_null_residue_name() const;
		const int & get_residue_number() const;
		const boost::optional<char> & get_opt_insert_code() const;
	};

	bool operator==(const residue_name &,
	                const residue_name &);

	std::string to_string(const residue_name &);
	std::ostream & operator<<(std::ostream &,
	                          const residue_name &);

	std::istream & operator>>(std::istream &,
	                          residue_name &);

	bool has_insert_code(const residue_name &);
	char get_insert_code(const residue_name &);
	std::string get_insert_code_string(const residue_name &);
	std::string make_residue_name_string_with_insert_or_space(const residue_name &);
	residue_name make_residue_name(const std::string &);
	residue_name make_residue_name_with_non_insert_char(const int &,
	                                                    const char &,
	                                                    const char &arg_non_insert_char = ' ');


}

#endif
