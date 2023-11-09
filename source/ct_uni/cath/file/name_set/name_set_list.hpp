/// \file
/// \brief The name_set_list class header

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

#ifndef CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_NAME_SET_NAME_SET_LIST_HPP
#define CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_NAME_SET_NAME_SET_LIST_HPP

#include "cath/file/file_type_aliases.hpp"
#include "cath/file/name_set/name_set.hpp"

namespace cath::file {

	/// \brief A list of name_set objects
	class name_set_list final {
	private:
		/// \brief The vector of name_sets
		name_set_vec name_sets;

	public:
		/// \brief A iterator type alias as part of making this a range over name_sets
		using iterator = name_set_vec::iterator;

		/// \brief A const_iterator type alias as part of making this a range over name_sets
		using const_iterator = name_set_vec::const_iterator;

		/// \brief Ctor
		name_set_list() noexcept = default;
		explicit name_set_list(name_set_vec);
		explicit name_set_list(const size_t &);

		[[nodiscard]] bool   empty() const;
		[[nodiscard]] size_t size() const;

		name_set & operator[](const size_t &);
		const name_set & operator[](const size_t &) const;

		iterator begin();
		iterator end();
		[[nodiscard]] const_iterator begin() const;
		[[nodiscard]] const_iterator end() const;
	};

	/// \brief Ctor from vector of name_sets
	inline name_set_list::name_set_list(name_set_vec prm_name_sets ///< The vector of name_sets from which to construct this name_set_list
	                                    ) : name_sets{ std::move( prm_name_sets ) } {
	}

	/// \brief Ctor from a number of empty name_sets to create
	inline name_set_list::name_set_list(const size_t &prm_size
	                                    ) : name_sets{ name_set_vec( prm_size ) } {
	}

	/// \brief Return whether this is empty
	inline bool name_set_list::empty() const {
		return name_sets.empty();
	}

	/// \brief Return the number of name_sets
	inline size_t name_set_list::size() const {
		return name_sets.size();
	}

	/// \brief Return the name_set stored at the specified index
	inline name_set & name_set_list::operator[](const size_t &prm_index ///< The index of the name_set to return
	                                            ) {
		return name_sets[ prm_index ];
	}

	/// \brief Return the name_set stored at the specified index
	inline const name_set & name_set_list::operator[](const size_t &prm_index ///< The index of the name_set to return
	                                                  ) const {
		return name_sets[ prm_index ];
	}

	/// \brief Standard begin() method, as part of making this into a range over the name_sets
	inline auto name_set_list::begin() -> iterator {
		return std::begin( name_sets );
	}

	/// \brief Standard begin() method, as part of making this into a range over the name_sets
	inline auto name_set_list::end() -> iterator {
		return std::end( name_sets );
	}

	/// \brief Standard const begin() method, as part of making this into a range over the name_sets
	inline auto name_set_list::begin() const -> const_iterator {
		return ::std::cbegin( name_sets );
	}

	/// \brief Standard const end() method, as part of making this into a range over the name_sets
	inline auto name_set_list::end() const -> const_iterator {
		return ::std::cend  ( name_sets );
	}

	name_set_list build_name_set_list(str_vec,
	                                  str_vec = str_vec{},
	                                  str_opt_vec = str_opt_vec{});

	std::string to_string(const name_set_list &);

	std::ostream & operator<<(std::ostream &,
	                          const name_set_list &);

	bool all_have_specified_id(const name_set_list &);
	bool all_have_domain_name_from_regions(const name_set_list &);

	str_vec get_names_from_acq(const name_set_list &);
	str_vec get_domain_or_specified_or_from_acq_names(const name_set_list &);

	str_vec get_alignment_html_names(const name_set_list &);
	str_vec get_multi_ssap_alignment_file_names(const name_set_list &);
	str_vec get_protein_list_names(const name_set_list &);
	str_vec get_supn_json_names(const name_set_list &);
	str_vec get_supn_pdb_file_names(const name_set_list &);
	str_vec get_viewer_names(const name_set_list &);

	void add_specified_ids(name_set_list &,
	                       str_vec);
	name_set_list add_specified_ids_copy(name_set_list,
	                                     str_vec);

	void add_domain_names_from_regions(name_set_list &,
	                                   str_opt_vec);
	name_set_list add_domain_names_from_regions_copy(name_set_list,
	                                                 str_opt_vec);

} // namespace cath::file

#endif // CATH_TOOLS_SOURCE_CT_UNI_CATH_FILE_NAME_SET_NAME_SET_LIST_HPP
