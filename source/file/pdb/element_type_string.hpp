/// \file
/// \brief The element_type_string class header

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

#ifndef _CATH_TOOLS_SOURCE_FILE_PDB_ELEMENT_TYPE_STRING_H
#define _CATH_TOOLS_SOURCE_FILE_PDB_ELEMENT_TYPE_STRING_H

#include <boost/utility/string_ref.hpp>

#include "common/string/string_parse_tools.hpp"
#include "file/pdb/coarse_element_type.hpp"

#include <string>

namespace cath {
	namespace file {

		/// \brief Hold an element type string along with a string_ref to the non-whitespace
		///        part of the string.
		///
		/// The point of this class is to bundle the string_ref and the string to which it refers together
		/// to ensure that the copy/move construction/assignment is done correctly.
		/// In all four cases, that requires the string_ref to be reset to point to the new string.
		///
		/// \todo If there are any other uses for bundling a string_ref with its string, then
		///       consider generalising this class.
		///
		/// \todo Come C++17 replace boost::string_ref with std::string_view
		class element_type_string final {
		private:
			/// \brief The untrimmed string describing the element type of this atom.
			///
			/// This is whitespace-untrimmed so that the correct whitespace can be returned whilst writing
			std::string element_type_untrimmed;

			/// \brief The trimmed string describing the element type of this atom.
			///
			/// \todo Come C++17, replace boost::string_ref with std::string_view
			boost::string_ref element_type;

			element_type_string(const element_type_string &&,
			                    const ptrdiff_t &,
			                    const size_t &);

		public:
			element_type_string(const std::string &);
			element_type_string(const std::string &&);

			element_type_string(const element_type_string &);
			element_type_string(element_type_string &&) noexcept;
			element_type_string & operator=(const element_type_string &);
			element_type_string & operator=(element_type_string &&) noexcept;

			const std::string & get_element_type_untrimmed() const;
			const boost::string_ref & get_element_type() const;
		};

		/// \brief Private delegation-ctor used to implement the move-ctor properly
		inline element_type_string::element_type_string(const element_type_string &&arg_rhs,   ///< The element_type_string from which to move-construct
		                                                const ptrdiff_t            &arg_start, ///< The start offset of the string_ref
		                                                const size_t               &arg_length ///< The length of the string_ref
		                                                ) : element_type_untrimmed{ std::move( arg_rhs.element_type_untrimmed ) },
		                                                    element_type          {
		                                                    	std::next(
		                                                    		element_type_untrimmed.data(),
		                                                    		arg_start
		                                                    	),
		                                                    	arg_length
		                                                    } {
		}

		/// \brief Constructor from lvalue string
		inline element_type_string::element_type_string(const std::string &arg_string ///< The source string
		                                                ) : element_type_untrimmed{ arg_string                                             },
		                                                    element_type          { common::dumb_trim_string_ref( element_type_untrimmed ) } {
		}

		/// \brief Constructor from rvalue string
		inline element_type_string::element_type_string(const std::string &&arg_string ///< The source string
		                                                ) : element_type_untrimmed{ std::move( arg_string )                                },
		                                                    element_type          { common::dumb_trim_string_ref( element_type_untrimmed ) } {
		}

		/// \brief Copy-ctor
		inline element_type_string::element_type_string(const element_type_string &arg_rhs ///< The element_type_string from which to copy-construct
		                                                ) : element_type_untrimmed{ arg_rhs.element_type_untrimmed },
		                                                    element_type          {
		                                                    	std::next(
		                                                    		element_type_untrimmed.data(),
		                                                    		std::distance(
		                                                    			arg_rhs.element_type_untrimmed.data(),
		                                                    			arg_rhs.element_type.data()
		                                                    		)
		                                                    	),
		                                                    	arg_rhs.element_type.length()
		                                                    } {
		}

		/// \brief Move-ctor
		///
		/// A bit trickier because it seems the C++ standard doesn't guarantee that
		/// moving a string won't invalidate iterators/pointers/references into the old one
		inline element_type_string::element_type_string(element_type_string &&arg_rhs ///< The element_type_string from which to move-construct
		                                                ) noexcept : element_type_string{
		                                                              	std::move( arg_rhs ),
		                                                              	std::distance(
		                                                              		arg_rhs.element_type_untrimmed.data(),
		                                                              		arg_rhs.element_type.data()
		                                                              	),
		                                                              	arg_rhs.element_type.length()
		                                                              } {
		}

		/// \brief Copy-assignment operator
		inline element_type_string & element_type_string::operator=(const element_type_string &arg_rhs ///< The element_type_string from which to copy-assign
		                                                            ) {
			element_type_untrimmed = arg_rhs.element_type_untrimmed;
			element_type           = boost::string_ref{
				std::next(
					element_type_untrimmed.data(),
					std::distance(
						arg_rhs.element_type_untrimmed.data(),
						arg_rhs.element_type.data()
					)
				),
				arg_rhs.element_type.length()
			};

			// Return a reference to this element_type_string
			return *this;
		}

		/// \brief Move-assignment operator
		///
		/// A bit trickier because it seems the C++ standard doesn't guarantee that
		/// moving a string won't invalidate iterators/pointers/references into the old one
		inline element_type_string & element_type_string::operator=(element_type_string &&arg_rhs ///< The element_type_string from which to move-assign
		                                                            ) noexcept {
			// Grab the start and length of the string_ref
			const ptrdiff_t start = std::distance(
				arg_rhs.element_type_untrimmed.data(),
				arg_rhs.element_type.data()
			);
			const size_t length = arg_rhs.element_type.length();

			// Move the string
			element_type_untrimmed = std::move( arg_rhs.element_type_untrimmed );

			// Update the string_ref from the newly moved string with the start and length 
			element_type           = boost::string_ref{
				std::next(
					element_type_untrimmed.data(),
					start
				),
				length
			};

			// Return a reference to this element_type_string
			return *this;
		}

		/// \brief Getter for the original trimmed string
		inline const std::string & element_type_string::get_element_type_untrimmed() const {
			return element_type_untrimmed;
		}

		/// \brief Getter for the trimmed element type string (as a string_ref)
		inline const boost::string_ref & element_type_string::get_element_type() const {
			return element_type;
		}

		/// \brief Get the coarse_element_type corresponding to the specified trimmed element string
		inline coarse_element_type get_coarse_element_type(const boost::string_ref &arg_trimmed_element_str ///< The trimmed element string (as it appears in PDB ATOM records)
		                                                   ) {
			if ( arg_trimmed_element_str.empty() || arg_trimmed_element_str.length() > 2 ) {
				return coarse_element_type::NON_CORE;
			}
			else if ( arg_trimmed_element_str.front() == 'C' ) {
				if      ( arg_trimmed_element_str == "CA" ) {
					return coarse_element_type::CARBON_ALPHA;
				}
				else if ( arg_trimmed_element_str == "C"  ) {
					return coarse_element_type::CARBON;
				}
				else if ( arg_trimmed_element_str == "CB" ) {
					return coarse_element_type::CARBON_BETA;
				}
				else {
					return coarse_element_type::NON_CORE;
				}
			}
			else if ( arg_trimmed_element_str == "N" ) {
				return coarse_element_type::NITROGEN;
			}
			else if ( arg_trimmed_element_str == "O" ) {
				return coarse_element_type::OXYGEN;
			}
			else {
				return coarse_element_type::NON_CORE;
			}
		}

		/// \brief Get the coarse_element_type corresponding to the specified element_type_string
		inline coarse_element_type get_coarse_element_type(const element_type_string &arg_element_type_string ///< The element_type_string to query
		                                                   ) {
			return get_coarse_element_type( arg_element_type_string.get_element_type() );
		}

	} // namespace file
} // namespace cath

#endif
