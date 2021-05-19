/// \file
/// \brief The options_block class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK_OPTIONS_BLOCK_HPP
#define _CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK_OPTIONS_BLOCK_HPP

#include <filesystem>
#include <string_view>

#include <boost/program_options.hpp>

#include "cath/common/type_aliases.hpp"

namespace cath {
	namespace opts {

		/// \brief Provide ABC-interface to make it simple to define reusable blocks of options
		///
		/// The options_block class should probably only ever be used from the executable_options
		/// class and any derived classes.
		///
		/// The standard operation of an options_block class is:
		///  * it is constructed,
		///  * its do_add_options_to_description() is called with an options_description
		///    (to which it adds its options, including references to its data members),
		///  * the client code uses the options description to parse some options and
		///  * the options_block is queried for the options it has received.
		///
		/// WARNING: this naturally leads to class that have two responsibilities:
		///  * handling the program options
		///  * storing, validating etc the options
		/// ...which may be OK for simpler options_blocks but they should be split as soon as they begin to show the strain.
		class options_block {
		  private:
			/// \brief A standard do_clone() method to act as a virtual copy-ctor
			///
			/// This is a pure virtual function (so must be overridden by any concrete, derived classes).
			[[nodiscard]] virtual std::unique_ptr<options_block> do_clone() const = 0;

			/// \brief Provide a name for the options block, as used in the options description text
			///
			/// This is a pure virtual function (so must be overridden by any concrete, derived classes).
			[[nodiscard]] virtual std::string do_get_block_name() const = 0;

			/// \brief Add the block's visible options to the provided options_description
			///
			/// This is a pure virtual function (so must be overridden by any concrete, derived classes).
			virtual void do_add_visible_options_to_description( boost::program_options::options_description &,
			                                                    const size_t & ) = 0;

			/// \brief Add the block's hidden options to the provided options_description
			///
			/// This is a virtual function (which may or may not be overridden by concrete, derived classes).
			virtual void do_add_hidden_options_to_description( boost::program_options::options_description &, const size_t & );

			/// \brief Identify any conflicts that make the currently stored options invalid
			///
			/// This is a pure virtual function (so must be overridden by any concrete, derived classes).
			///
			/// These methods should only catch absolute conflicts that could never be acceptable.
			/// In particular, they should not reject inadequately specifying options for normal work
			/// if that might be reasonable when the user requests help.
			///
			/// \returns A string describing the conflict in the options or an empty string if there's none
			[[nodiscard]] virtual str_opt do_invalid_string( const boost::program_options::variables_map & ) const = 0;

			/// \brief Pure virtual method with which each concrete options_block must define its full list of options names
			///
			/// This is a pure virtual function (so must be overridden by any concrete, derived classes).
			[[nodiscard]] virtual str_view_vec do_get_all_options_names() const = 0;

		  public:
			options_block() = default;
			[[nodiscard]] std::unique_ptr<options_block> clone() const;
			virtual ~options_block() noexcept = default;

			options_block(const options_block &) = default;
			options_block(options_block &&) noexcept = default;
			options_block & operator=(const options_block &) = default;
			options_block & operator=(options_block &&) noexcept = default;

			void add_visible_options_to_description(boost::program_options::options_description &,
			                                        const size_t &);
			void add_hidden_options_to_description(boost::program_options::options_description &,
			                                       const size_t &);

			boost::program_options::options_description get_all_options_description(const size_t &);
			boost::program_options::options_description get_visible_options_description(const size_t &);
			boost::program_options::options_description get_hidden_options_description(const size_t &);

			[[nodiscard]] str_opt invalid_string( const boost::program_options::variables_map & ) const;

			[[nodiscard]] str_view_vec get_all_options_names() const;

			static bool is_acceptable_output_file(const ::std::filesystem::path &);
			static bool is_acceptable_input_file(const ::std::filesystem::path &,
			                                     const bool & = false);
			static bool is_acceptable_output_dir(const ::std::filesystem::path &);
			static bool is_acceptable_input_dir(const ::std::filesystem::path &);
			static bool is_acceptable_executable(const ::std::filesystem::path &);

			/// \brief A string to use to separate (valid values and their descriptions) from each other
			static constexpr ::std::string_view SUB_DESC_SEPARATOR = "\n   ";

			/// \brief The string to use to separate valid values' names from their descriptions
			static constexpr ::std::string_view SUB_DESC_PAIR_SEPARATOR = " - ";
		};

		bool specifies_option(const boost::program_options::variables_map &,
		                      const std::string &);

		str_vec specified_options(const boost::program_options::variables_map &,
		                          const str_vec &);

		str_vec specified_options_from_block(const boost::program_options::variables_map &,
		                                     const options_block &);

		bool specifies_any_of_options(const boost::program_options::variables_map &,
		                              const str_vec &);

		bool specifies_options_from_block(const boost::program_options::variables_map &,
		                                  const options_block &);

		/// \brief Function to make options_block meet the Clonable concept (used in ptr_container)
		///
		/// NOTE: Don't call this yourself. Call the object's clone() method instead because that returns a
		///       smart pointer rather than the raw pointer this has to return to meet the Clonable concept.
		///
		/// This gets the smart pointer from the clone() method and then calls release on it.
		///
		/// \returns A raw pointer to a new copy of the options_block argument, with the same dynamic type.
		///          The caller is responsible for deleting this new object.
		inline options_block * new_clone(const options_block &prm_options_block ///< The options_block to clone
		                                 ) {
			return prm_options_block.clone().release();
		}
	} // namespace opts
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK_OPTIONS_BLOCK_HPP
