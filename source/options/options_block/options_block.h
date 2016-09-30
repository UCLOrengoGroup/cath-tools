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

#ifndef OPTIONS_BLOCK_H_INCLUDED
#define OPTIONS_BLOCK_H_INCLUDED

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "common/type_aliases.h"

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
			virtual std::unique_ptr<options_block> do_clone() const = 0;

			/// \brief Provide a name for the options block, as used in the options description text
			///
			/// This is a pure virtual function (so must be overridden by any concrete, derived classes).
			virtual std::string do_get_block_name() const = 0;

			/// \brief Add the block's visible options to the provided options_description
			///
			/// This is a pure virtual function (so must be overridden by any concrete, derived classes).
			virtual void do_add_visible_options_to_description(boost::program_options::options_description &) = 0;

			/// \brief Add the block's hidden options to the provided options_description
			///
			/// This is a virtual function (which may or may not be overridden by concrete, derived classes).
			virtual void do_add_hidden_options_to_description(boost::program_options::options_description &);

			/// \brief Identify any conflicts that make the currently stored options invalid
			///
			/// This is a pure virtual function (so must be overridden by any concrete, derived classes).
			///
			/// These methods should only catch absolute conflicts that could never be acceptable.
			/// In particular, they should not reject inadequately specifying options for normal work
			/// if that might be reasonable when the user requests help.
			///
			/// \returns A string describing the conflict in the options or an empty string if there's none
			virtual opt_str do_invalid_string() const = 0;

		public:
			std::unique_ptr<options_block> clone() const;
			virtual ~options_block() noexcept = default;

			boost::program_options::options_description get_all_options_description(const size_t &);
			boost::program_options::options_description get_visible_options_description(const size_t &);
			boost::program_options::options_description get_hidden_options_description();

			opt_str invalid_string() const;

			static bool is_acceptable_output_file(const boost::filesystem::path &);
			static bool is_acceptable_input_file(const boost::filesystem::path &,
			                                     const bool &arg_allow_empty = false);
			static bool is_acceptable_output_dir(const boost::filesystem::path &);
			static bool is_acceptable_input_dir(const boost::filesystem::path &);
			static bool is_acceptable_executable(const boost::filesystem::path &);

			static const std::string SUB_DESC_SEPARATOR;
			static const std::string SUB_DESC_PAIR_SEPARATOR;
		};

		/// \brief Function to make options_block meet the Clonable concept (used in ptr_container)
		///
		/// NOTE: Don't call this yourself. Call the object's clone() method instead because that returns a
		///       smart pointer rather than the raw pointer this has to return to meet the Clonable concept.
		///
		/// This gets the smart pointer from the clone() method and then calls release on it.
		///
		/// \returns A raw pointer to a new copy of the options_block argument, with the same dynamic type.
		///          The caller is responsible for deleting this new object.
		inline options_block * new_clone(const options_block &arg_options_block ///< The options_block to clone
										 ) {
			return arg_options_block.clone().release();
		}
	}
}

#endif
