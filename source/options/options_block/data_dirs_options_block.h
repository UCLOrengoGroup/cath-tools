/// \file
/// \brief The data_dirs_options_block class header

/// \copyright
/// CATH Binaries - Protein structure comparison tools such as SSAP and SNAP
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

#ifndef DATA_DIRS_OPTIONS_BLOCK_H_INCLUDED
#define DATA_DIRS_OPTIONS_BLOCK_H_INCLUDED

#include <boost/math/constants/constants.hpp>

#include "options/options_block/data_option.h"
#include "options/options_block/options_block.h"
#include "display/display_spec/display_spec.h"

namespace cath { namespace opts { class pdbs_acquirer; } }
namespace cath { class display_colour_list; }

namespace cath {
	namespace opts {

		/// \brief Represent the different file types that a data_dirs_options_block manages
		enum class data_file {
			PDB,  ///< Enum value to denote PDB  files
			DSSP, ///< Enum value to denote DSSP files
			WOLF, ///< Enum value to denote wolf files
			SEC   ///< Enum value to denote sec  files
		};

		/// \brief TODOCUMENT
		using data_file_vec = std::vector<data_file>;

		/// \brief TODOCUMENT
		using data_file_set = std::set<data_file>;

		/// \brief TODOCUMENT
		using data_file_path_pair = std::pair<data_file, boost::filesystem::path>;

		/// \brief TODOCUMENT
		using data_file_path_map = std::map<data_file, boost::filesystem::path>;

		/// \brief Manage the program options required to find various different types of file from a name
		///        (eg find a DSSP file from the name 1c0pA01)
		///
		/// Note that there are two different uses of the word 'path' around this code:
		///  * the Boost filesystem type that is used to store the location of a file or directory
		///  * a prioritised list of directories to search for some sort of file
		///
		/// For each file type, options are provided to specify
		///  - the prefix, which is prepended to the name when constructing the filename,
		///  - the suffix, which is appended to the name when constructing the filename and
		///  - the path, which is a colon-separated list of directories through which to search for the file (in descending order of preference).
		///
		/// \todo Separate this into two classes: a class that stores data_dirs options and an options_block
		///       that holds a reference to an object of the other class and configures it using program options.
		///
		/// \todo Add a (non-member, non-friend) function to do the reverse operation: ie use a data_dirs_options_block and data_file to deduce a
		///       sensible name from a file (by stripping of the path and then any suitable leading prefix or trailing suffix that
		///       don't leave an empty name). For example, this would generate the name 1c0pA01 from a dssp filename "/cath/data/current/dssp/1c0pA01.dssp").
		class data_dirs_options_block final : public cath::opts::options_block {
		private:
			using data_file_str_map = std::map<data_file, std::string>;
			using data_file_str_pair = data_file_str_map::value_type;

			using data_option_str_map = std::map<detail::data_option, std::string>;
			using data_option_str_pair = data_option_str_map::value_type;

			using file_option_str_map_map = std::map<data_file, data_option_str_map>;

			static const file_option_str_map_map DATA_FILE_TYPE_OPTION_DEFAULTS;
			static const data_file_str_map       DATA_FILE_NAMES;
			static const data_option_str_map     DATA_OPTION_SUFFIXES;
			static const data_option_str_map     DATA_OPTION_DESCRIPTION_START;
			static const data_option_str_map     DATA_OPTION_DESCRIPTION_END;

			/// \brief Store the actual values by data_file and then data_option
			file_option_str_map_map values;

			/// \brief A root directory from which the other directories can be constructed
			boost::filesystem::path cath_root_dir;

			virtual std::unique_ptr<options_block> do_clone() const override final;
			virtual std::string do_get_block_name() const override final;
			virtual void do_add_visible_options_to_description(boost::program_options::options_description &) override final;
			virtual void do_add_hidden_options_to_description(boost::program_options::options_description &) override final;
			virtual opt_str do_invalid_string() const override final;

			void set_value_of_option_and_data_file(const detail::data_option &,
			                                       const data_file &,
			                                       const std::string &);

			std::string get_value_of_option_and_data_file(const detail::data_option &,
			                                              const data_file &) const;

		public:
			data_dirs_options_block();
			virtual ~data_dirs_options_block() noexcept = default;

			void set_path_of_data_file(const data_file &,
			                           const std::string &);

			static std::string get_name_of_data_file(const data_file &);

			boost::filesystem::path get_cath_root_dir() const;
			std::string get_path_of_data_file(const data_file &) const;
			std::string get_prefix_of_data_file(const data_file &) const;
			std::string get_suffix_of_data_file(const data_file &) const;

			static const data_file_str_map DEFAULT_SUBDIR_NAME;
			static const std::string       PO_CATH_ROOT_DIR;
		};

		data_dirs_options_block build_data_dirs_options_block_of_path(const path_vec &);

		data_dirs_options_block build_data_dirs_options_block_of_dir(const boost::filesystem::path &);

		path_vec get_path_of_data_file(const data_dirs_options_block &,
		                                                           const data_file &);

		boost::filesystem::path find_file(const data_dirs_options_block &,
		                                  const data_file &,
		                                  const std::string &);

		boost::filesystem::path find_file(const path_vec &,
		                                  const std::string &);

		path_vec split_path_into_directories(const std::string &);

		std::string join_directories_into_path(const path_vec &);

		std::ostream & operator<<(std::ostream &,
		                          const data_file &);

		size_t str_length_of_data_file(const data_file &);

		size_t max_data_file_str_length();
	}
}

#endif