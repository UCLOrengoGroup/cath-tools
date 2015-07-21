/// \file
/// \brief The cath_check_pdb_options class header

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

#ifndef CATH_CHECK_PDB_OPTIONS_H_INCLUDED
#define CATH_CHECK_PDB_OPTIONS_H_INCLUDED

#include <boost/filesystem.hpp>

#include "options/executable/executable_options.h"
#include "options/options_block/check_pdb_options_block.h"
#include "options/options_block/misc_help_version_options_block.h"

#include <vector>

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		class cath_check_pdb_options final : public executable_options {
		private:
			using super = executable_options;

			/// \brief TODOCUMENT
			misc_help_version_options_block the_misc_options_block;

			/// \brief TODOCUMENT
			check_pdb_options_block         the_check_pdb_options_block;

			virtual std::string do_get_program_name() const override final;
			virtual boost::program_options::positional_options_description get_positional_options() override final;
			virtual std::string do_update_error_or_help_string(const boost::program_options::options_description &) const override final;

			static std::string get_help_prefix_string();

		public:
			cath_check_pdb_options();
			virtual ~cath_check_pdb_options() noexcept = default;

			boost::filesystem::path get_pdb_file() const;
			bool get_permit_no_atoms() const;

			static const std::string PROGRAM_NAME;
		};

		void check_pdb_file(const boost::filesystem::path &,
		                    const bool &);
	}
}

#endif
