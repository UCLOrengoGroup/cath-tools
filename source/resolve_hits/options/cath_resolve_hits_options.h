/// \file
/// \brief The cath_resolve_hits_options class header

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

#ifndef CATH_RESOLVE_HITS_OPTIONS_H_INCLUDED
#define CATH_RESOLVE_HITS_OPTIONS_H_INCLUDED

#include "options/executable/executable_options.h"
#include "resolve_hits/options/resolve_hits_input_options_block.h"

#include <iosfwd>

namespace cath {
	namespace rslv {

		/// \brief TODOCUMENT
		class cath_resolve_hits_options final : public opts::executable_options {
		private:
			using super = opts::executable_options;

			static const std::string STANDARD_USAGE_ERROR_STRING;

			/// \brief TODOCUMENT
			resolve_hits_input_options_block the_resolve_hits_input_options_block;

			virtual std::string do_get_program_name() const override final;
			virtual boost::program_options::positional_options_description get_positional_options() override final;
			virtual opt_str do_get_error_or_help_string() const override final;

			virtual std::string do_get_help_prefix_string() const override final;
			virtual std::string do_get_help_suffix_string() const override final;
			virtual std::string do_get_overview_string() const override final;

			// void check_ok_to_use() const;

		public:
			cath_resolve_hits_options();
			virtual ~cath_resolve_hits_options() noexcept = default;

			const resolve_hits_input_spec & get_resolve_hits_input_spec() const;

			static const std::string PROGRAM_NAME;
		};

	}
}

#endif
