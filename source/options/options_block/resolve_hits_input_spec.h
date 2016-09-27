/// \file
/// \brief The resolve_hits_input_spec class header

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

#ifndef RESOLVE_HITS_INPUT_SPEC_H_INCLUDED
#define RESOLVE_HITS_INPUT_SPEC_H_INCLUDED

#include <boost/optional.hpp>

#include "common/path_type_aliases.h"
#include "common/type_aliases.h"

namespace cath {
	namespace opts {

		class resolve_hits_input_spec final {
		private:
			/// \brief TODOCUMENT
			opt_path input_file;

			/// \brief TODOCUMENT
			bool read_from_stdin          = DEFAULT_READ_FROM_STDIN;

			/// \brief TODOCUMENT
			bool input_hits_are_grouped    = DEFAULT_INPUT_HITS_ARE_GROUPED;

			/// \brief TODOCUMENT
			bool domtblout                = DEFAULT_DOMTBLOUT;

			/// \brief TODOCUMENT
			bool hmmeraln                 = DEFAULT_HMMERALN;

			/// \brief TODOCUMENT
			bool raw_score_is_evalue      = DEFAULT_RAW_SCORE_IS_EVALUE;

			/// \brief TODOCUMENT
			opt_path domtblout_input_file;

			/// \brief TODOCUMENT
			bool apply_cath_rules         = DEFAULT_APPLY_CATH_RULES;

			// scoring numbers:
			// double longer_preference;
			// double better_preference;

			// shrinking numbers
		public:
			/// \brief TODOCUMENT
			static constexpr bool DEFAULT_READ_FROM_STDIN       = false;

			/// \brief TODOCUMENT
			static constexpr bool DEFAULT_INPUT_HITS_ARE_GROUPED = false;

			/// \brief TODOCUMENT
			static constexpr bool DEFAULT_DOMTBLOUT             = false;

			/// \brief TODOCUMENT
			static constexpr bool DEFAULT_HMMERALN              = false;

			/// \brief TODOCUMENT
			static constexpr bool DEFAULT_RAW_SCORE_IS_EVALUE   = false;

			/// \brief TODOCUMENT
			static constexpr bool DEFAULT_APPLY_CATH_RULES      = false;

			const opt_path & get_input_file() const;
			const bool & get_read_from_stdin() const;
			const bool & get_input_hits_are_grouped() const;
			const bool & get_domtblout() const;
			const bool & get_hmmeraln() const;
			const bool & get_raw_score_is_evalue() const;
			const bool & get_apply_cath_rules() const;

			resolve_hits_input_spec & set_input_file(const boost::filesystem::path &);
			resolve_hits_input_spec & set_read_from_stdin(const bool &);
			resolve_hits_input_spec & set_input_hits_are_grouped(const bool &);
			resolve_hits_input_spec & set_domtblout(const bool &);
			resolve_hits_input_spec & set_hmmeraln(const bool &);
			resolve_hits_input_spec & set_raw_score_is_evalue(const bool &);
			resolve_hits_input_spec & set_apply_cath_rules(const bool &);
		};

		opt_str get_invalid_description(const resolve_hits_input_spec &);

	}
}

#endif
