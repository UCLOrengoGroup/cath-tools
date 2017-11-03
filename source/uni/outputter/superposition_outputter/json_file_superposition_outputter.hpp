/// \file
/// \brief The json_file_superposition_outputter class header

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

#ifndef _CATH_TOOLS_SOURCE_OUTPUTTER_SUPERPOSITION_OUTPUTTER_JSON_FILE_SUPERPOSITION_OUTPUTTER_H
#define _CATH_TOOLS_SOURCE_OUTPUTTER_SUPERPOSITION_OUTPUTTER_JSON_FILE_SUPERPOSITION_OUTPUTTER_H

#include <boost/filesystem/path.hpp>

#include "common/json_style.hpp"
#include "display/options/display_spec.hpp"
#include "outputter/superposition_outputter/superposition_outputter.hpp"

namespace cath {
	namespace opts {

		/// \brief A superposition_outputter to write
		class json_file_superposition_outputter final : public superposition_outputter {
		private:
			/// \brief The file to which the JSON representing the superposition should be written
			boost::filesystem::path output_file;

			/// \brief The style in which the JSON should be written
			const common::json_style the_json_style = DEFAULT_JSON_STYLE;

			std::unique_ptr<superposition_outputter> do_clone() const final;
			void do_output_superposition(const sup::superposition_context &,
			                             std::ostream &,
			                             const boost::string_ref &) const final;
			bool do_involves_display_spec() const final;

		public:
			explicit json_file_superposition_outputter(const boost::filesystem::path &,
			                                           const common::json_style & = DEFAULT_JSON_STYLE);

			/// \brief The default style to use for outputting the JSON if it isn't specified
			static constexpr common::json_style DEFAULT_JSON_STYLE = common::json_style::PRETTY;
		};

	} // namespace opts
} // namespace cath

#endif
