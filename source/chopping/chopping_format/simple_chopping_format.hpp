/// \file
/// \brief The simple_chopping_format class header

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

#ifndef _CATH_TOOLS_SOURCE_CHOPPING_CHOPPING_FORMAT_SIMPLE_CHOPPING_FORMAT_H
#define _CATH_TOOLS_SOURCE_CHOPPING_CHOPPING_FORMAT_SIMPLE_CHOPPING_FORMAT_H

#include <boost/utility/string_ref.hpp>

#include "chopping/chopping_format/chopping_format.hpp"

namespace cath {
	namespace chop {

		/// \brief TODOCUMENT
		class simple_chopping_format final : public chopping_format {
		private:
			virtual std::unique_ptr<chopping_format> do_clone() const override final;

			virtual bool do_represents_fragments() const override final;

			virtual domain do_parse_domain(const std::string &) const override final;

		public:
			region parse_segment(const boost::string_ref &) const;

			residue_name parse_residue(const boost::string_ref &) const;
		};

	} // namespace chop
} // namespace cath

#endif
