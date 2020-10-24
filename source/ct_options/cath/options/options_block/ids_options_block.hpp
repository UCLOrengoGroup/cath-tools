/// \file
/// \brief The ids_options_block class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK_IDS_OPTIONS_BLOCK_HPP
#define _CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK_IDS_OPTIONS_BLOCK_HPP

#include <boost/filesystem/path.hpp>

#include "cath/common/type_aliases.hpp"
#include "cath/options/options_block/options_block.hpp"

#include <fstream>
#include <string>

namespace cath {
	namespace opts {

		/// \brief TODOCUMENT
		class ids_options_block final : public options_block {
		private:
			/// \brief TODOCUMENT
			using super = options_block;

			/// \brief TODOCUMENT
			str_vec ids;

			std::unique_ptr<options_block> do_clone() const final;
			std::string do_get_block_name() const final;
			void do_add_visible_options_to_description(boost::program_options::options_description &,
			                                           const size_t &) final;
	//		void do_add_hidden_options_to_description(boost::program_options::options_description &,
	//		                                          const size_t &) final;
			str_opt do_invalid_string(const boost::program_options::variables_map &) const final;
			str_vec do_get_all_options_names() const final;

		public:
			const str_vec & get_ids() const;
			
			static const std::string PO_ID;
		};
		
		// str_vec get_all_ids(const ids_options_block &);
		bool ids_specified(const ids_options_block &);
		std::string get_id_a(const ids_options_block &);
		std::string get_id_b(const ids_options_block &);
	} // namespace opts
} // namespace cath

#endif // _CATH_TOOLS_SOURCE_CT_OPTIONS_CATH_OPTIONS_OPTIONS_BLOCK_IDS_OPTIONS_BLOCK_HPP
