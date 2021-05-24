/// \file
/// \brief The old_ssap_options_block class header

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

#ifndef _CATH_TOOLS_SOURCE_CT_UNI_CATH_SSAP_OPTIONS_OLD_SSAP_OPTIONS_BLOCK_HPP
#define _CATH_TOOLS_SOURCE_CT_UNI_CATH_SSAP_OPTIONS_OLD_SSAP_OPTIONS_BLOCK_HPP

#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>

#include "cath/alignment/common_residue_selection_policy/common_residue_select_min_score_policy.hpp"
#include "cath/common/path_type_aliases.hpp"
#include "cath/common/type_aliases.hpp"
#include "cath/options/options_block/options_block.hpp"
#include "cath/structure/protein/protein_source_file_set/protein_file_combn.hpp"
#include "cath/superposition/io/sup_pdbs_script_policy.hpp"

namespace cath::opts {

	/// \brief Transitional options_block to (largely) replicate the old SSAP options
	///
	/// As a first step to improving the old SSAP option handling, all of the old SSAP options
	/// have been moved into this one single block and then tidied up.
	///
	/// These options are currently clumped together into one options block simply because
	/// they happened to all come from the old SSAP, not because this is a coherent, reusable
	/// group.
	///
	/// \todo Eventually, all of these options should be replaced or moved out into separate
	///       reusable, coherent groups of options.
	class old_ssap_options_block final : public options_block {
	private:
		using super = options_block;

		static constexpr bool                        DEF_BOOL     { false                                       }; ///< (Default) default for bool_switch options
		static constexpr sup::sup_pdbs_script_policy DEF_SCRIPT   { sup::sup_pdbs_script_policy::LEAVE_RAW_PDBS }; ///< (Default) default for whether to write a short Rasmol script with superposed PDBs
		static constexpr protein_file_combn          DEF_PROT_SRCS{ protein_file_combn::PDB_AND_CALC            }; ///< (Default) default for protein_source_files options
		static constexpr double                      DEF_REFAST   { 65.0                                        }; ///< Default maximum fast SSAP score to trigger running a second fast SSAP with looser cutoffs
		static constexpr double                      DEF_RESLOW   { 75.0                                        }; ///< Default maximum (best) fast SSAP score to trigger running a slow SSAP
		static constexpr double                      DEF_FILE_SC  { 0.0                                         }; ///< Default minimum final SSAP score for outputting alignment/superposition files
		static constexpr double                      DEF_SUP      {
			align::common_residue_select_min_score_policy::MIN_CUTOFF
		}; 

		str_vec                     names;                                        ///< The names of the structures to compare

		bool                        debug                        = DEF_BOOL;      ///< Whether debug level run information has been requested
		::std::filesystem::path     output_filename;                              ///< A file to which all output should be written, or empty if stdout should be used

		::std::filesystem::path     clique_file;                                  ///< A file from which to read clique information
		::std::filesystem::path     domin_file;                                   ///< A file from which to read domin information

		double                      max_score_to_fast_ssap_rerun = DEF_REFAST;    ///< Maximum fast SSAP score to trigger running a second fast SSAP with looser cutoffs
		double                      max_score_to_slow_ssap_rerun = DEF_RESLOW;    ///< Maximum (best) fast SSAP score to trigger running a slow SSAP
		bool                        slow_ssap_only               = DEF_BOOL;      ///< Whether to only run a slow SSAP (and skip all fast SSAPs)

		bool                        use_local_ssap_score         = DEF_BOOL;      ///< Use local score normalised over smallest protein
		bool                        write_all_scores             = DEF_BOOL;      ///< Whether to output all SSAP scores, rather than just the best
		protein_file_combn          protein_source_files         = DEF_PROT_SRCS; ///< The files from which to read the protein (WOLF_SEC or PDB_DSSP_SEC)

		::std::filesystem::path     superposition_dir;                            ///< A directory to which a superposition should be written, or empty if none should be written
		::std::filesystem::path     alignment_dir                = ".";           ///< A directory to which the alignment file should be written
		double                      min_score_for_writing_files  = DEF_FILE_SC;   ///< Minimum final SSAP score for outputting alignment/superposition files
		double                      min_score_for_superposition  = DEF_SUP;       ///< Minimum residue-pair score for inclusion in superposition calculation
		sup::sup_pdbs_script_policy write_rasmol_script          = DEF_SCRIPT;    ///< Whether to write a Rasmol superposition script file
		bool                        write_xml_sup                = DEF_BOOL;      ///< Whether to write an XML superposition file

		[[nodiscard]] std::unique_ptr<options_block> do_clone() const final;
		[[nodiscard]] std::string                    do_get_block_name() const final;
		void do_add_visible_options_to_description(boost::program_options::options_description &,
		                                           const size_t &) final;
		void do_add_hidden_options_to_description(boost::program_options::options_description &,
		                                          const size_t &) final;
		[[nodiscard]] str_opt do_invalid_string( const boost::program_options::variables_map & ) const final;
		[[nodiscard]] str_view_vec do_get_all_options_names() const final;

	  public:
		[[nodiscard]] bool        protein_names_specified() const;
		[[nodiscard]] std::string get_protein_name_a() const;
		[[nodiscard]] std::string get_protein_name_b() const;

		[[nodiscard]] bool get_debug() const;

		[[nodiscard]] bool                    get_output_to_file() const;
		[[nodiscard]] ::std::filesystem::path get_output_filename() const;

		[[nodiscard]] path_opt get_opt_clique_file() const;
		[[nodiscard]] path_opt get_opt_domin_file() const;

		[[nodiscard]] double get_max_score_to_fast_ssap_rerun() const;
		[[nodiscard]] double get_max_score_to_slow_ssap_rerun() const;
		[[nodiscard]] bool   get_slow_ssap_only() const;

		[[nodiscard]] bool                                           get_use_local_ssap_score() const;
		[[nodiscard]] bool                                           get_write_all_scores() const;
		[[nodiscard]] std::unique_ptr<const protein_source_file_set> get_protein_source_files() const;

		[[nodiscard]] path_opt                    get_opt_superposition_dir() const;
		[[nodiscard]] ::std::filesystem::path     get_alignment_dir() const;
		[[nodiscard]] double                      get_min_score_for_writing_files() const;
		[[nodiscard]] double                      get_min_score_for_superposition() const;
		[[nodiscard]] sup::sup_pdbs_script_policy get_write_rasmol_script() const;
		[[nodiscard]] bool                        get_write_xml_sup() const;

		old_ssap_options_block & set_write_rasmol_script(const sup::sup_pdbs_script_policy &);

		static constexpr char PO_CHAR_OUT_FILE = 'o';

		// clang-format off
		static constexpr ::std::string_view PO_NAME                 { "name"                    }; ///< The option name for the names option

		static constexpr ::std::string_view PO_DEBUG                { "debug"                   }; ///< The option name for the debug option
		static constexpr ::std::string_view PO_OUT_FILE             { "outfile"                 }; ///< The option name for the output_filename option

		static constexpr ::std::string_view PO_CLIQUE_FILE          { "clique-file"             }; ///< The option name for the clique_file option
		static constexpr ::std::string_view PO_DOMIN_FILE           { "domin-file"              }; ///< The option name for the domin_file option

		static constexpr ::std::string_view PO_MAX_SCORE_TO_REFAST  { "max-score-to-fast-rerun" }; ///< The option name for the max_score_to_fast_ssap_rerun option
		static constexpr ::std::string_view PO_MAX_SCORE_TO_RESLOW  { "max-score-to-slow-rerun" }; ///< The option name for the max_score_to_slow_ssap_rerun option
		static constexpr ::std::string_view PO_SLOW_SSAP_ONLY       { "slow-ssap-only"          }; ///< The option name for the slow_ssap_only option

		static constexpr ::std::string_view PO_LOC_SSAP_SCORE       { "local-ssap-score"        }; ///< The option name for the use_local_ssap_score option
		static constexpr ::std::string_view PO_ALL_SCORES           { "all-scores"              }; ///< The option name for the write_all_scores option
		static constexpr ::std::string_view PO_PROTEIN_SOURCE_FILES { "prot-src-files"          }; ///< The option name for the protein_source_files option

		static constexpr ::std::string_view PO_SUPN_DIR             { "supdir"                  }; ///< The option name for the superposition_dir option
		static constexpr ::std::string_view PO_ALIGN_DIR            { "aligndir"                }; ///< The option name for the alignment_dir option
		static constexpr ::std::string_view PO_MIN_OUT_SCORE        { "min-score-for-files"     }; ///< The option name for the min_score_for_writing_files option
		static constexpr ::std::string_view PO_MIN_SUP_SCORE        { "min-sup-score"           }; ///< The option name for the min_score_for_superposition option
		static constexpr ::std::string_view PO_RASMOL_SCRIPT        { "rasmol-script"           }; ///< The option name for the write_rasmol_script option
		static constexpr ::std::string_view PO_XML_SUP              { "xmlsup"                  }; ///< The option name for write_xml_sup option
		// clang-format on
	};

	bool has_clique_file(const old_ssap_options_block &);
	::std::filesystem::path get_clique_file(const old_ssap_options_block &);
	bool has_domin_file(const old_ssap_options_block &);
	::std::filesystem::path get_domin_file(const old_ssap_options_block &);
	bool has_superposition_dir(const old_ssap_options_block &);
	::std::filesystem::path get_superposition_dir(const old_ssap_options_block &);

} // namespace cath::opts

#endif // _CATH_TOOLS_SOURCE_CT_UNI_CATH_SSAP_OPTIONS_OLD_SSAP_OPTIONS_BLOCK_HPP
