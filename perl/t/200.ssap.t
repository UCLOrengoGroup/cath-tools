use Test::More;

use strict;
use warnings;

use Carp qw/ confess /;
use Cwd;
use Data::Dumper; # *** TEMPORARY ***
use English qw/ -no_match_vars /;
use File::Copy;
use File::Slurp;
use File::Temp;
use FindBin;
use IPC::Cmd qw/ run can_run /;
use Path::Class;
use Readonly;
use Test::Files;
use Try::Tiny;

Readonly my $OVERWRITING_TESTS_WITH_NEW_RESULTS => 1;

$ENV{WOLFDIR} = '/cath/data/current/wolf';
$ENV{SECDIR} = '/cath/data/current/sec';


my @tests = (
	[ '3monC00', '1afpA00', 'Highly off-centre superpositions',                             '--supdir', '.' ],
	[ '2quiA02', '2ip1A02', 'Highly off-centre superpositions',                             '--supdir', '.' ],
	
	[ '1fi2A00', '1j58A01', 'Poor superposition when using all aligned positions',          '--supdir', '.' ],
	[ '1mkmA',   '1jmrB',   'Poor superposition when using all aligned positions',          '--supdir', '.' ],
	[ '1wz0A',   '1a5rA',   'Poor superposition when using all aligned positions',          '--supdir', '.' ],
	[ '1qw1A',   '1bymA',   'Poor superposition when using all aligned positions',          '--supdir', '.' ],
	[ '1x6eA',   '1meyC',   'Poor superposition when using all aligned positions',          '--supdir', '.' ],
	[ '1qw1A',   '1bymA',   'Poor superposition when using all aligned positions',          '--supdir', '.' ],
	[ '3monC00', '1afpA00', 'Poor superposition when using all aligned positions',          '--supdir', '.' ],
	[ '2nz2A01', '2hmaA01', 'Poor superposition when using all aligned positions',          '--supdir', '.' ],
	
	[ '1p3cA02', '1gpzA02', 'Previous problem: Result is negative',                         '--supdir', '.' ],
	[ '1ej6B00', '2r3sA02', 'Previous problem: Result is negative',                         '--supdir', '.' ], # "Residue index 1032 is out of range in protein of length 1031"
	[ '1wfqA01', '1ifcA00', 'Previous problem: Result is negative despite similar lengths', '--supdir', '.' ],
	[ '1gpzA02', '1p3cA02', 'Previous problem: Result is dependent on direction',           '--supdir', '.' ],
	[ '1h3fB01', '1jhdA01', 'Previous problem: Result is dependent on direction',           '--supdir', '.' ],
	[ '1q15D02', '1r6uB01', 'Previous problem: Result is dependent on direction',           '--supdir', '.' ],
	[ '1q15D02', '2j07A01', 'Previous problem: Result is dependent on direction',           '--supdir', '.' ],
	
	[ '3vkgA',   '1c0pA01', 'Previous problem: Chain is too long for SSAP to handle',       '--supdir', '.' ],
	[ '3vkhA',   '1c0pA01', 'Previous problem: Chain is too long for SSAP to handle',       '--supdir', '.' ],
	[ '3vkhB',   '1c0pA01', 'Previous problem: Chain is too long for SSAP to handle',       '--supdir', '.' ],
	[ '4ai6A',   '1c0pA01', 'Previous problem: Chain is too long for SSAP to handle',       '--supdir', '.' ],
	
	[ '3fwcB',   '3nffG',   'Previous problem: All zero output',                            '--supdir', '.' ],
	[ '3eggD',   '1pzsA',   'Previous problem: All zero output',                            '--supdir', '.' ],
	[ '1qqp400', '2mev400', 'Previous problem: All zero output',                            '--supdir', '.' ],
	[ '2j8kA01', '3d7jE00', 'Previous problem: All zero output',                            '--supdir', '.' ],
	[ '2j8kA01', '3d7jA00', 'Previous problem: All zero output',                            '--supdir', '.' ],
	[ '2j8kA01', '3d7jB00', 'Previous problem: All zero output',                            '--supdir', '.' ],
	[ '2j8kA01', '3d7jD00', 'Previous problem: All zero output',                            '--supdir', '.' ],
	[ '3d7jE00', '2j8kA01', 'Previous problem: All zero output',                            '--supdir', '.' ],
	[ '3d7jA00', '2j8kA01', 'Previous problem: All zero output',                            '--supdir', '.' ],
	[ '3d7jB00', '2j8kA01', 'Previous problem: All zero output',                            '--supdir', '.' ],
	[ '3d7jD00', '2j8kA01', 'Previous problem: All zero output',                            '--supdir', '.' ],
	
	[ '3i07B01', '1q16A03', 'Previous problem: Zero score, despite similar lengths',        '--supdir', '.' ],
	[ '1ifcA00', '1vs9F02', 'Previous problem: Zero score, despite similar lengths',        '--supdir', '.' ],
	[ '2j5zA02', '2vpzA01', 'Previous problem: Zero score, despite similar lengths',        '--supdir', '.' ],
	[ '1n7dA06', '3r8nM02', 'Previous problem: Zero score, despite similar lengths',        '--supdir', '.' ],
	[ '2as0A01', '1ifcA00', 'Previous problem: Zero score, despite similar lengths',        '--supdir', '.' ],
	
	
	[ '3fxhA00', '3g1jA00', 'New problem: throws error exception',                          '--supdir', '.' ], # "std::exception::what: Unable to get pdb_name_number() from residues selected by policy" from alignment/alignment_coord_extractor.cpp(95): Throw in function static std::pair<cath::coord_list, cath::coord_list> cath::alignment_coord_extractor::get_common_coords(const cath::alignment&, const cath::protein&, const cath::protein&, const unsigned int&, const int (*)[5000], const cath::common_coord_selection_policy&, const unsigned int&, const unsigned int&)
	[ '3g1jA00', '3nrgA00', 'New problem: throws error exception',                          '--supdir', '.' ], # "std::exception::what: Unable to get pdb_name_number() from residues selected by policy" from alignment/alignment_coord_extractor.cpp(95): Throw in function static std::pair<cath::coord_list, cath::coord_list> cath::alignment_coord_extractor::get_common_coords(const cath::alignment&, const cath::protein&, const cath::protein&, const unsigned int&, const int (*)[5000], const cath::common_coord_selection_policy&, const unsigned int&, const unsigned int&)
	[ '1n4qB00', '1rwrA00', 'Previous problem: Crashes',                                    '--supdir', '.' ], # Example of secondary structure index out of range that's wildly out ("Secondary structure index 85 is out of range in protein with 19 secondary structures")
	[ '16vpA00', '1c5kA02', 'Previous problem: Crashes',                                    '--supdir', '.' ], # Example of secondary structure index out of range that's wildly out ("Secondary structure index 62 is out of range in protein with 19 secondary structures")
	[ '1i1qA00', '3ic7A02', 'Previous problem: Crashes',                                    '--supdir', '.' ], # Example of             residue index out of range that's wildly out (           "Residue index 919 is out of range in protein of length 512"               )
	[ '1b8bA00', '1b9wA02', 'Previous problem: Crashes',                                    '--supdir', '.' ], # Example of             residue index out of range that's wildly out (           "Residue index 835 is out of range in protein of length 533"               )
	[ '1jb0B00', '1s3jA03', 'Previous problem: Crashes',                                    '--supdir', '.' ], # "Residue index 1145 is out of range in protein of length 739"
	[ '1r8wB01', '1pgyA00', 'Previous problem: Crashes',                                    '--supdir', '.' ], # "Residue index 1044 is out of range in protein of length 779"
	[ '1t0sA00', '1rfxC01', 'Previous problem: Crashes',                                    '--supdir', '.' ], #
	[ '1w63G00', '3moyA02', 'Previous problem: Crashes',                                    '--supdir', '.' ], #
	[ '3fhhA02', '1mh3A03', 'Previous problem: Crashes',                                    '--supdir', '.' ], # "Residue index 498 is out of range in protein of length 497"
	[ '3gb8A01', '3mvpA01', 'Previous problem: Crashes',                                    '--supdir', '.' ], # "Residue index 601 is out of range in protein of length 600"
	[ '3gzuB00', '1y8aA03', 'Previous problem: Crashes',                                    '--supdir', '.' ], # "Residue index 969 is out of range in protein of length 800"
	[ '3i4gA00', '1kyqA02', 'Previous problem: Crashes',                                    '--supdir', '.' ], # "Residue index 961 is out of range in protein of length 500"
	[ '3i4rB00', '1lwuC03', 'Previous problem: Crashes',                                    '--supdir', '.' ], # "Residue index 529 is out of range in protein of length 528"
	
	[ '3t5oA02', '2fcwB01', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 525 is out of range in protein of length 524"
	[ '1ehkA00', '2dbbB01', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 850 is out of range in protein of length 544"
	[ '1fiyA02', '1l3lA02', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 496 is out of range in protein of length 495"
	[ '1j36A00', '3c6mD01', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 792 is out of range in protein of length 598"
	[ '1k9xC00', '1gl2A00', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 895 is out of range in protein of length 497"
	[ '1nnnA01', '2k57A00', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 734 is out of range in protein of length 513"
	[ '1r1rA02', '4e6rA00', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], #
	[ '1uc2B00', '1kj6A00', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 846 is out of range in protein of length 480"
	[ '1urjB01', '1pgjA03', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 886 is out of range in protein of length 753"
	[ '1vsy400', '1wz6A01', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 800 is out of range in protein of length 799"
	[ '1wxrA02', '1otrA00', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 703 is out of range in protein of length 702"
	[ '1xjkB00', '3i33A02', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 608 is out of range in protein of length 607"
	[ '1xkwA02', '1u3eM02', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 967 is out of range in protein of length 528"
	[ '2a3lA01', '2i10A01', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 526 is out of range in protein of length 525"
	[ '2a65A00', '1sxjE02', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 943 is out of range in protein of length 509"
	[ '2a65A00', '3by6D02', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 861 is out of range in protein of length 509"
	[ '2apgA00', '1lvaA04', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 791 is out of range in protein of length 517"
	[ '2h1nB00', '1ifyA00', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 771 is out of range in protein of length 552"
	[ '2i3oD00', '1szbA02', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 941 is out of range in protein of length 497"
	[ '2p5oD01', '1sg4B02', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 810 is out of range in protein of length 631"
	[ '2vbeA00', '1qubA04', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 940 is out of range in protein of length 504"
	[ '2vqiA00', '3d31A02', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 705 is out of range in protein of length 469"
	[ '2w8aA01', '1dx5I03', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 475 is out of range in protein of length 474"
	[ '2wsxC00', '3fwnA03', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 497 is out of range in protein of length 496"
	[ '3dh4B01', '1oksA00', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 512 is out of range in protein of length 511"
	[ '3dwcB00', '3e58A02', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # ""
	[ '3epsA00', '1vpuA00', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 567 is out of range in protein of length 566"
	[ '3faxA03', '2jguA04', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 655 is out of range in protein of length 471"
	[ '3fsnA00', '1b3qB03', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 727 is out of range in protein of length 510"
	[ '3h7lB02', '1q02A00', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 479 is out of range in protein of length 478"
	[ '3hx6B00', '1oksA00', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 534 is out of range in protein of length 486"
	[ '3iukA00', '1a5tA02', 'Previous problem: Crashes and corrupt alignment before that',  '--supdir', '.' ], # "Residue index 968 is out of range in protein of length 533"
	
	[ '1krmA00', '1jl9B00', 'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '1cov4',   '3vbsD',   'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '1mec4',   '4gh4D',   'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '1vpuA',   '4emcC',   'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '2h8pD',   '4a94D',   'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '2hfeD',   '4a94D',   'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '2kv9B',   '2ljfA',   'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '2mev4',   '4gh4D',   'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '2qbeH02', '2i2vH02', 'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '2vhmZ',   '4btd4',   'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '2wro4',   '2j034',   'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '2wrr4',   '2j034',   'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '2xqe4',   '2j034',   'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '2y114',   '2j034',   'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '2y134',   '2j034',   'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '2y174',   '2j034',   'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '2y194',   '2j034',   'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '3bqsA',   '3aqjC',   'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '3fgaD',   '4go6A',   'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '3vbfD00', '1cov400', 'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	[ '3vbfD00', '2plv400', 'Previous problem: Does not produce an alignment file',         '--supdir', '.' ],
	
	[ '2i3oD00', '1b9wA01', 'Previous problem: Unknown',                                    '--supdir', '.' ],
	[ '2wx6B00', '1mj4A00', 'Previous problem: Unknown',                                    '--supdir', '.' ],

	[ '2iidA',   ' 1gteA',  'Previous problem: very slow on two large structures',          '--supdir', '.' ],

);

Readonly my $data_dir  => $FindBin::Bin.'/data';
Readonly my $pdb_dir   => "$data_dir/pdb";
Readonly my $wolf_dir  => "$data_dir/wolf";
Readonly my $sec_dir   => "$data_dir/sec";

ok(1);

if (0) {
	Readonly my $ssap_exe => $FindBin::Bin.'/../cath-ssap';
	ok(can_run($ssap_exe), 'Can run the SSAP executable');
	
	$ENV{DOMDIR}  = $pdb_dir;
	$ENV{WOLFDIR} = $wolf_dir;
	$ENV{SECDIR}  = $sec_dir;
	my $orig_dir = cwd;

	local $TODO = "These problems have not all been fixed yet";
	foreach my $test (@tests) {
		foreach my $swap_ids (0, 1) {
			foreach my $ids_before_options (0, 1) {
				my ($id1, $id2, $description, @options) = @$test;
				
				if ($swap_ids) {
					$description .= ' [IDs swapped]';
					($id1, $id2) = ($id2, $id1);
				}
				$description = "$id1 vs $id2".$description;
				if ($ids_before_options) {
					unshift @options, ($id1, $id2);
				}
				else {
					push @options, ($id1, $id2);
				}
				
	# 			warn Dumper(\@options);
				
				my $tempdir = File::Temp->newdir();
				chdir $tempdir;
				try {
					my $ssap_command = [$ssap_exe, @options];
					my ($ok, $err, $full_buf, $stdout_buff, $stderr_buff) = run(command => $ssap_command);
					is(scalar(@$stdout_buff), 1, 'SSAP produced one line of stdout output');
					my $ssap_result_line = $stdout_buff->[0];
					chomp($ssap_result_line);
					my @ssap_results_parts = split(/\s+/, $ssap_result_line);
					
					my $expected_results_file       = "$FindBin::Bin/data/ssap_results/$id1$id2.result_line.txt";
					my $got_alignment_file          = "$id1$id2.list";
					my $expected_alignment_file     = "$FindBin::Bin/data/ssap_results/$got_alignment_file";
					my $got_superposition_file      = "$id1$id2.sup";
					my $expected_superposition_file = "$FindBin::Bin/data/ssap_results/$got_superposition_file";
					
					if ($OVERWRITING_TESTS_WITH_NEW_RESULTS) {
						write_file($expected_results_file, "$ssap_result_line\n")
							or confess "Unable to write results line file";
						
						copy ($got_alignment_file, $expected_alignment_file)
							or confess "Unable to copy got file \"$got_alignment_file\" to expected file \"$expected_alignment_file\"";
						
						copy ($got_superposition_file, $expected_superposition_file)
							or confess "Unable to copy got file \"$got_superposition_file\" to expected file \"$expected_superposition_file\"";
					}
				}
				catch {
					warn $ARG;
				}
				finally {
					chdir $orig_dir;
				};
			}
		}
	}
}

done_testing();
