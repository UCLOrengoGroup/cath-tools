use Test::More tests => 2 + 34 * 34 * 4;

use strict;
use warnings;

use Carp qw/ cluck confess /;
use Cwd;
# use Data::Dumper; # *** TEMPORARY ***
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

Readonly my $OVERWRITING_TESTS_WITH_NEW_RESULTS => 0;

Readonly my $ids_to_test_filename => $FindBin::Bin.'/long_pairwise_domains.txt';
my @ids_to_test = file($ids_to_test_filename)->slurp();
chomp(@ids_to_test);
foreach my $id_to_test (@ids_to_test) {
	$id_to_test =~ s/\s.*//g;
}
@ids_to_test = sort(@ids_to_test);
my $num_ids_to_test = scalar(@ids_to_test);
is($num_ids_to_test, 34, 'There are the correct number of IDs to test');

Readonly my $data_dir  => $FindBin::Bin.'/../t/data';
Readonly my $pdb_dir   => "$data_dir/pdb";
Readonly my $wolf_dir  => "$data_dir/wolf";
Readonly my $sec_dir   => "$data_dir/sec";

# Readonly my $ssap_exe => '/usr/local/svn/trunk/bin/SSAP';
# Readonly my $ssap_exe => '/usr/local/svn/trunk/bin/cathedral_pairwise';
Readonly my $ssap_exe => $FindBin::Bin.'/../release/SSAP';
# Readonly my $ssap_exe => $FindBin::Bin.'/../debug/SSAP';

ok(can_run($ssap_exe), 'Can run the SSAP executable');

# 1krmA00 1jl9B00 does	 not produce an alignment
# 3ncyW02 1b3qB01 does not produce an alignment

$ENV{DOMDIR}  = $pdb_dir;
$ENV{WOLFDIR} = $wolf_dir;
$ENV{SECDIR}  = $sec_dir;
my $orig_dir = cwd;

for (my $id_ctr1 = 0; $id_ctr1 < $num_ids_to_test; ++$id_ctr1) {
	my $id1 = $ids_to_test[$id_ctr1];
	for (my $id_ctr2 = 0; $id_ctr2 < $num_ids_to_test; ++$id_ctr2) {
		my $id2 = $ids_to_test[$id_ctr2];
		
		my $tempdir = File::Temp->newdir(CLEANUP => 0);
		chdir $tempdir;
		try {
			my $ssap_command = [
				$ssap_exe,
				'--xmlsup',
				$id1,
				$id2
			];
			warn "About to run \"".join(" ", @$ssap_command)."\" in dir ".cwd()."\n";
			my ($ok, $err, $full_buf, $stdout_buff, $stderr_buff) = run(command => $ssap_command);
# 			die Dumper([$stdout_buff, $stderr_buff]);
			is(scalar(@$stdout_buff), 1, "$id1 vs $id2 : SSAP produced one line of stdout output");
			my $ssap_result_line = $stdout_buff->[0];
			chomp($ssap_result_line);
			my @ssap_results_parts = split(/\s+/, $ssap_result_line);
			
			my $expected_results_file           = "$FindBin::Bin/data/ssap_results/$id1$id2.result_line.txt";
			my $got_alignment_file              = "$id1$id2.list";
			my $expected_alignment_file         = "$FindBin::Bin/data/ssap_results/$got_alignment_file";
			my $got_xml_superposition_file      = "$id1$id2.superpose.xml";
			my $expected_xml_superposition_file = "$FindBin::Bin/data/ssap_results/$got_xml_superposition_file";
			
			is(        "$ssap_result_line\n",       file($expected_results_file)->slurp(),  "$id1 vs $id2 : SSAP result line matches previous");
			
			SKIP: {
				skip "Pair $id1 vs $id2 fail to output alignment", 2 if ("$id1:$id2" eq "1krmA00:1jl9B00" || "$id1:$id2" eq "3ncyW02:1b3qB01");
				compare_ok($got_alignment_file,         $expected_alignment_file,               "$id1 vs $id2 : SSAP alignment compares OK");
				compare_ok($got_xml_superposition_file, $expected_xml_superposition_file,       "$id1 vs $id2 : XML superposition file compares OK");
				
				if ($OVERWRITING_TESTS_WITH_NEW_RESULTS) {
	# 				if (!-e $expected_results_file) {
						write_file($expected_results_file, "$ssap_result_line\n")
							or cluck "Unable to write results line file : $OS_ERROR";
	# 				}
					
	# 				if (!-e $expected_alignment_file) {
					copy ($got_alignment_file, $expected_alignment_file)
						or cluck "Unable to copy got file \"$got_alignment_file\" to expected file \"$expected_alignment_file\" : $OS_ERROR";
	# 				}
					
	# 				if (!-e $expected_xml_superposition_file) {
					copy ($got_xml_superposition_file, $expected_xml_superposition_file)
						or cluck "Unable to copy got file \"$got_xml_superposition_file\" to expected file \"$expected_xml_superposition_file\" : $OS_ERROR";
	# 				}
				}
			}
		}
		catch {
			die $ARG;
		}
		finally {
			chdir $orig_dir;
		};
	}
}

