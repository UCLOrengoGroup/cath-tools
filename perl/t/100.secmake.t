use Test::More;

use strict;
use warnings;

use File::Temp;
use FindBin;
use IPC::Cmd qw/ run can_run /;
use Path::Class;
use Readonly;
use Test::Files;

ok(1);

# Readonly my $ids_to_test_filename => $FindBin::Bin.'/full_id_list.txt';
# my @ids_to_test = file($ids_to_test_filename)->slurp();
# chomp(@ids_to_test);
# is(scalar(@ids_to_test), 59, 'There are the correct number of IDs to test');
# 
# Readonly my $data_dir  => $FindBin::Bin.'/data';
# Readonly my $pdb_dir   => "$data_dir/pdb";
# Readonly my $wolf_dir  => "$data_dir/wolf";
# Readonly my $sec_dir   => "$data_dir/sec";
# Readonly my $grath_dir => "$data_dir/grath";
# 
# Readonly my $secmake_exe => $FindBin::Bin.'/../build_release/secmake_src/secmake';
# ok(can_run($secmake_exe), 'Can run the secmake executable');
# 
# is(1, 0, 'Should not have to pass --usewolf flag as a default');
# 
# foreach my $id_to_test (@ids_to_test) {
# 	my $tmp_grath_output = File::Temp->new( SUFFIX => '.gth' );
# 	my $tmp_sec_output   = File::Temp->new( SUFFIX => '.sec' );
# 	my $secmake_command = [
# 		$secmake_exe,
# 		'--pdb',   "$pdb_dir/$id_to_test",
# 		'--wolf',  "$wolf_dir/$id_to_test.wolf",
# 		'--sec',   $tmp_sec_output,
# 		'--grath', $tmp_grath_output,
# 		'--id',    $id_to_test,
# 		'--usewolf',
# 	];
# 	my ($ok, $err, $full_buf, $stdout_buff, $stderr_buff) = run(command => $secmake_command);
# 	is (join("\n", @$stdout_buff), '', 'Empty stdout');
# 	is (join("\n", @$stderr_buff), '', 'Empty stderr');
# 	is($ok, 1, 'secmake runs OK');
# 	compare_ok($tmp_sec_output,   "$sec_dir/$id_to_test.sec",   'Sec   files compare OK' );
# 	compare_ok($tmp_grath_output, "$grath_dir/$id_to_test.gth", 'Grath files compare OK' );
# }

done_testing();