package TestsFor::Ssap;

use base 'Test::Class';
use strict;
use warnings;

use Test::More;
use Try::Tiny;
use Ssap::TestJobComparison;
use Test::Benchmark;
use Path::Class;

local $Test::Benchmark::VERBOSE = 1;

my $ssap_exe_bootstrap = (map { chomp; file($_)->absolute } `which SSAP`)[0];
my $ssap_exe_test      = (map { chomp; file($_)->absolute } `which cath-ssap`)[0];

diag( "" );
diag( "SSAP (bootstrap): $ssap_exe_bootstrap" );
diag( "SSAP (test):      $ssap_exe_test" );

sub single_ssap : Test(1) {
	my ($self, $report) = @_;
	
	my $test = Ssap::TestJobComparison->new(
		job => { id1 => '1oaiA00', id2 => '3frhA01' },
		test_envs => [
			{ ssap_exe => $ssap_exe_bootstrap },
			{ ssap_exe => $ssap_exe_test },
		],
	);

	$test->run_tests();

	my ($bootstrap_result, $test_result) = $test->list_all_results();
	
	is_deeply( $test_result->ssap_result, $bootstrap_result->ssap_result, 'test algorithm matches bootstrap' );
}

sub check_single_slow_ssap : Test(1) {
	my ($self, $report) = @_;
	
	my $job = Ssap::TestJob->new(
		id1 => '1sxjC03',
		id2 => '3frhA01',
	);
	
	my $expected = $job->run_ssap( ssap_exe => $ssap_exe_bootstrap );
	my $got      = $job->run_ssap( ssap_exe => $ssap_exe_test );

	is_deeply( $got->{result}, $expected->{result}, '1sxjC03 vs 3frhA01 ssap results match' );
}


sub time_many_ssaps : Test(1) {
	my ($self, $report) = @_;

	# 1: short vs short ; homologue
	# 2: short vs long  ; homologue
	# 3: long  vs long  ; homologue
	# 4: short vs short ; non-homologue
	# 5: short vs long  ; non-homologue
	# 6: long  vs long  ; non-homologue
	
	my $job1 = Ssap::TestJob->new(
		id1 => '1oaiA00',
		id2 => '3frhA01',
	);

	my $runs = 20;

	my $run_test = sub {
		my $job = shift;
		$job->run_ssap( ssap_exe => $ssap_exe_test );
	};

	my $run_bootstrap = sub {
		my $job = shift;
		$job->run_ssap( ssap_exe => $ssap_exe_bootstrap );
	};

	TODO: {
		local $TODO = "This is currently taking longer than it should...";
		is_faster( $runs, sub { $run_bootstrap->($job1) }, sub { $run_test->($job1) }, 
			"short vs short (homologues) x $runs" );
	}
}


sub reverse_ssap : Test(9) {
	my ($self, $report) = @_;

	my $job = Ssap::TestJob->new(
			id1 => '1oaiA00', 
			id2 => '3frhA01'
		);
	
	my $expected = $job->run_ssap( ssap_exe => $ssap_exe_test );

	$job->swap_ids;
	my $got = $job->run_ssap( ssap_exe => $ssap_exe_test );
	
	# (Unpack,) swap the IDs and lengths (and then pack again)
	$got->{result} = Ssap::Result->unpack($got->{result})->swapped_copy()->pack();
	
	foreach my $field ( qw / id1 id2 length1 length2 ssap_score equivalent_residues percent_overlap percent_identity rmsd / ) {
		is( $got->{result}->{$field}, $expected->{result}->{$field}, "reversed ids have same $field" );
	}

}

sub check_result {
	my ($self, $result1, $result2) = @_;
	
}

1;


