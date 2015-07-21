package Ssap::TestJobComparison;

use Ssap::Moose;
use Ssap::Types qw/ ArrayRef TestJob ArrayOfTestEnv ArrayOfTestJobResults /;
use namespace::autoclean;

has 'job' => (
	is => 'ro',
	isa => TestJob,
	coerce => 1,
	required => 1,
);

has 'test_envs' => (
	traits => [ 'Array' ],
	is => 'ro',
	isa => ArrayOfTestEnv,
	required => 1,
	coerce => 1,
	handles => {
		list_all_test_envs => 'elements',
	}
);

has results => (
	traits => [ 'Array' ],
	is => 'ro',
	isa => ArrayOfTestJobResults,
	handles => {
		list_all_results => 'elements',
		get_result => 'get',
		add_result => 'push',
	}
);

sub run_tests {
	my $self = shift;

	my $test_job = $self->job;
	
	my @working_dirs;
	foreach my $test_env ($self->list_all_test_envs) {
		# run ssap
		my $result = $test_job->run_ssap( test_env => $test_env );
		push @working_dirs, $test_env->working_dir;
	}
	
	# check that the working dirs have the same files	
}

__PACKAGE__->meta->make_immutable;

1;
