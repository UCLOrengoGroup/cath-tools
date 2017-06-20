package Ssap::TestJob::Result;

use Ssap::Moose;
use Ssap::Types qw/ SsapResult /;
use namespace::autoclean;

has 'ssap_result' => (
	is => 'ro',
	isa => SsapResult,
	required => 1,
	coerce => 1,
);

has 'duration' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'command' => (
	is => 'ro',
	isa => 'Str',
	required => 1,
);

has 'stdout' => (
	is => 'ro',
	isa => 'Str',
);

has 'stderr' => (
	is => 'ro',
	isa => 'Str',
);

has 'bootstrap' => (
	is => 'ro',
	isa => 'Bool',
	required => 1,
);


__PACKAGE__->meta->make_immutable;
1;

