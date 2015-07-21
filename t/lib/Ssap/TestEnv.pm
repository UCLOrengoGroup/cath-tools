package Ssap::TestEnv;

use Ssap::Moose;
use File::Temp qw/ tempdir /;
use Ssap::Types qw/ AbsoluteFile HashRef Dir /;

BEGIN {
  my $bin_dir = defined $ENV{CATH_BINARIES_BIN_DIR} ? $ENV{CATH_BINARIES_BIN_DIR} : undef;

  if ($bin_dir) {
    warn "Prepending $bin_dir to PATH for testing...\n";
    $ENV{'PATH'} = "$bin_dir:".$ENV{'PATH'};
  }
}

has 'ssap_exe' => (
	is => 'ro',
	isa => AbsoluteFile,
	coerce => 1,
	required => 1,
);

has 'env' => (
	traits => [ 'Hash' ],
	is => 'ro',
	isa => HashRef,
	default => sub { \%ENV },
	handles => {
		has_env => 'exists',
		get_env => 'get',
		set_env => 'set',
	}
);

has 'working_dir' => (
	is => 'ro',
	isa => Dir,
	default => sub { Path::Class::Dir->new( tempdir() ) },
);

1;
