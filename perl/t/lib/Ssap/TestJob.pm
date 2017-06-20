package Ssap::TestJob;

use Ssap::Moose;
use Ssap::Types qw/ File Dir /;
use MooseX::Params::Validate;
use JSON::Any;
use FindBin;
use Benchmark::Timer;
use Ssap::TestJob::Result;
use namespace::autoclean;


has 'debug_flag' => ( is => 'ro', isa => 'Bool', default => sub { 0 } );

has 'json_parser' => (
	is => 'ro',
	default => sub { JSON::Any->new() },
);

has 'bootstrap_dir' => (
	is => 'ro',
	isa => Dir,
	lazy_build => 1,
);

sub _build_bootstrap_dir {
	my $self = shift;
	my $package_name = __PACKAGE__;
	$package_name =~ s{::}{__}xms;
	my $bootstrap_dir = Path::Class::Dir->new( $FindBin::Bin, 'data', $package_name );
	if ( !-d $bootstrap_dir ) {
		$bootstrap_dir->mkpath;
	}
	return $bootstrap_dir;
}

has 'default_options' => (
	is => 'ro',
	isa => 'Str',
	default => '',
);

has 'id1' => (
	is => 'rw',
	isa => 'Str',
	writer => 'set_id1',
);

has 'id2' => (
	is => 'rw',
	isa => 'Str',
	writer => 'set_id2',
);

sub get_bootstrap_data {
	my ($self, $key) = @_;
	
	my $result_file = $self->bootstrap_dir->file( $key . '.json' );
	my $result_text = $result_file->slurp();

	my $result_data = $self->json_parser->from_json( $result_text );

	return $result_data;
}

sub set_bootstrap_data {
	my ($self, $key, $data) = @_;

	my $result_file = $self->bootstrap_dir->file( $key, '.json' );
	my $result_text = $self->json_parser->to_json( $data );
	
	$result_file->spew( $result_text );
	
	return 1;
}

sub run_ssap {
	my ($self, %params) = validated_hash( \@_, 
		bootstrap => { isa => 'Bool', default => 0 },
		ssap_exe  => { isa => File, coerce => 1 },
	);

	my $ssap_exe = $params{ssap_exe};
	
	my $ssap_com = join( " ", $ssap_exe, $self->default_options, $self->id1, $self->id2 );
	my $t = Benchmark::Timer->new;
	my $tag = 'ssap_timer';

	if ( $self->debug_flag ) {
		print STDERR "# Running `$ssap_com` ... ";
	}

	$t->start($tag);
	my $ssap_out = `$ssap_com`;
	$t->stop($tag);
	
	my $t_diff = $t->result( $tag );

	if ( $self->debug_flag ) {
		print STDERR " [$t_diff]\n";
	}

	confess "! Error: failed to get STDOUT while running command: `$ssap_com`"
		unless $ssap_out;

	my $ssap_result = Ssap::Result->new_from_string( $ssap_out );

	return Ssap::TestJob::Result->new(
		ssap_result => $ssap_result,
		duration    => $t_diff,
		command     => $ssap_com,
		stdout 	    => $ssap_out,
		bootstrap   => $params{bootstrap},
	);
}

sub swap_ids {
	my $self = shift;
	my $id1 = $self->id1;
	my $id2 = $self->id2;
	$self->set_id1( $id2 );
	$self->set_id2( $id1 );
}

1;

