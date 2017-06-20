package Ssap::Result;

use Moose;
use MooseX::Storage;

with Storage( 'format' => 'JSON', 'io' => 'File' );

my $VERSION = '0.01';

has 'id1'                 => ( is => 'ro', isa => 'Str' );
has 'id2'                 => ( is => 'ro', isa => 'Str' );
has 'length1'             => ( is => 'ro', isa => 'Int' );
has 'length2'             => ( is => 'ro', isa => 'Int' );
has 'ssap_score'          => ( is => 'ro', isa => 'Num' );
has 'equivalent_residues' => ( is => 'ro', isa => 'Int' );
has 'percent_overlap'     => ( is => 'ro', isa => 'Str' );
has 'percent_identity'    => ( is => 'ro', isa => 'Str' );
has 'rmsd'                => ( is => 'ro', isa => 'Str' );

sub new_from_string {
	shift if $_[0] eq __PACKAGE__;
	my $ssap_line = shift;

	my ($id1, $id2, $len1, $len2, $ssap, $equiv_res, $per_ov, $per_id, $rmsd) = split(/\s+/, $ssap_line); 
	return __PACKAGE__->new(
		id1                  => $id1,
		id2                  => $id2,
		length1              => $len1,
		length2              => $len2,
		ssap_score           => $ssap,
		equivalent_residues  => $equiv_res,
		percent_overlap      => $per_ov,
		percent_identity     => $per_id,
		rmsd                 => $rmsd,
	);
}

sub swapped_copy {
	my $self = shift;
	
	#warn "swapped_copy.self: " . $self->dump();
	
	return __PACKAGE__->new(
		id1                  => $self->id2,
		id2                  => $self->id1,
		length1              => $self->length2,
		length2              => $self->length1,
		ssap_score           => $self->ssap_score,
		equivalent_residues  => $self->equivalent_residues,
		percent_overlap      => $self->percent_overlap,
		percent_identity     => $self->percent_identity,
		rmsd                 => $self->rmsd,
	);
}

no Moose;

1;
