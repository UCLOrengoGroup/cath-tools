package Ssap::Types;

use base 'MooseX::Types::Combine';

__PACKAGE__->provide_types_from(qw/ Ssap::MyTypes MooseX::Types::Moose MooseX::Types::Path::Class /);

1;
