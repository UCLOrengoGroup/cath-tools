package Ssap::Moose;

use Moose ();
use Moose::Exporter;
use Ssap::Types;
use namespace::autoclean;

Moose::Exporter->setup_import_methods(
	also => 'Moose',
);

1;
