package Ssap::MyTypes;

use MooseX::Types -declare => [
	qw(
		TestJob
		TestEnv
        	TestJobResult
		SsapResult
		ArrayOfTestEnv
		ArrayOfTestJobResults
		AbsoluteFile
	)
];

use MooseX::Types::Path::Class qw/ File Dir /;
use MooseX::Types::Moose qw/ ArrayRef HashRef /;

subtype AbsoluteFile,
	as File;

coerce AbsoluteFile,
	from File,
		via { $_->absolute };

class_type TestJob,       { class => 'Ssap::TestJob' };
class_type TestEnv,       { class => 'Ssap::TestEnv' };
class_type TestJobResult, { class => 'Ssap::TestJob::Result' };
class_type SsapResult,    { class => 'Ssap::Result' };

coerce TestJob,
	from ArrayRef,
		via { require Ssap::TestJob; Ssap::TestJob->new( id1 => $_->[0], id2 => $_->[1] ) },
	from HashRef,
		via { require Ssap::TestJob; Ssap::TestJob->new( $_ ) };

coerce SsapResult,
	from HashRef,
		via { require Ssap::Result; Ssap::Result->new( %$_ ) };	

subtype ArrayOfTestEnv,
	as ArrayRef[TestEnv];

coerce ArrayOfTestEnv,
	from ArrayRef[HashRef],
		via { require Ssap::TestEnv; [ map { Ssap::TestEnv->new($_) } @$_ ] };

subtype ArrayOfTestJobResults,
	as ArrayRef[TestJobResult];

coerce ArrayOfTestJobResults,
	from ArrayRef[HashRef],	
		via { require Ssap::TestJob::Result; [ map { Ssap::TestJob::Result->new($_) } @$_ ] };


1;
