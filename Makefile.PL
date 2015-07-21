use strict;
use warnings;
use inc::Module::Install;
use Module::Install::External;
use 5.008;

name 'cath_binaries';
perl_version '5.008';
# all_from 'lib/Cath.pm';
license 'perl';
version '0.0.1';

test_requires 'Test::Class';

tests_recursive();

WriteAll();
