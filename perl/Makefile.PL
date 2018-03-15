use strict;
use warnings;
use inc::Module::Install;
#use Module::Install::External;
use 5.008;

name 'cath-tools';
perl_version '5.008';
# all_from 'lib/Cath.pm';
license 'perl';
version '0.0.1';

test_requires 'File::Slurp';
test_requires 'File::Temp';
test_requires 'IPC::Cmd';
test_requires 'Path::Class';
test_requires 'Readonly';
test_requires 'Test::Class';
test_requires 'Test::Files';
test_requires 'Try::Tiny';

tests_recursive();

WriteAll();
