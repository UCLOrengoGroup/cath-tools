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

# TODO:
# these dependencies could probably be cut down without a huge amount of effort 
# if this causes a significant problem to anyone.
# The following list was based on:
# `grep -R -h -oP 'use\s+([^\s;])' ./t | sort | uniq`
test_requires 'Benchmark::Timer';
test_requires 'File::Copy';
test_requires 'File::Slurp';
test_requires 'IPC::Cmd';
test_requires 'JSON::Any';
test_requires 'Moose';
test_requires 'Moose::Exporter';
test_requires 'MooseX::Params::Validate';
test_requires 'MooseX::Storage';
test_requires 'MooseX::Types';
test_requires 'MooseX::Types::Moose';
test_requires 'MooseX::Types::Path::Class';
test_requires 'namespace::autoclean';
test_requires 'Path::Class';
test_requires 'Readonly';
#test_requires 'Test::Benchmark'; # not installing cleanly
test_requires 'Test::Class::Load';
test_requires 'Test::Files';
test_requires 'Test::More';
test_requires 'Try::Tiny';

tests_recursive();

WriteAll();
