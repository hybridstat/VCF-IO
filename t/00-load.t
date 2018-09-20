#!perl -T
use 5.10.0;
use strict;
use warnings;
use Test::More;
use Test::File;

plan tests => 4;

BEGIN {
    use_ok( 'VCF::IO' ) || print "Bail out!\n";
    use_ok( 'VCF::Header' ) || print "Bail out!\n";
    use_ok( 'VCF::Record' ) || print "Bail out!\n";
    use_ok( 'VCF::_Helper' ) || print "Bail out!\n";
}

diag( "Testing VCF::IO $VCF::IO::VERSION, Perl $], $^X" );
