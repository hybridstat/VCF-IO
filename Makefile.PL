use 5.10.0;
use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME             => 'VCF::IO',
    AUTHOR           => q{Panagiotis Moulos <pmoulos@hybridstat.gr>},
    VERSION_FROM     => 'lib/VCF/IO.pm',
    ABSTRACT_FROM    => 'lib/VCF/IO.pm',
    LICENSE          => 'mit',
    PL_FILES         => {},
    MIN_PERL_VERSION => '5.006',
    CONFIGURE_REQUIRES => {
        'ExtUtils::MakeMaker' => '0',
    },
    BUILD_REQUIRES => {
        'Test::More' => '0',
        'Test::File' => '0',
        'Test::Exception' => '0'
    },
    PREREQ_PM => {
        'Tie::IxHash::Easy'      => '0',
        'IO::Uncompress::Gunzip' => '0',
        'List::Util' => '0'
    },
    dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean => { FILES => 'VCF-IO-*' },
);
