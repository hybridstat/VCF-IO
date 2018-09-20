#!perl -T
use 5.10.0;
use strict;
use warnings;
use Test::More;
use Test::Exception;

use VCF::IO;
use VCF::Header;
use VCF::Record;
use VCF::_Helper;

plan tests => 1;

BEGIN {
    subtest 'Check VCF file parsing' => sub {
        plan tests => 8;
        
        my $t1 = "./t/t1.vcf"; # SNPs
        my $t2 = "./t/t2.vcf"; # Structural
        
        my $vcf1 = VCF::IO->new({
            'file' => $t1
        });
        my $vcf2 = VCF::IO->new({
            'file' => $t2
        });
        
        dies_ok(sub { $vcf1->validate; },"SNP VCF file just validate fails");
        lives_ok(sub { $vcf1->parse; },"SNP VCF file parse ok");
        lives_ok(sub {
            $vcf1->parse;
            $vcf1->validate;
        },"SNP VCF file parse, validate ok");
        lives_ok(sub { $vcf1->parse_and_validate; },
            "SNP VCF file parse and validate ok");
        
        dies_ok(sub { $vcf2->validate; },
			"Structural VCF file just validate fails");
        lives_ok(sub { $vcf2->parse; },"Structural VCF file parse ok");
        lives_ok(sub {
            $vcf2->parse;
            $vcf2->validate;
        },"Structural VCF file parse, validate ok");
        lives_ok(sub { $vcf2->parse_and_validate; },
            "Structural VCF file parse and validate ok");
    };
}
