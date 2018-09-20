#!perl -T
use 5.10.0;
use strict;
use warnings;
use Test::More;
use Test::Exception;
use Test::File;

use VCF::IO;
use VCF::Header;
use VCF::Record;
use VCF::_Helper;

plan tests => 5;

BEGIN {
    subtest 'Check header parsing' => sub {
        plan tests => 9;
        
        my $t1 = "./t/t1.vcf"; # SNPs
        my $t2 = "./t/t2.vcf"; # Structural
        
        my $h1 = VCF::Header->new({
            'file' => $t1
        });
        my $h2 = VCF::Header->new({
            'file' => $t2
        });
        
        dies_ok(sub { $h1->validate; },"SNP header just validate fails");
        lives_ok(sub { $h1->parse; },"SNP header parse ok");
        lives_ok(sub {
            $h1->parse;
            $h1->validate;
        },"SNP header parse, validate ok");
        lives_ok(sub { $h1->parse_and_validate; },
            "SNP header parse and validate ok");
        
        dies_ok(sub { $h2->validate; },"Structural header just validate fails");
        lives_ok(sub { $h2->parse; },"Structural header parse ok");
        lives_ok(sub {
            $h2->parse;
            $h2->validate;
        },"Structural header parse, validate ok");
        lives_ok(sub { $h2->parse_and_validate; },
            "Structural header parse and validate ok");

        is($h1->parse_sample_names->[0],("Sample"),"Sample name parsing ok");
    };
    
    subtest 'Check header values' => sub {
        plan tests => 8;
        
        my $t1 = "./t/t1.vcf"; # SNPs
        my $t2 = "./t/t2.vcf"; # Structural
        
        my $h1 = VCF::Header->new({
            'file' => $t1
        });
        my $h2 = VCF::Header->new({
            'file' => $t2
        });
        
        $h1->parse_and_validate;
        $h2->parse_and_validate;
        
        my $hh1 = $h1->get_header;
        my $hh2 = $h2->get_header;
        
        is($hh1->{"fileformat"},"VCFv4.2","SNP fileformat ok");
        is($hh1->{"contig"}->[0]->{"length"},"249250621","SNP contig ok");
        is($hh1->{"FORMAT"}->[0]->{"ID"},"GT","SNP FORMAT GT ok");
        is(scalar @{$hh1->{"INFO"}},47,"SNP INFO fields ok");
        
        is($hh2->{"fileformat"},"VCFv4.1","Structural fileformat ok");
        is(scalar @{$hh2->{"columns"}},8,"Structural column fields ok");
        is(scalar @{$hh2->{"ALT"}},5,"Structural ALT fields ok");
        is(scalar @{$hh2->{"INFO"}},25,"Structural INFO fields ok");
    };
    
    subtest 'Check header record add' => sub {
        plan tests => 2;
        
        my $t1 = "./t/t1.vcf"; # SNPs
        my $h1 = VCF::Header->new({
            'file' => $t1
        });
        $h1->parse_and_validate;
        
        lives_ok(sub {
            my %attr = (
                "ID" => "NEWANN",
                "Number" => ".",
                "Type" => "String",
                "Description" => "New annotation"
            );
            my $type = "INFO";
            $h1->add(\%attr,$type);
        },"Successful header record add");
        
        dies_ok(sub {
            my %attr = (
                "ID" => "NEWANN",
                "Number" => ".",
                "Type" => "IllegalType",
                "Description" => "New annotation"
            );
            my $type = "INFO";
            $h1->add(\%attr,$type);
        },"Unsuccessful header record add");
    };
    
    subtest 'Check header record change' => sub {
        plan tests => 4;
        
        my $t1 = "./t/t1.vcf"; # SNPs
        my $h1 = VCF::Header->new({
            'file' => $t1
        });
        $h1->parse_and_validate;
        
        lives_ok(sub {
            my $with = {
             'id' => "DP",
             'key' => "Description",
             'value' => "Total read depth at the genomic locus"
            };
            $h1->change("INFO",$with);
        },"Successful header INFO record change");
        is($h1->get_header->{"INFO"}->[1]->{"Description"},
            "Total read depth at the genomic locus",
            "Header INFO record change confirmed");
        
        lives_ok(sub { $h1->change("fileDate","20180918") },
            "Successful header fileDate record change");
        is($h1->get_header->{"fileDate"},"20180918",
            "Header fileDate record change confirmed");
    };
    
    subtest 'Check header writing' => sub {
        plan tests => 4;
        
        my $t1 = "./t/t1.vcf"; # SNPs
        my $h1 = VCF::Header->new({
            'file' => $t1
        });
        $h1->parse_and_validate;
        
        lives_ok(sub { 
            $h1->output("./t/t1h.tmp");
            $h1->write;
            $h1->done;
        },"Successful header write");
        file_readable_ok('./t/t1h.tmp',"Test file t1h.tmp readable");
        file_min_size_ok('./t/t1h.tmp',5000,"t1h.tmp min size ok");
        file_max_size_ok('./t/t1h.tmp',15000,"t1h.tmp max size ok");
        
        unlink("./t/t1h.tmp")
    };
}

