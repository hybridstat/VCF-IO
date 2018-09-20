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

plan tests => 3;

BEGIN {
    subtest 'Check module initialization with file' => sub {
        plan tests => 6;
        
        my $t1 = "./t/t1.vcf";
        
        my $vcf1 = VCF::IO->new({
            'file' => $t1
        });
        my $h1 = VCF::Header->new({
            'file' => $t1
        });
        my $r1 = VCF::Record->new({
            'file' => $t1
        });
        my $vcf2 = VCF::IO->new($t1);
        my $h2 = VCF::Header->new($t1);
        my $r2 = VCF::Record->new($t1);
        
        isa_ok($vcf1,'VCF::IO',"Init VCF::IO from params ok");
        isa_ok($h1,'VCF::Header',"Init VCF::Header from params ok");
        isa_ok($r1,'VCF::Record',"Init VCF::Record from params ok");
        isa_ok($vcf2,'VCF::IO',"Init VCF::IO from file ok");
        isa_ok($h2,'VCF::Header',"Init VCF::Header from file ok");
        isa_ok($r2,'VCF::Record',"Init VCF::Record from file ok");
    };
    
    subtest 'Check module initialization with data' => sub {
        plan tests => 2;
        
        # Define a simple header
        my $h = {
			'fileformat' => 'VCFv4.2',
			'fileDate' => '20170630',
			'reference' => 'hg19.fasta',
			'contig' => [{
				'ID' => '1',
				'length' => '249250621'
			}],
			'INFO' => [{
				'ID' => 'NS',
				'Number' => '1',
				'Type' => 'Integer',
				'Description' => '"Number of samples with data"'
			},{
				'ID' => 'DP',
				'Number' => '1',
				'Type' => 'Integer',
				'Description' => '"Total read depth at the locus"'
			},{
				'ID' => 'AF',
				'Number' => 'A',
				'Type' => 'Float',
				'Description' => '"Estimated allele frequency in the range (0,1]"'
			}],
			'FORMAT' => [{
				'ID' => 'GT',
				'Number' => '1',
				'Type' => 'String',
				'Description' => '"Genotype"'
			},{
				'ID' => 'DP',
				'Number' => '1',
				'Type' => 'Integer',
				'Description' => '"Read Depth"'
			}],
			'other' => {
				'SnpEffVersion' => '"4.2 (build 2015-12-05), by Pablo Cingolani"'
			},
			'columns' => ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO',
				'FORMAT','Sample']
		};
		
		# Define two simple records
		my $r = [{
			'QUAL' => '371.317',
			'INFO' => {
				'AB' => ['0'],
				'ABP' => ['0'],
				'AC' => ['2'],
				'AF' => ['1']
			},
			'ALT' => ['G'],
			'ID' => ['.'],
			'SAMPLE' => [{
				'name' => 'Sample',
				'order' => 0,
				'GT' => ['1/1'],
				'GQ' => ['102.247'],
				'DP' => ['12'],
				'DPR' => ['12','12']
			}],
			'REF' => 'A',
			'FORMAT' => ['GT','GQ','DP','DPR'],
			'CHROM' => '1',
			'POS' => 69511,
			'FILTER' => ['.']
		},{
			'QUAL' => '261.217',
			'INFO' => {
			'AB' => ['0'],
			'ABP' => ['0'],
			'AC' => ['3'],
			'AF' => ['1']
			},
			'ALT' => ['T'],
			'ID' => ['.'],
			'SAMPLE' => [{
				'name' => 'Sample',
				'order' => 0,
				'GT' => ['1/1'],
				'GQ' => ['100.247'],
				'DP' => ['13'],
				'DPR' => ['13','13']
			}],
			'REF' => 'C',
			'FORMAT' => ['GT','GQ','DP','DPR'],
			'CHROM' => '1',
			'POS' => 12345,
			'FILTER' => ['.']
		}];


        my $ho = VCF::Header->new({
            'data' => $h
        });
        my $ro = VCF::Record->new({
            'data' => $r
        });
        
        isa_ok($ho,'VCF::Header',"Init VCF::Header from data ok");
        isa_ok($ro,'VCF::Record',"Init VCF::Header from data ok");
    };
    
    subtest 'Check module initialization fail with faulty data' => sub {
        plan tests => 2;
        
        # Define a simple wrong header
        my $h = {
			'fileDate' => '20170630',
			'reference' => 'hg19.fasta',
			'contig' => [{
				'ID' => '1',
				'length' => '249250621'
			}],
			'INFO' => [{
				'ID' => 'NS',
				'Number' => '1',
				'Type' => 'Integer',
				'Description' => '"Number of samples with data"'
			},{
				'ID' => 'DP',
				'Number' => '1',
				'Type' => 'Integer',
				'Description' => '"Total read depth at the locus"'
			},{
				'ID' => 'AF',
				'Number' => 'A',
				'Type' => 'Float',
				'Description' => '"Estimated allele frequency in the range (0,1]"'
			}],
			'FORMAT' => [{
				'ID' => 'GT',
				'Number' => '1',
				'Type' => 'String',
				'Description' => '"Genotype"'
			},{
				'ID' => 'DP',
				'Number' => '1',
				'Type' => 'Integer',
				'Description' => '"Read Depth"'
			}],
			'other' => {
				'SnpEffVersion' => '"4.2 (build 2015-12-05), by Pablo Cingolani"'
			},
			'columns' => ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO',
				'FORMAT','Sample']
		};
		
		# Define two simple wrong records
		my $r = [{
			'QUAL' => 'A',
			'INFO' => {
				'AB' => ['0'],
				'ABP' => ['0'],
				'AC' => ['2'],
				'AF' => ['1']
			},
			'ALT' => ['G'],
			'ID' => ['.'],
			'SAMPLE' => [{
				'name' => 'Sample',
				'order' => 0,
				'GT' => ['1/1'],
				'GQ' => ['102.247'],
				'DP' => ['12'],
				'DPR' => ['12','12']
			}],
			'REF' => 'A',
			'FORMAT' => ['GT','GQ','DP','DPR'],
			'CHROM' => '1 X',
			'POS' => 69511,
			'FILTER' => ['.']
		},{
			'QUAL' => '261.217',
			'INFO' => {
			'AB' => ['0'],
			'ABP' => ['0'],
			'AC' => ['3'],
			'AF' => ['1']
			},
			'ALT' => ['T'],
			'ID' => ['.'],
			'SAMPLE' => [{
				'name' => 'Sample',
				'order' => 0,
				'GT' => ['A'],
				'GQ' => ['100.247'],
				'DP' => ['13'],
				'DPR' => ['13','13']
			}],
			'REF' => 'C',
			'FORMAT' => ['GT','GQ','DP','DPR'],
			'CHROM' => '1',
			'POS' => 12345,
			'FILTER' => ['.']
		}];


        my $ho = VCF::Header->new({
            'data' => $h
        });
        my $ro = VCF::Record->new({
            'data' => $r
        });
        
        dies_ok(sub {$ho->validate},"Init VCF::Header from wrong data failed");
        dies_ok(sub {$ro->validate_all},
			"Init VCF::Record from wrong data failed");
    };
}

