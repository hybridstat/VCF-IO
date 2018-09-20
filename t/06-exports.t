#!perl -T
use 5.10.0;
use strict;
use warnings;
use Test::More;
use Test::Exception;
use Test::File;

use VCF::IO qw(parse_vcf_header validate_vcf_header 
    parse_and_validate_vcf_header parse_vcf_sample_names write_vcf_header
    parse_vcf_records validate_vcf_records parse_and_validate_vcf_records
    write_vcf_records parse_vcf_file parse_and_validate_vcf_file 
    write_vcf_file);
use VCF::Header;
use VCF::Record;
use VCF::_Helper;

plan tests => 2;

BEGIN {
    subtest 'Check SNP VCF file parsing' => sub {
        plan tests => 21;
        
        my $t1 = "./t/t1.vcf";
        
        my ($h1,$s1,$r1);

		lives_ok(sub {
			$h1 = &parse_vcf_header($t1);
		},"SNP exported header parsing ok");
		lives_ok(sub {
			$h1 = &validate_vcf_header($h1);
		},"SNP exported header validation ok");
		lives_ok(sub {
			$h1 = &parse_and_validate_vcf_header($t1);
		},"SNP exported header parsing and validation ok");
		lives_ok(sub {
			$s1 = &parse_vcf_sample_names($t1);
		},"SNP exported sample name parsing ok");
		lives_ok(sub {
			&write_vcf_header($h1,"./t/t1eh.txt");
		},"SNP exported header writing ok");
		file_readable_ok('./t/t1eh.txt',"Test file t1eh.txt readable");
        file_min_size_ok('./t/t1eh.txt',10000,"t1eh.txt min size ok");
        file_max_size_ok('./t/t1eh.txt',15000,"t1eh.txt max size ok");
		unlink('./t/t1eh.txt');
		
		lives_ok(sub {
			$r1 = &parse_vcf_records($t1)
		},"SNP exported record parsing ok");
		lives_ok(sub {
			$r1 = &validate_vcf_records($r1,$h1);
		},"SNP exported record validation ok");
		lives_ok(sub {
			$r1 = &parse_and_validate_vcf_records($t1,$h1);
		},"SNP exported record parsing and validation ok");
		lives_ok(sub {
			&write_vcf_records($r1,"./t/t1er.txt");
		},"SNP exported record writing ok");
		file_readable_ok('./t/t1er.txt',"Test file t1er.txt readable");
        file_min_size_ok('./t/t1er.txt',35000,"t1er.txt min size ok");
        file_max_size_ok('./t/t1er.txt',45000,"t1er.txt max size ok");
		unlink('./t/t1er.txt');

		lives_ok(sub {
			($h1,$r1) = &parse_vcf_file($t1);
		},"SNP exported file parsing ok");
		lives_ok(sub {
			($h1,$r1) = &parse_and_validate_vcf_file($t1);
		},"SNP exported file parsing and validation ok");
		lives_ok(sub {
			&write_vcf_file($h1,$r1,"./t/t1ehv.vcf");
		},"SNP exported file writing ok");
		file_readable_ok('./t/t1ehv.vcf',"Test file t1ehv.vcf readable");
        file_min_size_ok('./t/t1ehv.vcf',50000,"t1ehv.vcf min size ok");
        file_max_size_ok('./t/t1ehv.vcf',60000,"t1ehv.vcf max size ok");
		unlink('./t/t1ehv.vcf');
    };
    
    subtest 'Check Structural VCF file parsing' => sub {
        plan tests => 21;
        
        my $t2 = "./t/t2.vcf";
        
        my ($h2,$s2,$r2);

		lives_ok(sub {
			$h2 = &parse_vcf_header($t2);
		},"SNP exported header parsing ok");
		lives_ok(sub {
			$h2 = &validate_vcf_header($h2);
		},"SNP exported header validation ok");
		lives_ok(sub {
			$h2 = &parse_and_validate_vcf_header($t2);
		},"SNP exported header parsing and validation ok");
		lives_ok(sub {
			$s2 = &parse_vcf_sample_names($t2);
		},"SNP exported sample name parsing ok");
		lives_ok(sub {
			&write_vcf_header($h2,"./t/t2eh.txt");
		},"SNP exported header writing ok");
		file_readable_ok('./t/t2eh.txt',"Test file t2eh.txt readable");
        file_min_size_ok('./t/t2eh.txt',3000,"t2eh.txt min size ok");
        file_max_size_ok('./t/t2eh.txt',4000,"t2eh.txt max size ok");
		unlink('./t/t2eh.txt');
		
		lives_ok(sub {
			$r2 = &parse_vcf_records($t2)
		},"SNP exported record parsing ok");
		lives_ok(sub {
			$r2 = &validate_vcf_records($r2,$h2);
		},"SNP exported record validation ok");
		lives_ok(sub {
			$r2 = &parse_and_validate_vcf_records($t2,$h2);
		},"SNP exported record parsing and validation ok");
		lives_ok(sub {
			&write_vcf_records($r2,"./t/t2er.txt");
		},"SNP exported record writing ok");
		file_readable_ok('./t/t2er.txt',"Test file t2er.txt readable");
        file_min_size_ok('./t/t2er.txt',40000,"t2er.txt min size ok");
        file_max_size_ok('./t/t2er.txt',50000,"t2er.txt max size ok");
		unlink('./t/t2er.txt');

		lives_ok(sub {
			($h2,$r2) = &parse_vcf_file($t2);
		},"SNP exported file parsing ok");
		lives_ok(sub {
			($h2,$r2) = &parse_and_validate_vcf_file($t2);
		},"SNP exported file parsing and validation ok");
		lives_ok(sub {
			&write_vcf_file($h2,$r2,"./t/t2ehv.vcf");
		},"SNP exported file writing ok");
		file_readable_ok('./t/t2ehv.vcf',"Test file t2ehv.vcf readable");
        file_min_size_ok('./t/t2ehv.vcf',40000,"t2ehv.vcf min size ok");
        file_max_size_ok('./t/t2ehv.vcf',55000,"t2ehv.vcf max size ok");
		unlink('./t/t2ehv.vcf');
    };
}
