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

plan tests => 6;

BEGIN {
    subtest 'Check record parsing' => sub {
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
        
        my $r1 = VCF::Record->new({
            'file' => $t1
        });
        my $r2 = VCF::Record->new({
            'file' => $t2
        });
        
        my $hh1 = $h1->get_header;
        my $hh2 = $h2->get_header;
        
        dies_ok(sub { $r1->validate_all($hh1); },
            "SNP records just validate fails");
        lives_ok(sub { $r1->parse_all; },"SNP records parse ok");
        lives_ok(sub {
            $r1->parse_all;
            $r1->validate_all($hh1);
        },"SNP records parse, validate ok");
        lives_ok(sub { $r1->parse_and_validate_all($hh1); },
            "SNP records parse and validate ok");
        
        dies_ok(sub { $r2->validate_all($hh2); },
            "SNP records just validate fails");
        lives_ok(sub { $r2->parse_all; },"SNP records parse ok");
        lives_ok(sub {
            $r2->parse_all;
            $r2->validate_all($hh2);
        },"SNP records parse, validate ok");
        lives_ok(sub { $r2->parse_and_validate_all($hh2); },
            "SNP records parse and validate ok");
    };    
    
    subtest 'Check record values' => sub {
        plan tests => 21;
        
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
        
        my $r1 = VCF::Record->new({
            'file' => $t1
        });
        my $r2 = VCF::Record->new({
            'file' => $t2
        });
        
        my $hh1 = $h1->get_header;
        my $hh2 = $h2->get_header;
        
        $r1->parse_and_validate_all($hh1);
        $r2->parse_and_validate_all($hh2);
        
        my $rr1 = $r1->get_records->[0];
        my $rr2 = $r2->get_records->[0];
        
        is($r1->length,29,"SNP record length ok");
        is($r2->length,216,"Structural record length ok");
        
        is($rr1->{"CHROM"},"1","SNP record CHROM ok");
        is($rr1->{"POS"},69511,"SNP record CHROM ok");
        is($rr1->{"REF"},"A","SNP record REF ok");
        is($rr1->{"ALT"}->[0],"G","SNP record ALT ok");
        is($rr1->{"INFO"}->{"AB"}->[0],"0","SNP record INFO AB ok");
        is($rr1->{"INFO"}->{"AC"}->[0],"2","SNP record INFO AC ok");
        is($rr1->{"FORMAT"}->[0],"GT","SNP record FORMAT GT ok");
        is($rr1->{"FORMAT"}->[1],"GQ","SNP record FORMAT GQ ok");
        is($rr1->{"SAMPLE"}->[0]->{"GT"}->[0],"1/1","SNP record SAMPLE GT ok");
        is($rr1->{"SAMPLE"}->[0]->{"DP"}->[0],"12","SNP record SAMPLE DP ok");
        is($rr1->{"_id"},"1^^^69511^^^A^^^G","SNP record _id ok");
                
        is($rr2->{"CHROM"},"1","Structural record CHROM ok");
        is($rr2->{"POS"},1,"Structural record POS ok");
        is($rr2->{"REF"},"N","Structural record REF ok");
        is($rr2->{"ALT"}->[0],"<DEL>","Structural record ALT ok");
        is($rr2->{"INFO"}->{"SVTYPE"}->[0],"DEL",
            "Structural record INFO SVTYPE ok");
        is($rr2->{"INFO"}->{"END"}->[0],"2857518",
            "Structural record INFO END ok");
        is(scalar(@{$rr2->{"FORMAT"}}),0,"Structural record FORMAT ok");
        is($rr2->{"_id"},"1^^^1^^^N^^^<DEL>","Structural record _id ok");
    };
    
    subtest 'Check record key-value attribute add' => sub {
        plan tests => 3;
        
        my $t1 = "./t/t1.vcf"; # SNPs
        
        my $h1 = VCF::Header->new({
            'file' => $t1
        });
        $h1->parse_and_validate;
        
        my $r1 = VCF::Record->new({
            'file' => $t1
        });
        my $hh1 = $h1->get_header;
        $r1->parse_and_validate_all($hh1);
        my $rr1 = $r1->get_records->[0];
        
        my $type = "INFO";
        my %attr1 = (
            'ID' => "MYANN",
            'Number' => ".",
            'Type' => "String",
            'Description' => "My super annotation"
        );
        my %attr2 = (
            'ID' => "MYSCORE",
            'Number' => "1",
            'Type' => "Float",
            'Description' => "My ultra score"
        );
        $h1->add(\%attr1,$type);
        $h1->add(\%attr2,$type);
        
        $hh1 = $h1->get_header;
        
        lives_ok(sub {
            my %kv1 = (
                'MYANN' => ["A|geneA|dis1^dis2^dis3|phen1^phen2^phen3"],
                'MYSCORE' => [0.458]
            );
            $rr1 = $r1->add_key_value_attr(\%kv1,$type,$hh1,$rr1);
        },"SNP record key-value add ok");
        is($rr1->{"INFO"}->{"MYSCORE"}->[0],0.458,"SNP record new MYSCORE ok");
        
        dies_ok(sub {
            my %kv2 = (
                'MYANN' => ["A|geneA|dis1^dis2^dis3|phen1^phen2^phen3"],
                'MYSCORE' => [0.458,0.789]
            );
            $rr1 = $r1->add_key_value_attr(\%kv2,$type,$hh1,$rr1);
        },"SNP record key-value add safe fail ok");
    };
    
    subtest 'Check record key attribute add' => sub {
        plan tests => 6;
        
        my $t1 = "./t/t1.vcf"; # SNPs
        
        my $h1 = VCF::Header->new({
            'file' => $t1
        });
        $h1->parse_and_validate;
        
        my $r1 = VCF::Record->new({
            'file' => $t1
        });
        my $hh1 = $h1->get_header;
        $r1->parse_and_validate_all($hh1);
        my $rr1 = $r1->get_records->[0];
        
        my $type1 = "FILTER";
		my $type2 = "FORMAT";
		my %fil1 = (
			'ID' => "q10",
			'Description' => "q10 filter"
		);
		my %fil2 = (
			'ID' => "p20",
			'Description' => "p20 filter"
		);
		my %for1 = (
			'ID' => "BA",
			'Number' => "A",
			'Type' => "Integer",
			'Description' => "Format 1"
		);
		my %for2 = (
			'ID' => "KA",
			'Number' => "1",
			'Type' => "Float",
			'Description' => "Format 2"
		);
		$h1->add(\%fil1,$type1);
		$h1->add(\%fil2,$type1);
		$h1->add(\%for1,$type2);
		$h1->add(\%for2,$type2);
        
        $hh1 = $h1->get_header;
        
        lives_ok(sub {
            my @k1 = ("q10","p20");
			my @f1 = ("BA","KA");
			$rr1 = $r1->add_key_attr(\@k1,$type1,$hh1,$rr1);
			$rr1 = $r1->add_key_attr(\@f1,$type2,$hh1,$rr1);
        },"SNP record key add ok");
        is($rr1->{"FILTER"}->[0],"q10","SNP record new FILTER ok");
        
        dies_ok(sub {
            my $rr2 = $r1->get_records->[0];
			my @k1 = ("q10","p30");
			$rr2 = $r1->add_key_attr(\@k1,$type1,$hh1,$rr2);
        },"SNP record key add safe fail 1 ok");
        
        dies_ok(sub {
            my $rr2 = $r1->get_records->[0];
			my @k1 = ("q10","q10");
			$rr2 = $r1->add_key_attr(\@k1,$type1,$hh1,$rr2);
        },"SNP record key add safe fail 2 ok");
        
        dies_ok(sub {
            my $rr2 = $r1->get_records->[0];
			my @f1 = ("BA","DA");
			$rr2 = $r1->add_key_attr(\@f1,$type2,$hh1,$rr2);
        },"SNP record key add safe fail 3 ok");
        
        dies_ok(sub {
            my $rr2 = $r1->get_records->[0];
			my @f1 = ("BA","BA");
			$rr2 = $r1->add_key_attr(\@f1,$type2,$hh1,$rr2);
        },"SNP record key add safe fail 4 ok");
    };
    
    subtest 'Check new record add' => sub {
        plan tests => 5;
        
        my $t1 = "./t/t1.vcf"; # SNPs
        
        my $h1 = VCF::Header->new({
            'file' => $t1
        });
        $h1->parse_and_validate;
        
        my $r1 = VCF::Record->new({
            'file' => $t1
        });
        my $hh1 = $h1->get_header;
        $r1->parse_and_validate_all($hh1);
        
        lives_ok(sub {
			my $good_rec = {
				'ALT' => ['A'],
				'POS' => 752776,
				'FORMAT' => ['GT','GQ','DP','DPR','RO','QR','AO','QA','GL'],
				'FILTER' => ['.'],
				'SAMPLE' => [{
					'name' => 'Sample',
					'order' => 0,
					'GT' => ['0/1'],
				  'GQ' => ['122.1854'],
				  'DP' => ['12'],
				  'DPR' => ['12','12'],
				  'RO' => ['0'],
				  'QR' => ['0'],
				  'AO' => ['10'],
				  'QA' => ['544'],
				  'GL' => ['-32.3168','-1.70927','0']
				}],
				'REF' => 'T',
				'CHROM' => '1',
				'ID' => ['.'],
			  'INFO' => {
				  'AB' => ['0'],
				  'ABP' => ['0'],
				  'AC' => ['3'],
				  'AF' => ['1'],
				  'AN' => ['2'],
				  'AO' => ['11'],
				  'CIGAR' => ['1X'],
				  'DP' => ['10'],
				  'DPB' => ['9'],
				  'DPRA' => ['0'],
				  'EPP' => ['10.543'],
				  'EPPR' => ['0'],
				  'GTI' => ['0'],
				  'LEN' => ['1'],
				  'MEANALT' => ['1'],
				  'MQM' => ['61'],
				  'MQMR' => ['0'],
				  'NS' => ['1'],
				  'NUMALT' => ['1'],
				  'ODDS' => ['19.3844'],
				  'PAIRED' => ['1'],
				  'PAIREDR' => ['0'],
				  'PAO' => ['0'],
				  'PQA' => ['0'],
				  'PQR' => ['0'],
				  'PRO' => ['0'],
				  'QA' => ['344'],
				  'QR' => ['0'],
				  'RO' => ['0'],
				  'RPL' => ['2'],
				  'RPP' => ['9.04217'],
				  'RPPR' => ['0'],
				  'RPR' => ['7'],
				  'RUN' => ['1'],
				  'SAF' => ['9'],
				  'SAP' => ['22.5536'],
				  'SAR' => ['0'],
				  'SRF' => ['0'],
				  'SRP' => ['0'],
				  'SRR' => ['0'],
				  'TYPE' => ['snp'],
				  'technology.ILLUMINA' => ['1'],
				  'ANN' => [
					'G|upstream_gene_variant|MODIFIER|FAM87B|ENSG00000177757|'.
					'transcript|ENST00000326734|lincRNA||n.-30A>G|||||30|',
					'G|intron_variant|MODIFIER|RP11-206L10.10|ENSG00000240453|'.
					'transcript|ENST00000435300|processed_transcript|1/1|'.
					'n.340+32T>C||||||'
				]},
				'QUAL' => '1063.016'
			};
			$r1->add($good_rec,$hh1);
		},"Add record ok");
		is($r1->length,30,"New record length ok");
		is($r1->get_records->[29]->{"QUAL"},"1063.016","New record QUAL ok");
		is($r1->get_records->[29]->{"POS"},752776,"New record POS ok");
		
		dies_ok(sub {
			my $bad_rec = {
				'ALT' => ['2'],
				'POS' => "752776",
				'FORMAT' => ['GT','GQ','DP','DPR','RO','QR','AO','QA','GL'],
				'FILTER' => ['.'],
				'SAMPLE' => [{
				  'name' => 'Mourlas_Outsource',
				  'order' => 0,
				  'GT' => ['0/1'],
				  'GQ' => ['122.1854'],
				  'DP' => ['12'],
				  'DPR' => ['12','12'],
				  'RO' => ['0'],
				  'QR' => ['0'],
				  'AO' => ['10'],
				  'QA' => ['544'],
				  'GL' => ['-32.3168','-1.70927','0']
				}],
				'REF' => 'T',
				'CHROM' => '1',
				'ID' => ['.'],
				'INFO' => {
					'AB' => ['0'],
					'ABP' => ['0'],
					'AC' => ['3'],
					'AF' => ['1'],
					'AN' => ['2'],
					'AO' => ['11'],
					'CIGAR' => ['1X'],
					'DP' => ['10'],
					'DPB' => ['9'],
					'DPRA' => ['0'],
					'EPP' => ['10.543'],
					'EPPR' => ['0'],
					'GTI' => ['0'],
					'LEN' => ['1'],
					'MEANALT' => ['1'],
					'MQM' => ['61'],
					'MQMR' => ['0'],
					'NS' => ['1'],
					'NUMALT' => ['1'],
					'ODDS' => ['19.3844'],
					'PAIRED' => ['1'],
					'PAIREDR' => ['0'],
					'PAO' => ['0'],
					'PQA' => ['0'],
					'PQR' => ['0'],
					'PRO' => ['0'],
					'QA' => ['344'],
					'QR' => ['0'],
					'RO' => ['0'],
					'RPL' => ['2'],
					'RPP' => ['9.04217'],
					'RPPR' => ['0'],
					'RPR' => ['7'],
					'RUN' => ['1'],
					'SAF' => ['9'],
					'SAP' => ['22.5536'],
					'SAR' => ['0'],
					'SRF' => ['0'],
					'SRP' => ['0'],
					'SRR' => ['0'],
					'TYPE' => ['snp'],
					'technology.ILLUMINA' => ['1'],
					'ANN' => [
						'G|upstream_gene_variant|MODIFIER|FAM87B|'.
						'ENSG00000177757|transcript|ENST00000326734|lincRNA|',
						'|n.-30A>G|||||30|G|intron_variant|MODIFIER|'.
						'RP11-206L10.10|ENSG00000240453|transcript|'.
						'ENST00000435300|processed_transcript|'.
						'1/1|n.340+32T>C||||||'
					]},
				'QUAL' => '1063.016'
			};
			$r1->add($bad_rec,$hh1);
		},"Add record safe fail ok");
    };
    
    subtest 'Check record writing' => sub {
        plan tests => 4;
        
        my $t1 = "./t/t1.vcf"; # SNPs
        
        my $h1 = VCF::Header->new({
            'file' => $t1
        });
        $h1->parse_and_validate;
        
        my $r1 = VCF::Record->new({
            'file' => $t1
        });
        my $hh1 = $h1->get_header;
        $r1->parse_and_validate_all($hh1);
        
        lives_ok(sub { 
            $r1->output("./t/t1r.tmp");
            $r1->write;
            $r1->done;
        },"Successful records write");
        file_readable_ok('./t/t1r.tmp',"Test file t1r.tmp readable");
        file_min_size_ok('./t/t1r.tmp',35000,"t1r.tmp min size ok");
        file_max_size_ok('./t/t1r.tmp',45000,"t1r.tmp max size ok");
        
        unlink("./t/t1r.tmp");
    };
}

