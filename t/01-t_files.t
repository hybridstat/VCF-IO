#!perl -T
use 5.10.0;
use strict;
use warnings;
use Test::More;
use Test::File;

plan tests => 5;

BEGIN {
    file_readable_ok('./t/t1.vcf',"Test file t1.vcf readable");
    file_readable_ok('./t/t2.vcf',"Test file t2.vcf readable");
    file_not_exists_ok('./t/foo.vcf',"Test file foo.vcf does not exist");
    file_line_count_is('./t/t1.vcf',180,"Test file t1.vcf #lines ok");
    file_line_count_is('./t/t2.vcf',250,"Test file t2.vcf #lines ok");
}

