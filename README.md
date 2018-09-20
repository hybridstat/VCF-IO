# VCF::IO
A Perl module for parsing, manipulating and writing VCF files

## TL;DR

Read, validate, manipulate and write VCF files, headers and records, either in
an OO manner or via exported functions

### Quick installation

Install the package

```
git clone https://github.com/hybridstat/VCF-IO.git
cd VCF-IO
sudo sh build.sh
```

Start playing around

```
use VCF::IO;

# OO usage
my $vcf = VCF::IO->new({
    'file' => 'my_variants.vcf'
});
$vcf->parse_and_validate;

# Non-OO usage
my ($header,$records) = &parse_and_validate_vcf_file($file);
```

## Introduction

VCF::IO is a Perl module which provides facilities to read, manipulate and 
write Variant Call Format (VCF) files which are usually the output of algorithms
that detect variations in genomes as compared to a reference genome. The
explanation of the procedure is beyond the scope of this introduction and some
(or a bit more) knowledge of Bioinformatics is required to use the module.

VCF::IO can be used to easily parse, validate, manipulate and write VCF files in
a Perl friendly way. It provides methods for extensive VCF record parsing,
supporting and respecting both small and structural variant VCF specifications
as well as creating, validating and appending new records, or modifying existing
ones. A lot of attention is paid to the validation of VCF header and
meta-information as well as the VCF records themselves. 

VCF::IO supports VCF formats 4.1, 4.2, 4.3 which are the most widely used
versions at this point. Support for most up-to-date versions is on the way.

Format specifications can be found in the following links

| Version | Specifications |
| ------- | ----------------------------------------------------------------- |
| 4.1     | [VCFv4.1 specs](https://samtools.github.io/hts-specs/VCFv4.1.pdf) |
| 4.2     | [VCFv4.2 specs](https://samtools.github.io/hts-specs/VCFv4.2.pdf) |
| 4.3     | [VCFv4.3 specs](https://samtools.github.io/hts-specs/VCFv4.3.pdf) |

Unlike other Perl (or other language) packages/scripts/suites used for VCF
processing and downstream bioinformatics analysis, VCF::IO focuses a lot on VCF
header and record validation, as a lot of effort and manual labor is being put 
into bringing the extremely flexible VCF format in a form that can be used by
other tools. As the VCF format is very flexible and sometimes the specifications
are ambiguous or incomplete, the validation routines of the module are doing as
much validation as possible. This may sometimes be even not desirable and cause]
problems. For example, older versions of the FreeBayes variant caller output a
FORMAT field which does not respect the specs 
(https://github.com/ekg/freebayes/pull/323). For this reason, VCF::IO allows
for multiple validation levels with varying strictness.

## Quick examples

Some quick examples following illutrating basic functionalities. For a full 
documentation either use the ```perldoc``` function after installation:

```perldoc VCF::IO```

or read the [Wiki](https://github.com/hybridstat/VCF-IO/wiki) pages (essentially 
converted from POD documentation using ```Pod::Markdown```)

### Parse and validate a VCF file

```
use VCF::IO;

# OO usage
my $vcf = VCF::IO->new({
    'file' => 'my_variants.vcf'
});
$vcf->parse_and_validate;

# Non-OO usage
my ($header,$records) = &parse_and_validate_vcf_file($file);
```

### Parse and validate VCF header

```
use VCF::Header;
    
# Parse from file
my $header = VCF::Header->new($file);
$header->parse; # Simple parsing
$header->parse_and_validate; # Parsing and strict validation

# Parse from data
my $hdata = $header->get_header;
my $new_header = VCF::Header->new({
    'data' => $hdata,
    'validation' => "relaxed"
});
```

### Parse and validate VCF records

```
use VCF::Header;
use VCF::Record;

# Parse from file
my $record = VCF::Record->new($file);

# Simple parsing
$record->parse; 

# Parsing and strict validation - requires VCF::Header
$header->parse_and_validate->get_header;
$record->parse_and_validate($header); 

# Parse from data
my $rdata = $record->get_records;
my $new_record = VCF::Record->new({
    'data' => $rdata,
    'validation' => "relaxed"
});
```

### Write a VCF file

```
use VCF::IO qw(parse_and_validate_vcf_file write_vcf_file);
 
# OO usage
my $vcf = VCF::IO->new({
    'file' => "my_variants.vcf"
});
$vcf->parse_and_validate;

# Do some manipulation...
my $records = $vcf->get_records;
my $oblivion = shift @{$records};
$vcf->set_records($records);

# Write new file
$vcf->output("my_filtered_variants.vcf");
$vcf->write;
$vcf->done;

# Non OO usage
my ($header,$records) = &parse_and_validate_vcf_file("my_variants.vcf");

# Do some manipulation
my $oblivion = shift @{$records};
$vcf->set_records($records);

# Write new file
&write_vcf_file($header,$records,"my_filtered_variants.txt");
```

## Further work

Although the package is functional at this first present state, there are
a lot of functionalities to add both regarding VCF validation and
manipulatio as well as documentation (related _TODO_ tags). Stay tuned!
