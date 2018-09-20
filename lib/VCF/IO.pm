package VCF::IO;

use 5.10.0;
use strict;
use warnings;
use utf8;

use Carp;
use File::Basename;
use File::Spec;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Tie::IxHash::Easy;

use VCF::Header;
use VCF::Record;
use VCF::_Helper;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(parse_vcf_file parse_and_validate_vcf_file write_vcf_file
    parse_vcf_header validate_vcf_header parse_and_validate_vcf_header
    parse_vcf_sample_names write_vcf_header
    parse_vcf_records validate_vcf_records parse_and_validate_vcf_records
    write_vcf_records);

use vars qw($helper);

BEGIN {
    $helper = VCF::_Helper->new;
    binmode STDOUT, ":utf8";
    $|=1;
}

use Data::Dumper qw(Dumper);

=head1 NAME

VCF::IO - A Perl module for parsing, manipulating and writing VCF files.

=head1 VERSION

Version 0.01

=cut

our $MODNAME = "VCF::IO";
our $VERSION = '0.01';

=head1 SYNOPSIS

VCF::IO is used to easily parse, validate, manipulate and write VCF files in
a Perl friendly way. It provides methods for extensive VCF record parsing,
supporting and respecting both small and structural variant VCF specifications
as well as creating, validating and appending new records, or modifying existing
ones. A lot of attention is paid to the validation of VCF header and
meta-information as well as the VCF records themselves. As the VCF format is
very flexible and sometimes the specifications are ambiguous or incomplete, the
validation routines of the module are doing as much validation as possible.
This may cause problems. For example, older versions of the FreeBayes variant
caller output a FORMAT field which does not respect the specs 
(https://github.com/ekg/freebayes/pull/323). For this reason, VCF::IO allows
for multiple validation levels with varying strictness.

There are three validation levels in VCF::IO (and its children modules,
VCF::Header and VCF::Record):
 
=over

=item *

B<basic>: When validation is C<'basic'>, only very essential validation
procedures are performed both for header and records, such as requiring that
fields are not missing (e.g. consecutive tabs) and basic content checks (e.g.
that the required attributes do exist in the VCF header and meta-info lines).

=item *

B<relaxed>: When validation is C<'relaxed'>, additional checks are performed in
both VCF header and records. For example, while when validation is C<'basic'>
the chromosome name is not checked for naming consistency according to the 
specs, with C<'relaxed'> validation, this check is performed. However, slight
deviations from the specs are allowed (e.g. chromosome may not exist in the
##contig meta-info lines in the case of relaxed record validation, or the 
##contig meta-info lines may not exist at all in the case of the header). The
relaxed validation may prove useful in cases of heavily customized applications
which rely on the basics of the VCF format though.

=item *
 
B<'strict'>: When validation is C<'strict'>, VCF::IO performs as many validation
steps as possible in order to comply with each VCF format version specs.
 
=back

If validation parameter is not provided, the strict validation is performed by
default.

VCF::IO also provides functions that allow the usage of the module in an Object
Oriented or a more script-like manner.


    use VCF::IO;

    # OO usage
    my $vcf = VCF::IO->new({
        'file' => 'my_variants.vcf',
        'validation' => 'relaxed'
    });
    $vcf->parse_and_validate;
    
    # Non-OO usage
    my ($header,$records) = &parse_and_validate_vcf_file($file);

=head1 EXPORT

parse_vcf_file
parse_and_validate_vcf_file
write_vcf_file
parse_vcf_header
validate_vcf_header
parse_and_validate_vcf_header
parse_vcf_sample_names
write_vcf_header
parse_vcf_records
validate_vcf_records
parse_and_validate_vcf_records
write_vcf_records

=cut

# All the exported functions that use the interfaces should live in this section

=head2 parse_vcf_file

Simple non-OO parsing of a VCF file, returning two references: the first is a
hash reference containing the VCF header and meta-info line and the second is
an array reference pointing to an array of hashes with the actual VCF records. 
This function does not perform any essential validations.

    my $file = 'my_variants.vcf';
    my ($header,$records) = &parse_vcf_file($file);

=cut

sub parse_vcf_file {
    my $file = $_[0];
    my $header = &parse_vcf_header($file);
    my $records = &parse_vcf_records($file,$header);
    return($header,$records);
}

=head2 parse_and_validate_vcf_file

Simple non-OO parsing and validation of a VCF file, returning the same 
(validated) output as C<parse_vcf_file>.

    my $file = 'my_variants.vcf';
    my ($header,$records) = &parse_and_validate_vcf_file($file);

=cut

sub parse_and_validate_vcf_file {
    my $file = $_[0];
    my $header = &parse_and_validate_vcf_header($file);
    my $records = &parse_and_validate_vcf_records($file,$header);
    return($header,$records);
}

=head2 write_vcf_file

Simple non OO writing of a VCF file given a VCF header and some VCF records. 
The first argument must be a VCF::Header object or a valid hash with VCF header 
elements such as the one returned by C<VCF::IO::get_header> or 
C<VCF::Header->get_header>. The second argument must be a VCF::Record object or 
a valid hash with VCF record elements such as the one returned by 
C<VCF::IO::get_records> or C<VCF::Record->get_records>. The third argument must
be a filename to write or a valid filehandle. No validation of the contents is
performed prior to writing! Use other package facilities for this.

    my $file = 'my_variants.vcf';
    my ($header,$records) = &parse_vcf_file($file);
    &write_vcf_file($header,$records,"my_test_out.vcf);

=cut

sub write_vcf_file {
    my ($header,$records,$file) = @_;
    
    # Init the objects if not objects
    $header = VCF::Header->new({'data' => $header}) 
        if (ref($header) ne "VCF::Header");
    $records = VCF::Record->new({'data' => $records})
        if (ref($records) ne "VCF::Record");
    
    # Initialize the filehandle and write the header
    $header->output($file);
    $header->write;
    $header->done;
    # Then write the records
    $records->output($file,"append");
    $records->write;
    $records->done;
    
    return(0);
}

=head2 parse_vcf_header

Simple non-OO parsing of a VCF file header and meta-info lines. This function
does not perform any essential validations.

    my $file = 'my_variants.vcf';
    my $header = parse_vcf_header($file);

=cut

sub parse_vcf_header {
    my $vcf = $_[0];
    
    # Use the OO interface
    my $header = VCF::Header->new($vcf);
    $header->parse;
    
    return($header->get_header);
}

=head2 validate_vcf_header

Simple non-OO validating of a VCF file header and meta-info lines. This function
performs validations according to the supported VCF specifications. It can be 
executed with the result of C<parse_vcf_header>.

    my $file = 'my_variants.vcf';
    my $header = parse_vcf_header($file);
    $header = validate_vcf_header($header);

=cut

sub validate_vcf_header {
    my $h = $_[0];
    
    croak "\nThe provided header hash is not a hash reference!\n"
        if (!$helper->_is_hash_ref($h));
    
    # Use the OO interface
    my $header = VCF::Header->new({'data' => $h});
    $header->validate;
    
    return($header->get_header);
}

=head2 parse_and_validate_vcf_header

Simple non-OO parsing and validation of a VCF file header and meta-info lines.

    my $file = 'my_variants.vcf';
    my $header = parse_and_validate_vcf_header($file);

=cut

sub parse_and_validate_vcf_header {
    my $vcf = $_[0];
    
    # Use the OO interface
    my $header = VCF::Header->new($vcf);
    $header->parse_and_validate;
    
    return($header->get_header);
}

=head2 parse_vcf_sample_names

Simple non-OO parsing of sample names with analyzed variants included in a VCF
file.

    my $file = 'my_variants.vcf';
    my $sample_names = parse_vcf_sample_names($file);

=cut

sub parse_vcf_sample_names {
    my $vcf = $_[0];
    my $header = VCF::Header->new($vcf);
    return($header->parse_sample_names);
}

=head2 write_vcf_header

Simple non-OO writing of a VCF file header and meta-info lines. The first 
argument must be a VCF::Header object or a valid hash with VCF header elements
such as the one returned by C<VCF::IO::get_header> or C<VCF::Header->get_header>
The second argument must be a filename to write or a valid filehandle. No
validation of the contents is performed prior to writing! Use other package
facilities for this.

    write_vcf_header($header,'output.vcf');

=cut

sub write_vcf_header {
    my ($header,$file) = @_;
    
    # Init an object if not an object
    $header = VCF::Header->new({'data' => $header}) 
        if (ref($header) ne "VCF::Header");
    
    # Initialize the filehandle and write
    $header->output($file);
    $header->write;
    $header->done;
}

=head2 parse_vcf_records

Simple non-OO parsing of a VCF file records and return an array of hashes. This
function does not perform any essential record validations.

    my $file = 'my_variants.vcf';
    my $records = parse_vcf_records($file);

=cut

sub parse_vcf_records {
    my $vcf = $_[0];
    
    # Use the OO interface
    my $record = VCF::Record->new($vcf);
    $record->parse_all;
    
    return($record->get_records);
}

=head2 validate_vcf_records

Simple non-OO validating of a VCF file records according to the supported VCF
specifications. It can be executed with the result of C<parse_vcf_records>.
A valid VCF::Header object is required.

    my $file = 'my_variants.vcf';
    my $header = parse_and_validate_vcf_header($file);
    my $records = parse_vcf_records($file);
    $records = validate_vcf_records($records,$header);

=cut

sub validate_vcf_records {
    my ($r,$h) = @_;
    
    croak "\nThe provided records array is not an array reference!\n"
        if (!$helper->_is_array_ref($r));
    croak "\nNo VCF header provided! Cannot proceed with validation.\n"
        if (!$h);
    
    # Use the OO interface
    my $records = VCF::Record->new({'data' => $r});
    $records->validate_all($h);
    
    return($records->get_records);
}

=head2 parse_and_validate_vcf_records

Simple non-OO parsing and validation of a VCF file records. A VCF::Header 
object may be provided, otherwise it will be created from the input VCF file.

    my $file = 'my_variants.vcf';
    my $header = parse_and_validate_vcf_header($file);
    my $records = parse_and_validate_vcf_records($file,$header);
    
    # or just
    my $records = parse_and_validate_vcf_records($file);
    

=cut

sub parse_and_validate_vcf_records {
    my ($vcf,$header) = @_;
    
    # Use the OO interface
    if (!$header) {
		$header = VCF::Header->new($vcf);
		$header->parse_and_validate;
		$header = $header->get_header;
	}
    my $records = VCF::Record->new($vcf);
    $records->parse_and_validate_all($header);
    
    return($records->get_records);
}

=head2 write_vcf_records

Simple non-OO writing of a VCF file records. The first  argument must be a
VCF::Record object or a valid hash with VCF records such as the one returned by
C<VCF::IO::get_records> or C<VCF::Record->get_records>. The second argument
must be a filename to write or a valid filehandle. No validation of the
contents is performed prior to writing! Use other package facilities for this.

    write_vcf_records('output.vcf');

=cut

sub write_vcf_records {
    my ($records,$file) = @_;
    
    # Init an object if not an object
    $records = VCF::Record->new({'data' => $records}) 
        if (ref($records) ne "VCF::Record");
    
    # Initialize the filehandle and write
    $records->output($file);
    $records->write;
    $records->done;
}

################################################################################

=head1 SUBROUTINES/METHODS

=head2 new

The VCF::IO object constructor. It accepts a set of parameters that are
required to run the asigner and get the output. The parameters can also be a
single file name.

    my $vcf = VCF::IO->new({
        'file' => 'my_variant.vcf',
        'validation' => 'strict'
    });
    
    my $vcf = VCF::IO->new('my_variant.vcf');

=cut

sub new {
    my ($class,$params) = @_;
    my $self = {};
    
    # Bless so that we start using methods (including paramter validation)
    bless($self,$class);
    
    # Validate the input parameters
    $params = $self->_check_params($params);
    
    # If all OK, continue with object initialization
    $self->_init($params);
    return($self);
}

=head2 parse

Parses a VCF file without validation regardless of the validation parameter
initialized with.

    my $vcf = VCF::IO->new('my_variant.vcf');
    $vcf->parse;

=cut

sub parse {
    my $self = $_[0];
    
    my $file = $self->get('file');
    
    my $header = VCF::Header->new($file);
    $header->parse;
    $self->set_header($header->get_header);
    
    my $record = VCF::Record->new($file);
    $record->parse_all;
    $self->set_records($record->get_records);
    
    return($self);
}

=head2 validate

Validates a VCF file based on the validation strictness initialized with.
Parsing must have been performed using C<parse>.

    my $vcf = VCF::IO->new('my_variant.vcf');
    $vcf->parse;
    $vcf->validate;

=cut

sub validate {
    my $self = $_[0];
    
    croak "\nNo VCF header found! Have you called VCF::IO->parse first?\n"
		if (!$self->has_header);
	croak "\nNo VCF records found! Have you called VCF::IO->parse first?\n"
		if (!$self->has_header);
    
    my $header = VCF::Header->new({'data' => $self->get_header});
    $header->validate;
    $self->set_header($header->get_header);
    
    my $record = VCF::Record->new({'data' => $self->get_records});
    $record->validate_all($header->get_header);
    $self->set_records($record->get_records);
    
    return($self);
}

=head2 parse_and_validate

Parses and validated a VCF file.

    my $vcf = VCF::IO->new('my_variant.vcf');
    $vcf->parse_and_validate;

=cut

sub parse_and_validate {
    my $self = $_[0];
    
    my $file = $self->get('file');
    
    my $header = VCF::Header->new($file);
    $header->parse_and_validate;
    $self->set_header($header->get_header);
    
    my $record = VCF::Record->new($file);
    $record->parse_and_validate_all($header->get_header);
    $self->set_records($record->get_records);
    
    return($self);
}

=head2 write

VCF::IO writing to file. The C<output> method must be called first or a
valid filehandle must be provided.

    $vcf->write($fh);

=cut

sub write {
    my ($self,$fh) = @_;
    
    # Even if $self->output has been called with a file handle and a 
    # file handle is provided, the second precedes. If $self->output has not
    # been called, and $fh is not a valid file handle output will be written
    # to STDOUT
    $self->{'_output'} = $fh if ($fh && fileno($fh) != -1);
    
    my $header_data = $self->get_header;
    my $header = VCF::Record->new({'data' => $header_data});
    my $record_data = $self->get_records;
    my $record = VCF::Record->new({'data' => $record_data});
    
    # Initialize the filehandle and write the contents
    $header->write($fh);
    $record->write($fh);
    $self->done;
}

=head2 output

Opens a file for writing the parsed/manipulated VCF contents. Essentially a 
wraper over Perl's C<open>.

    $vcf->output('output.vcf');

=cut

sub output {
    my ($self,$output) = @_;
    croak "\nNo file provided for writing!\n" if (!$output);
    open(my $fh,">",$output) or croak "Could not open $output for writing!";
    $self->{'_output'} = $fh;
    return($self);
}

=head2 done

Closes the file where VCF content is written. Essentially a  wraper over Perl's
C<close>.

    $vcf->done;

=cut

sub done {
    my $self = $_[0];
    my $fh = $self->{'_output'};
    if (!$fh) {
        carp "\nNo filehandle to close! Returning...";
    }
    else {
        close($fh);
    }
    return(0);
}

=head2 has_header

Returns a true value if the VCF::IO object has a parsed header.

    $vcf->has_header;

=cut

sub has_header {
    my $self = $_[0];
    return(1) if (defined($self->{"_header"}));
    return(0);
}

=head2 has_records

Returns a true value if the VCF::IO object has parsed records.

    $vcf->has_records;

=cut

sub has_records {
    my $self = $_[0];
    return(1) if (defined($self->{"_records"}));
    return(0);
}

=head2 get_header

Returns a simple hash with the VCF header fields.

    my $header = $vcf->get_header;
    
TODO: Add documentation (the field names)

=cut

sub get_header {
    my $self = $_[0];
    return($self->{"_header"});
}

=head2 get_records

Returns a simple array of hashes with the VCF record fields. Each hash in the
array represents one VCF record.

    my $records = $vcf->get_records;
    
TODO: Add documentation (the field names)

=cut

sub get_records {
    my $self = $_[0];
    return($self->{"_records"});
}

=head2 get

VCF::IO object getter

    my $param_value = $vcf->get("param_name");

TODO: Validate that the parameter to be got is valid

=cut

sub get {
    my ($self,$name) = @_;
    return($self->{$name});
}

=head2 set_header

Sets the header content of the VCF file (a VCF::Header object) either after
parsing or constructed manually.
    
    # Code example
    
=cut

sub set_header {
    my ($self,$header) = @_;
    $self->{"_header"} = $header;
    return($self);
}

=head2 set_records

Sets the records content of the VCF file (an array of VCF::Record objects) 
either after parsing or constructed manually.
    
    # Code example
    
=cut

sub set_records {
    my ($self,$records) = @_;
    $self->{"_records"} = $records;
    return($self);
}

=head2 set

VCF::IO object setter

    $vcf->set("param_name","param_value");
    
TODO: Validate that the parameter to be set is valid
    
=cut

sub set {
    my ($self,$name,$value) = @_;
    $self->{$name} = $value;
    return($self);
}

sub _init {
    my ($self,$params) = @_;
    
    foreach my $p (keys(%{$params})) { 
        $self->set($p,$params->{$p}); 
    }
    $self->{"_header"} = undef;
    $self->{"_records"} = undef;
    $self->{"_output"} = undef;
    
    return($self);
}

sub _check_params {
    my ($self,$params) = @_;
    
    # $params can be a hash reference with parameters or an existing filename
    # If the second, create the $params with this file
    if (-f $params) {
        $params = {'file' => $params};
    }
    
    my @accept = ("file","validation");
    
    # Check fatal
    my $stop;
    $stop .= "--- Please specify VCF file to read/write ---\n" 
        if (!$params->{"file"});
        
    if ($stop) {
        print STDERR "\n$stop\n";
        croak "Type perldoc $MODNAME for help in usage.\n\n";
        exit;
    }
    
    # Check and warn for unrecognized parameters
    foreach my $p (keys(%{$params})) {
        carp "\nUnrecognized parameter : $p   --- Ignoring..." 
            if (!$helper->_smatch($p,@accept));
    }
    
    # Check if validation provided and in proper ranges
    if ($params->{"validation"}) {
        croak "\nvalidation parameter must be either 'strict', 'relaxed' or ".
            "'basic'!\n" if ($params->{"validation"} ne "strict" 
                && $params->{"validation"} ne "relaxed"
                && $params->{"validation"} ne "basic");
    }
    else { # Set to strict validation automatically
        $params->{"validation"} = "strict";
    }
    
    return($params);
}

=head1 AUTHOR

Panagiotis Moulos, C<< <pmoulos at hybridstat.gr> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-vcf-io at rt.cpan.org>, or through
the web interface at L<https://rt.cpan.org/NoAuth/ReportBug.html?Queue=VCF-IO>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc VCF::IO


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<https://rt.cpan.org/NoAuth/Bugs.html?Dist=VCF-IO>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/VCF-IO>

=item * CPAN Ratings

L<https://cpanratings.perl.org/d/VCF-IO>

=item * Search CPAN

L<https://metacpan.org/release/VCF-IO>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2018 Panagiotis Moulos.

This program is distributed under the MIT (X11) License:
L<http://www.opensource.org/licenses/mit-license.php>

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.


=cut

1; # End of VCF::IO
