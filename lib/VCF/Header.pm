package VCF::Header;

use 5.10.0;
use strict;
use warnings;
use utf8;

use Carp;
use File::Basename;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Tie::IxHash::Easy;

use VCF::_Helper;

use vars qw($helper);

BEGIN {
    $helper = VCF::_Helper->new;
    binmode STDOUT, ":utf8";
    $|=1;
}

use Data::Dumper qw(Dumper);

=head1 NAME

VCF::Header - A Perl module to parse and manipulate VCF file headers and meta
information lines.

=head1 VERSION

Version 0.01

=cut

our $MODNAME = 'VCF::Header';
our $VERSION = '0.01';

=head1 SYNOPSIS

VCF::Header provides methods for parsing, validating, manipulating and writing
VCF file headers and meta-information lines based on the VCF file 
specifications. A VCF::Header object can be initialized by providing a VCF file
as a single input, or by providing a hash with parameters with the following 
posssible keys: C<file>, C<data>, C<validation>.

C<file> is a valid, existing VCF file, C<data> is a hash with header data. To
see an example, just parse a VCF file and inspect the has returned by the
function C<get_header>. If both C<file> and C<data> are provided, C<file> will
be used. C<validation> sets the header validation strictness.

There are three validation levels in VCF::Header:
 
=over

=item *

B<basic>: When validation is C<'basic'>, only very essential validation
procedures are performed, such as basic content checks (e.g. that the required
attributes do exist in the VCF header and meta-info lines).

=item *

B<relaxed>: When validation is C<'relaxed'>, additional checks are performed.
For example, while when validation is C<'basic'> the ##contig meta-info lines
may not exist. The relaxed validation may prove useful in cases of heavily 
customized applications which rely on the basics of the VCF format though.

=item *
 
B<'strict'>: When validation is C<'strict'>, VCF::Header performs as many 
validation steps as possible in order to comply with each VCF format version 
specs.
 
=back

If validation parameter is not provided, the strict validation is performed by
default.

Usage:

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

=head1 SUBROUTINES/METHODS

=head2 new

The VCF::Header object constructor. It accepts an existing VCF file.

    my $header = VCF::Header->new({
        'file' => 'my_variant.vcf'
    });

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

VCF::Header object simple OO parsing of a VCF file header and meta-info lines. 
This method does not perform any essential validations. The user may need to 
add a header element and validate next. If validation is desired, C<'parse'>
should be followed by the C<validate> method.

    $header->parse;
    # If separate validation 
    # $header->validate;

=cut

sub parse {
    my $self = $_[0];
    
    # Open file
    my $vcf = $self->get("file");

    open(VCF_FILE,$vcf) or croak "\nFile $vcf not found!";
    
    # Get first line, should be fileformat
    my $line = <VCF_FILE>;
    $line =~ s/\r|\n$//g;
    $line =~ s/^##//g;
    my @ff = split("=",$line);
    
    # An iconic fileformat must be defined (w/o parsing to header hash) just for
    # initiating specs
    $ff[1] = "VCFv4.3" if ($ff[0] ne "fileformat");
    
    # Init header scaffolds according to fileformat
    my $header = $self->_init_header_scaffold($ff[1]);
    my @spec_fields = @{$self->_init_header_metanames($ff[1])};
        
    # Get the rest
    $header->{"fileformat"} = $ff[1];
    while (my $line=<VCF_FILE>) {
        last if ($line !~ m/^#/);
        
        $line =~ s/\r|\n$//g;
        if ($line =~ m/^##/) { # Meta-fields
            $line =~ s/^##//g;
            my ($key,@all_values) = split("=",$line);
            # "=" exists internally...
            my $value = join("=",@all_values);
            
            # Start parsing header attributes
            if ($value =~ m/^</) {
                # Check closing tag and issue warning
                carp "\nWarning: unclosed VCF ".$key." header found in line $.! ".
                    "Please check the VCF file integrity. Offending value: ".
                    $value if ($value !~ m/>$/);
                # Continue...
                $value =~ s/^<|>$//g;
                $value = $self->_hack_comma_in_quotes($value);
                my %attr;
                tie %attr, "Tie::IxHash::Easy";
                my @attrs = split(",",$value);
                for (my $i=0; $i<@attrs; $i++) {
                    $attrs[$i] =~ s/^\s+|\s+$//g;
                    my @valpair = split("=",$attrs[$i]);
                    $valpair[0] =~ s/^\s+|\s+$//g;
                    # Extreme case of Description containing "=".....
                    my $vl = scalar @valpair;
                    if ($vl > 2) {
                        for (my $j=1; $j<$vl; $j++) {
                            $valpair[$j] =~ s/^\s+|\s+$//g;
                        }
                        $valpair[1] = join("=",@valpair[1..($vl-1)]);
                    }
                    else {
                        $valpair[1] =~ s/^\s+|\s+$//g;
                    }
                    $valpair[1] = "\"".
                        $self->_unhack_comma_in_quotes($valpair[1])."\""
                        if $valpair[0] eq "Description";
                    $attr{$valpair[0]} = $valpair[1];
                }
                (grep {$_ eq $key} @spec_fields) ? 
                    (push(@{$header->{$key}},\%attr)) :
                    (push(@{$header->{"other"}->{$key}},\%attr));
            }
            else {
                (grep {$_ eq $key} @spec_fields) ? ($header->{$key} = $value) :
                    ($header->{"other"}->{$key} = $value);
            }
        }
        else { # VCF columns header
            my @field_heads = split(/\t/,$line);
            $header->{"columns"} = \@field_heads;
        }
    }
    
    close(VCF_FILE);
    
    $self->set_header($header);
    
    return($self);
}

=head2 validate

VCF::Header object validation of a VCF file header and meta-info lines. This
method performs header field validations according to each supported file
format specifications and according to validation level. It can be used when
parsing and validating at separate stages, e.g. when a meta-info line must be
added in the VCF header.

    $header->validate;

=cut

sub validate {
    my $self = $_[0];
    
    croak "\nNo VCF header information found! Have you called parse?\n"
        if (!$self->_has_header);
    
    my %header = %{$self->get_header};
    
    # Get the validation level
    my $level = $self->get("validation");
    
    # First check if the absolutely required meta lines and column names
    # are present
    my @rec_fields = ("fileformat","INFO","FORMAT","FILTER");
    my @rec_cols = ("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO");
    
    # Check meta lines
    my @metas = keys(%header);
    foreach my $r (@rec_fields) {
        croak "\nIncomplete VCF header structure: required field $r not found!"
            if (!(grep {$_ eq $r} @metas));
    }
    
    # Then check columns names
    foreach my $r (@rec_cols) {
        croak "\nIncomplete VCF data column names: required column $r not found!"
            if (!(grep {$_ eq $r} @{$header{"columns"}}));
    }
    
    # But also check tha NO other than our allowed field names are present.
    my $ff = $header{"fileformat"};
    my @spec_fields = @{$self->_init_header_metanames($ff)};
    # Add other and columns so as not crash later
    push(@spec_fields,"other");
    push(@spec_fields,"columns");
    if ($level eq "strict") {
        foreach my $r (@metas) {
            croak "\nUnkown VCF header structure: field $r not in allowed ".
                "names!" if (!(grep {$_ eq $r} @spec_fields));
        }
    }

    # Now we have to validate each field
    # $value is an array of internal hashes
    while (my($key,$value) = each(%header)) {
        if ($key eq "INFO") {
            foreach my $info (@{$value}) {
                $self->_validate_header_field_info($info,$ff);
            }
        }
        if ($key eq "FILTER") {
            foreach my $filter (@{$value}) {
                $self->_validate_header_field_filter($filter,$ff);
            }
        }
        if ($key eq "FORMAT") {
            foreach my $format (@{$value}) {
                $self->_validate_header_field_format($format,$ff);
            }
        }
        if ($key eq "ALT") {
            foreach my $alt (@{$value}) {
                $value = $self->_validate_header_field_alt($alt,$ff);
            }
        }
        if ($key eq "SAMPLE") {
            foreach my $sample (@{$value}) {
                $self->_validate_header_field_sample($sample,$ff);
            }
        }
        if ($key eq "META") {
            foreach my $meta (@{$value}) {
                $self->_validate_header_field_meta($meta,$ff);
            }
        }
        # else, no, we can't validate "other"
    }
    
    return(0);
}

=head2 parse_and_validate

VCF::Header object parsing of a VCF file header and meta-info lines. This
method performs header field validations according to each supported file
format specifications and according to validation level.

    $header->parse_and_validate;

=cut

sub parse_and_validate {
    my $self = $_[0];
    
    # Get the validation level
    my $level = $self->get("validation");
    
    # Open file
    my $vcf = $self->get("file");
    open(VCF_FILE,$vcf) or die "\nFile $vcf not found!";
    
    # Get first line, must be fileformat
    my $line = <VCF_FILE>;
    $line =~ s/\r|\n$//g;
    $line =~ s/^##//g;
    my @ff = split("=",$line);
    croak "\nThe fileformat field was not found on the input file! Please ".
        "check if it's a proper VCF file.\n\n" if (lc($ff[0]) ne "fileformat");
    
    # Check several cases, more might be added. In each case we check its
    # validity according to the VCF specs. If some is missing, we should later
    # fill with our data from the database.
    
    # Init header scaffolds according to fileformat
    my $header = $self->_init_header_scaffold($ff[1]);
    my @spec_fields = @{$self->_init_header_metanames($ff[1])};
        
    my @rec_fields = ("INFO","FORMAT","FILTER");
    my @rec_cols = ("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO");
    
    # Get the rest
    $header->{"fileformat"} = $ff[1];
    while (my $line=<VCF_FILE>) {
        last if ($line !~ m/^#/);
        
        $line =~ s/\r|\n$//g;
        if ($line =~ m/^##/) { # Meta-fields
            $line =~ s/^##//g;
            my ($key,@all_values) = split("=",$line);
            # "=" exists internally...
            my $value = join("=",@all_values);
            
            # Continue with parsing
            if ($value =~ m/^</) {
                # Validate the fields for which a spec exists
                my $attr;
                if ($key eq "INFO") {
                    $attr = $self->_validate_header_field($value,"info",$ff[1],
                        $.);
                }
                elsif ($key eq "FILTER") {
                    $attr = $self->_validate_header_field($value,"filter",
                        $ff[1],$.);
                }
                elsif ($key eq "FORMAT") {
                    $attr = $self->_validate_header_field($value,"format",
                        $ff[1],$.);
                }
                elsif ($key eq "ALT") {
                    $attr = $self->_validate_header_field($value,"alt",
                        $ff[1],$.);
                }
                elsif ($key eq "SAMPLE") {
                    $attr = $self->_validate_header_field($value,"sample",
                        $ff[1],$.);
                }
                elsif ($key eq "META") {
                    $attr = $self->_validate_header_field($value,"meta",
                        $ff[1],$.);
                }
                else {
                    $attr = $self->_validate_header_field($value);
                }
                (grep {$_ eq $key} @spec_fields) ? 
                    (push(@{$header->{$key}},$attr)) :
                    (push(@{$header->{"other"}->{$key}},$attr));
            }
            else {
                (grep {$_ eq $key} @spec_fields) ? ($header->{$key} = $value) :
                    ($header->{"other"}->{$key} = $value);
            }
        }
        else { # VCF columns header
            my @field_heads = split(/\t/,$line);
            foreach my $r (@rec_cols) {
                croak "\nIncomplete VCF columns, line ".$..": required ".
                    "column $r not found!\nOffending line content: $line"
                        if (!(grep {$_ eq $r} @field_heads));
            }
            # If not dead, all required columns are here, so add them to the
            # header information
            $header->{"columns"} = \@field_heads;
        }
    }
    
    close(VCF_FILE);
    
    $self->set_header($header);
    
    return($self);
}

=head 2 parse_sample_names

This function parses a VCF file up to the point where #CHROM is found and from
this line it returns an array reference with the available samples in the VCF
file

$header->parse_sample_names;

=cut

sub parse_sample_names {
    my $self = $_[0];
    my $vcf = $self->get("file");
    
    my $sample_names;
    open(VCF_FILE,$vcf) or die "\nFile $vcf not found!";
    while (my $line = <VCF_FILE>) {
        if ($line =~ m/^#CHROM/) {
            $line =~ s/\r|\n$//g;
            $sample_names = $self->_parse_sample_names_internal($line);
            last;
        }
    }
    close(VCF_FILE);
    
    return($sample_names);
}

=head2 add

Validates a new meta-info line and adds it to the VCF header object. The first
argument is a hash with the meta-info line elements to be added. The second is
the type of meta-info in B<lowercase>.
    
    my $add_info = {
        'ID' => "MYANN",
        'Number' => ".",
        'Type' => "String",
        'Description' => "My annotation elements"
    };
    $header->add($add_info,"info");
    
    my $add_form = {
        'ID' => "MYFOR",
        'Number' => 1,
        'Type' => "Float",
        'Description' => "My format elements"
    };
    $header->add($add_form,"format");
    
    my $add_filt = {
        'ID' => "MYFIL",
        'Description' => "My filter QUAL>20 & DP>100"
    };
    $header->add($add_filt,"filter");

=cut

sub add {
    my ($self,$raw_attr,$type) = @_;
    my $header = $self->get_header;
    
    # Checking for already existing keys
    my $header_ids = $self->_extract_header_ids($header);
    my @ex;
    if ($helper->_is_hash_ref($header_ids->{$type})) {
        @ex = keys(%{$header_ids->{$type}});
    }
    elsif ($helper->_is_array_ref($header_ids->{$type})) {
        @ex = @{$header_ids->{$type}};
    }
    my $k = $raw_attr->{"ID"};
    die "You are trying to add $type attribute $k which already exists! ".
        "Please use the change function!" if (grep{$_ eq $k} @ex);
    
    # Validate using validation subroutines
    my $attr;
    my $vcfv = $header->{"fileformat"};
    if (lc($type) eq "info") {
        $attr = $self->_prepare_new_header_field_info($raw_attr);
        $self->_validate_header_field_info($attr,$vcfv);
    }
    elsif (lc($type) eq "format") {
        $attr = $self->_prepare_new_header_field_format($raw_attr);
        $self->_validate_header_field_format($attr,$vcfv);
    }
    elsif (lc($type) eq "filter") {
        $attr = $self->_prepare_new_header_field_filter($raw_attr);
        $self->_validate_header_field_filter($attr,$vcfv);
    }
    elsif (lc($type) eq "alt") {
        $attr = $self->_prepare_new_header_field_alt($raw_attr);
        $self->_validate_header_field_alt($attr,$vcfv);
    }
    elsif (lc($type) eq "sample") {
        $attr = $self->_prepare_new_header_field_sample($raw_attr);
        $self->_validate_header_field_sample($attr,$vcfv);
    }
    elsif (lc($type) eq "meta") {
        $attr = $self->_prepare_new_header_field_meta($raw_attr);
        $self->_validate_header_field_meta($attr,$vcfv);
    }
    
    # If not dead, add to the header
    # FIXME: Check naming robustness,e.g. INFO - INFo, fileDate - fileDate
    #        Do we allow this?
    my $metanames = $self->_init_header_metanames($vcfv);
    my @array_meta = ("FILTER","INFO","FORMAT","ALT","SAMPLE","META","PEDIGREE",
        "contig");
    my %in_array_meta = map {$_ => 1} @array_meta;
    my @not_array_types = grep {not $in_array_meta{$_}} @{$metanames};
    if (grep {$_ eq $type} @array_meta) {
        push(@{$header->{$type}},$attr);
    }
    elsif (grep {$_ eq $type} @not_array_types) {
        $header->{$type} = $attr;
    }
    else {
        push(@{$header->{"other"}},$attr);
    }
    
    $self->set_header($header);
    
    return($self);
}

=head2 change

Change the attributes of a current meta-information. This function does B<not>
validate the VCF records if the change of a header element breaks a certain 
data type. This should be done by the respective VCF::Record facilities. Also,
this function will B<not> add new keys and values to an existing meta-info line.
For example if an I<INFO> record has I<ID>, I<Number>, I<Type> and 
I<Description> keys with respective values, a new potential key I<MyDescription>
will not be added as it may violate the VCF I<fileformat> attribute. Finally,
the I<fileformat> attribute cannot be changed.

The C<change> method must generally be used with much caution.
    
    my $with = {
        'id' => 'DP',
        'key' => 'Description',
        'value' => 'New description'
    };
    $header->change("INFO",$with);
    
    my $with = "20170630"
    $header->change("fileDate",$with);

=cut

# $with must be either:
# - a scalar in the case of simple meta-info lines (e.g. fileformat)
# - a hash with ID, $key, $value elements, e.g.
# $with = {
#     'id' => 'DP',
#     'key' => 'Description',
#     'value' => 'New description'
# };

sub change {
    my ($self,$what,$with) = @_;
    my $header = $self->get_header;
    
    if ($helper->_is_hash_ref($with)) {
        # It should be one of "FILTER","INFO","FORMAT","ALT","SAMPLE","META",
        # "PEDIGREE","contig"
        my @array_multi = ("FILTER","INFO","FORMAT","ALT","SAMPLE","META",
            "PEDIGREE","contig");
        my @existing = keys(%{$header});
        my %array_multi_map = map {$_ => 1} @array_multi;
        my @common = grep($array_multi_map{$_},@existing);
        if (!$helper->_smatch($what,@common)) {
            carp "\n$what attribute not found in the header! Returning...\n";
            return(0);
        }
        
        # $with must have id, key and value keys
        my @require_with = ("id","key","value");
        my @newk = keys(%{$with});
        foreach my $r (@require_with) {
            croak "\nRequired key $r for changing header attribute not found!\n"
                if (!$helper->_smatch($r,@newk));
        }
        
        # If all ok do the change work
        my $id = $with->{"id"};
        my $counter = 0;
        foreach my $h (@{$header->{$what}}) {
            if ($h->{"ID"} eq $with->{"id"}) {
                $header->{$what}->[$counter]->{$with->{"key"}} = 
                    $with->{"value"};
                last;
            }
            $counter++;
        }
    }
    elsif (!ref($with) && defined($with)) { # Should be scalar
        $header->{$what} = $with;
    }
    
    # Now, revalidate the header, if something misplaced, the process will die
    $self->set_header($header);
    $self->validate;
    
    # If all OK, return
    return($self);
}

=head2 write

VCF::Header writing to file. The C<output> method must be called first or a
valid filehandle must be provided.

    $header->write;

=cut

sub write {
    my ($self,$fh) = @_;
    
    # Even if $self->output has been called with a file handle and a 
    # file handle is provided, the second precedes. If $self->output has not
    # been called, and $fh is not a valid file handle output will be written
    # to STDOUT
    $self->{'_output'} = $fh if ($fh && fileno($fh) != -1);
    
    my $header = $self->get_header;
    my @header_lines = @{$self->_prepare_header_lines};
    my $column_names = 
        $self->_format_column_names(@{$self->_get_column_names});
    push(@header_lines,$column_names);
    
    # If file handle given, write to file otherwise to STDOUT
    my $fout = $self->{'_output'};
    ($fout) ? print $fout join("\n",@header_lines),"\n" :
        print STDOUT join("\n",@header_lines),"\n";
}

=head2 output

Opens a file for writing the parsed/manipulated VCF header. Essentially a wraper
over Perl's C<open>. The second argument determines the writing mode, which is
one of C<'new'> (default) or C<'append'>.

    $header->output('header.vcf');

=cut

sub output {
    my ($self,$output,$mode) = @_;
    
    croak "\nNo file provided for writing!\n" if (!$output);
    if ($mode) {
        croak "\nWriting mode must be either 'new' or 'append'!\n"
            if ($mode ne "new" && $mode ne "append");
        $mode = $mode eq "new" ? ">" : ">>";
    }
    else { # New/overwrite by default
        $mode = ">";
    }

    open(my $fh,$mode,$output) or croak "Could not open $output for writing!";
    $self->{'_output'} = $fh;
    return($self);
}

=head2 done

Closes the file where VCF header is written. Essentially a  wraper over Perl's
C<close>.

    $header->done;

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

=head2 get_header

Returns a simple hash with the VCF header fields.

    my $header = $header->get_header;
    
TODO: Add documentation (the field names)

=cut

sub get_header {
    my $self = $_[0];
    return($self->{"_content"});
}

=head2 get

VCF::Header object getter

    my $param_value = $header->get("param_name");
    
TODO: Validate that the parameter to be got is valid

=cut

sub get {
    my ($self,$name) = @_;
    return($self->{$name});
}

=head2 set_header

Sets the actual parsed header content of the VCF file (a hash) either 
after parsing or constructed manually.
    
    my $file = 'my_variants.vcf';
    my $my_header = &parse_vcf_header($file);
    $header->set_header($my_header);
    
=cut

sub set_header {
    my ($self,$header) = @_;
    $self->{"_content"} = $header;
    return($self);
}

=head2 set

VCF::Header object setter

    $header->set("param_name","param_value");
    
TODO: Validate that the parameter to be set is valid
    
=cut

sub set {
    my ($self,$name,$value) = @_;
    $self->{$name} = $value;
    return($self);
}

# Internal methods

sub _validate_header_field {
    my ($self,$value,$type,$vcfv,$ln) = @_;
    
    # Get validation level
    my $level = $self->get("validation");
    
    # Check closing tag and issue warning or die according to validation
    if ($level eq "basic") {
        carp "\nWarning: unclosed VCF ".uc($type)." header found in line $ln! ".
            "Please check the VCF file integrity. Offending value: $value" 
                if ($value !~ m/>$/);
    }
    elsif ($level eq "relaxed" || $level eq "strict") {
        croak "\nWarning: unclosed VCF ".uc($type)." header found in line ".
            "$ln! Please check the VCF file integrity. Offending value: $value" 
                if ($value !~ m/>$/);
    }
    
    $value =~ s/^<|>$//g;
    $value = $self->_hack_comma_in_quotes($value);
    
    my %attr;
    tie %attr, "Tie::IxHash::Easy";
    my @attrs = split(",",$value);
    for (my $i=0; $i<@attrs; $i++) {
        $attrs[$i] =~ s/^\s+|\s+$//g;
        my @valpair = split("=",$attrs[$i]);
        $valpair[0] =~ s/^\s+|\s+$//g;
        # Extreme case of Description containing "=".....
        my $vl = scalar @valpair;
        if ($vl > 2) {
            for (my $j=1; $j<$vl; $j++) {
                $valpair[$j] =~ s/^\s+|\s+$//g;
            }
            $valpair[1] = join("=",@valpair[1..($vl-1)]);
        }
        else {
            $valpair[1] =~ s/^\s+|\s+$//g;
        }
        $valpair[1] = "\"".$self->_unhack_comma_in_quotes($valpair[1])."\""
            if $valpair[0] eq "Description";
        $attr{$valpair[0]} = $valpair[1];
    }
    
    if ($type) {
        if ($type eq "info") {
            $self->_validate_header_field_info(\%attr,$vcfv,$value,$ln);
        }
        elsif ($type eq "format") {
            $self->_validate_header_field_format(\%attr,$vcfv,$value,$ln);
        }
        elsif ($type eq "filter") {
            $self->_validate_header_field_filter(\%attr,$vcfv,$value,$ln);
        }
        elsif ($type eq "alt") {
            $self->_validate_header_field_alt(\%attr,$vcfv,$value,$ln);
        }
        elsif ($type eq "sample") {
            $self->_validate_header_field_sample(\%attr,$vcfv,$value,$ln);
        }
        elsif ($type eq "meta") {
            $self->_validate_header_field_meta(\%attr,$vcfv,$value,$ln);
        }
    }
    
    return(\%attr);
}

# Two input possibilities:
# 1. $value and $ln provided: coming from file validation
# 2. $value and $ln not provided: coming from internal/API construction
sub _validate_header_field_info {
    my $self = $_[0];
    my %attr = %{$_[1]};
    my $vcfv = $_[2];
    my $value = $_[3];
    my $ln = $_[4];
    
    # Get the validation level
    my $level = $self->get("validation");
    
    my @req = ("ID","Number","Type","Description");
    my @types = ("Integer","Float","Flag","Character","String");
    
    my @hav = keys(%attr);
    # Invalid cases:
    my %check = %{$helper->_unique(@hav)};
    
    # Validate generic fields
    # %attr coming from file parsing, $value and $ln given
    if ($value && $ln) {
        # Validation always
        # 1. hav not in req
        foreach my $r (@req) {
            croak "\nCorrupted VCF INFO header, line ".$ln.": required key ".$r.
                " not found!\nOffending line content: $value"
                    if (!defined($check{$r}) || !$check{$r});
        }
        # 2. missing values in hav
        while (my ($k,$v) = each %attr) {
            croak "\nCorrupted VCF INFO header, line ".$ln.": key ".$k.
                " does not have a value!\nOffending line content: $value"
                    if (!defined($v) || length($v)==0);
        }
        
        # Validation in relaxed and strict
        if ($level eq "strict" || $level eq "relaxed") {
            # 3. Type not in ("Integer","Float","Flag","Character","String")
            croak "\nCorrupted VCF INFO header, line ".$ln.": key Type not ".
                "one of (Integer,Float,Flag,Character,String)!\nOffending ".
                "line content: $value"
                    if (!(grep {$_ eq $attr{"Type"}} @types));
                    
            # 4. Type is "Flag" but Number is not 0
            croak "\nCorrupted VCF INFO header, line ".$ln.": key Number is ".
                "not 0 when Type is Flag!\nOffending line content: $value"
                    if ($attr{"Type"} eq "Flag" && $attr{"Number"} ne '0');
        }
        
        # Validation in strict only
        if ($level eq "strict") {
            # 5. duplicates in hav
            while (my ($k,$v) = each %check) {
                croak "\nCorrupted VCF INFO header, line ".$ln.": key ".$k.
                    " found multiple times!\nOffending line content: $value"
                        if ($v > 1);
            }
            
            # 6. Description is not a string
            croak "\nCorrupted VCF INFO header, line ".$ln.": key Description ".
                "not a string!\nOffending line content: $value"
                    if (!defined($attr{"Description"})
                        || ref($attr{"Description"}));
        }

        # Validate version specific fields
        ## $expr_{n,i} defined two times because potentially coming from
        ## different VCF file formats
        my $expr_n = $self->_init_header_expr($vcfv);
        my $number = $expr_n->{"INFO"}->{"Number"};
        # 1. Number is not 0 or a positive integer or not in @numbers
        if ($level eq "strict" || $level eq "relaxed") {
            croak "\nCorrupted VCF INFO header, line ".$ln.": key Number not ".
                "an integer or a valid character!\nOffending line ".
                    "content: $value" if ($attr{"Number"} !~ m/$number/);
        }
        # 2. ID must follow a specific expression
        my $expr_i = $self->_init_header_expr($vcfv);
        if ($level eq "strict") {
            if (defined($expr_i->{"INFO"}->{"ID"})) {
                my $id = $expr_i->{"INFO"}->{"ID"};
                croak "\nCorrupted VCF INFO header, line ".$ln.": key ID not ".
                    "following naming specifications!\nOffending line ".
                        "content: ".$value if ($attr{"ID"} !~ m/$id/);
            }
        }
    }
    else { # Not given, adjust messages
        # Validation always
        # 1. hav not in req
        foreach my $r (@req) {
            croak "\nCorrupted VCF INFO header: required key $r not found!"
                if (!defined($check{$r}) || !$check{$r});
        }
        # 2. missing values in hav
        while (my ($k,$v) = each %attr) {
            croak "\nCorrupted VCF INFO header: key $k does not have a value!"
                if (!defined($v) || length($v)==0);
        }
        
        # Validation in relaxed and strict
        if ($level eq "strict" || $level eq "relaxed") {
            # 3. Type not in ("Integer","Float","Flag","Character","String")
            croak "\nCorrupted VCF INFO header: key Type not one of (Integer,".
                "Float,Flag,Character,String)!" 
                    if (!(grep {$_ eq $attr{"Type"}} @types));
                    
            # 4. Type is "Flag" but Number is not 0
            croak "\nCorrupted VCF INFO header: key Number is not 0 when Type ".
                "is Flag!" 
                    if ($attr{"Type"} eq "Flag" && $attr{"Number"} ne '0');
        }
        
        # Validation in strict only
        if ($level eq "strict") {
            # 5. duplicates in hav
            while (my ($k,$v) = each %check) {
                croak "\nCorrupted VCF INFO header: key $k found multiple ".
                    "times!" if ($v > 1);
            }
            
            # 6. Description is not a string
            croak "\nCorrupted VCF INFO header: key Description not a string!"
                if (!defined($attr{"Description"}) 
                    || ref($attr{"Description"}));
        }

        # Validate version specific fields
        my $expr_n = $self->_init_header_expr($vcfv);
        my $number = $expr_n->{"INFO"}->{"Number"};
        # 1. Number is not 0 or a positive integer or not in @numbers
        if ($level eq "strict" || $level eq "relaxed") {
            croak "\nCorrupted VCF INFO header: key Number not an integer or ".
                "a valid character!" if ($attr{"Number"} !~ m/$number/);
        }
        # 2. ID must follow a specific expression
        my $expr_i = $self->_init_header_expr($vcfv);
        if ($level eq "strict") {
            if (defined($expr_i->{"INFO"}->{"ID"})) {
                my $id = $expr_i->{"INFO"}->{"ID"};
                croak "\nCorrupted VCF INFO header: key ID not following ".
                    "naming specifications!" if ($attr{"ID"} !~ m/$id/);
            }
        }
    }
    
    return(0);
}

sub _validate_header_field_format {
    my $self = $_[0];
    my %attr = %{$_[1]};
    my $vcfv = $_[2];
    my $value = $_[3];
    my $ln = $_[4];
    
    # Get the validation level
    my $level = $self->get("validation");
    
    my @req = ("ID","Number","Type","Description");
    my @types = ("Integer","Float","Character","String");
    
    my @hav = keys(%attr);
    # Invalid cases:
    my %check = %{$helper->_unique(@hav)};
    
    # Validate generic fields
    if ($value && $ln) { # From file
        # Validation always
        # 1. hav not in req
        foreach my $r (@req) {
            croak "\nCorrupted VCF FORMAT header, line ".$ln.": required ".
                "key $r not found!\nOffending line content: $value"
                    if (!defined($check{$r}) || !$check{$r});
        }
        # 2. missing values in hav
        while (my ($k,$v) = each %attr) {
            croak "\nCorrupted VCF FORMAT header, line ".$ln.": key ".$k.
                " does not have a value!\nOffending line content: $value"
                    if (!defined($v) || length($v)==0);
        }
        
        # Validation in relaxed and strict
        if ($level eq "strict" || $level eq "relaxed") {
            # 3. Type not in ("Integer","Float","Character","String")
            croak "\nCorrupted VCF FORMAT header, line ".$ln.": key Type not ".
                "one of (Integer,Float,Character,String)!\nOffending line ".
                "content: ".$value if (!(grep {$_ eq $attr{"Type"}} @types));
        }
        
        # Validation in strict only
        if ($level eq "strict") {
            # 4. duplicates in hav
            while (my ($k,$v) = each %check) {
                croak "\nCorrupted VCF FORMAT header, line ".$ln.": key ".$k.
                    " found multiple times!\nOffending line content: $value"
                        if ($v > 1);
            }
            
            # 5. Description is not a string
            croak "\nCorrupted VCF FORMAT header, line ".$ln.": key ".
                "Description not a string!\nOffending line content: $value"
                    if (!defined($attr{"Description"})
                        || ref($attr{"Description"}));
        }
        
                
        # Validate version specific fields
        my $expr = $self->_init_header_expr($vcfv);
        my $number = $expr->{"FORMAT"}->{"Number"};
        # 1. Number is a positive integer or not in @numbers
        if ($level eq "strict" || $level eq "relaxed") {
            croak "\nCorrupted VCF FORMAT header, line ".$ln.": key Number ".
            "not an integer or a valid character!\nOffending line ".
            "content: $value " if ($attr{"Number"} !~ m/$number/);
        }
    }
    else { # From structure
        # Validation always
        # 1. hav not in req
        foreach my $r (@req) {
            croak "\nCorrupted VCF FORMAT header: required key $r not found!"
                if (!defined($check{$r}) || !$check{$r});
        }
        # 2. missing values in hav
        while (my ($k,$v) = each %attr) {
            croak "\nCorrupted VCF FORMAT header: key $k does not have a value!"
                if (!defined($v) || length($v)==0);
        }
        
        # Validation in relaxed and strict
        if ($level eq "strict" || $level eq "relaxed") {
            # 3. Type not in ("Integer","Float","Character","String")
            croak "\nCorrupted VCF FORMAT header: key Type not one of ".
                "(Integer,Float,Character,String)!"
                    if (!(grep {$_ eq $attr{"Type"}} @types));
        }
        
        # Validation in strict only
        if ($level eq "strict") {
            # 4. duplicates in hav
            while (my ($k,$v) = each %check) {
                croak "\nCorrupted VCF FORMAT header: key $k found multiple ".
                    "times!" if ($v > 1);
            }

            # 5. Description is not a string
            croak "\nCorrupted VCF FORMAT header: key Description not a string!"
                if (!defined($attr{"Description"}) 
                    || ref($attr{"Description"}));
        }
                            
        # Validate version specific fields
        my $expr = $self->_init_header_expr($vcfv);
        my $number = $expr->{"FORMAT"}->{"Number"};
        # 1. Number is a positive integer or not in @numbers
        if ($level eq "strict" || $level eq "relaxed") {
            croak "\nCorrupted VCF FORMAT header: key Number not an integer ".
                "or a valid character!" if ($attr{"Number"} !~ m/$number/);
        }
    }
    
    return(0);
}

sub _validate_header_field_filter {
    my $self = $_[0];
    my %attr = %{$_[1]};
    my $vcfv = $_[2];
    my $value = $_[3];
    my $ln = $_[4];
    
    # Get the validation level
    my $level = $self->get("validation");
    
    my @req = ("ID","Description");
    
    my @hav = keys(%attr);
    # Invalid cases:
    my %check = %{$helper->_unique(@hav)};
    
    if ($value && $ln) { # From file
        # Validation always
        # 1. hav not in req
        foreach my $r (@req) {
            croak "\nCorrupted VCF FILTER header, line ".$ln.": required ".
                "key $r not found!\nOffending line content: $value"
                    if (!defined($check{$r}) || !$check{$r});
        }
        # 2. missing values in hav
        while (my ($k,$v) = each %attr) {
            croak "\nCorrupted VCF FILTER header, line ".$ln.": key ".$k.
                " does not have a value!\nOffending line content: $value"
                    if (!defined($v) || length($v)==0);
        }
        
        # Validation in strict only
        if ($level eq "strict") {
            # 3. duplicates in hav
            while (my ($k,$v) = each %check) {
                croak "\nCorrupted VCF FILTER header, line ".$ln.": key ".$k.
                    " found multiple times!\nOffending line content: $value"
                        if ($v > 1);
            }
            
            # 4. Description is not a string
            croak "\nCorrupted VCF FILTER header, line ".$ln.": key ".
                "Description not a string!\nOffending line content: $value"
                    if (!defined($attr{"Description"})
                        || ref($attr{"Description"}));
        }
    }
    else { # From structure
        # Validation always
        # 1. hav not in req
        foreach my $r (@req) {
            croak "\nCorrupted VCF FILTER header: required key $r not found!"
                if (!defined($check{$r}) || !$check{$r});
        }
        # 2. missing values in hav
        while (my ($k,$v) = each %attr) {
            croak "\nCorrupted VCF FILTER header: key $k does not have a value!"
                if (!defined($v) || length($v)==0);
        }
        
        # Validation in strict only
        if ($level eq "strict") {
            # 3. duplicates in hav
            while (my ($k,$v) = each %check) {
                croak "\nCorrupted VCF FILTER header: key $k found multiple ".
                    "times!" if ($v > 1);
            }
        
            # 4. Description is not a string
            croak "\nCorrupted VCF FILTER header: key Description not a string!"
                if (!defined($attr{"Description"}) 
                    || ref($attr{"Description"}));
        }
    }
    
    return(0);
}

sub _validate_header_field_alt {
    my $self = $_[0];
    my %attr = %{$_[1]};
    my $vcfv = $_[2];
    my $value = $_[3];
    my $ln = $_[4];
    
    # Get the validation level
    my $level = $self->get("validation");
    
    my @req = ("ID","Description");
    my @types = ("DEL","INS","DUP","INV","CNV","DUP:TANDEM","DEL:ME","INS:ME");
    
    my @hav = keys(%attr);
    # Invalid cases:
    my %check = %{$helper->_unique(@hav)};
    
    if ($value && $ln) { # From file
        # Validation always
        # 1. hav not in req
        foreach my $r (@req) {
            croak "\nCorrupted VCF ALT header, line ".$ln.": required key ".$r.
                " not found!\nOffending line content: $value"
                    if (!defined($check{$r}) || !$check{$r});
        }
        # 2. missing values in hav
        while (my ($k,$v) = each %attr) {
            croak "\nCorrupted VCF ALT header, line ".$ln.": key ".$k.
                " does not have a value!\nOffending line content: $value"
                    if (!defined($v) || length($v)==0);
        }
        
        # Validation in relaxed and strict
        if ($level eq "strict" || $level eq "relaxed") {
            # FIXME: Cannot use this until we find a vcf file example with 
            # this info
            ## 3. Type not in ("DEL","INS","DUP","INV","CNV","DUP:TANDEM",
            ## "DEL:ME","INS:ME");
            #croak "\nCorrupted VCF ALT header, line ".$ln.": key Type not one ".
            #   "(Integer,Float,Character,String)!\nOffending line "
            #    "content: $value" if (!(grep {$_ eq $attr{"Type"}} @types));
        }
        
        # Validation in strict only
        if ($level eq "strict") {
            # 4. duplicates in hav
            while (my ($k,$v) = each %check) {
                croak "\nCorrupted VCF ALT header, line ".$ln.": key ".$k.
                    " found multiple times!\nOffending line content: $value"
                        if ($v > 1);
            }
            
            # 5. Description is not a string
            croak "\nCorrupted VCF ALT header, line ".$ln.": key Description ".
                "not a string!\nOffending line content: $value"
                    if (!defined($attr{"Description"})
                        || ref($attr{"Description"}));
        }
    }
    else { # From structure
        # Validation always
        # 1. hav not in req
        foreach my $r (@req) {
            croak "\nCorrupted VCF ALT header: required key $r not found!"
                if (!defined($check{$r}) || !$check{$r});
        }
        # 2. missing values in hav
        while (my ($k,$v) = each %attr) {
            croak "\nCorrupted VCF ALT header: key $k does not have a value!"
                if (!defined($v) || length($v)==0);
        }
        
        # Validation in relaxed and strict
        if ($level eq "strict" || $level eq "relaxed") {
            # FIXME: Cannot use this until we find a vcf file example with
            # this info
            ## 3. Type not in ("DEL","INS","DUP","INV","CNV","DUP:TANDEM",
            ## "DEL:ME","INS:ME");
            #croak "\nCorrupted VCF ALT header: key Type not one of (Integer,".
            #    "Float,Character,String)!" 
            #        if (!(grep {$_ eq $attr{"Type"}} @types));
        }
        
        # Validation in strict only
        if ($level eq "strict") {
            # 4. duplicates in hav
            while (my ($k,$v) = each %check) {
                croak "\nCorrupted VCF ALT header: key $k found multiple times!"
                    if ($v > 1);
            }
            
            # 5. Description is not a string
            croak "\nCorrupted VCF ALT header: key Description not a string!"
                if (!defined($attr{"Description"}) 
                    || ref($attr{"Description"}));
        }
    }
    
    return(0);
}

# FIXME: More validation is needed here. Another function validating META and
#        SAMPLE key pairs. If version <4.3 then just parse.
sub _validate_header_field_sample {
    my $self = $_[0];
    my %attr = %{$_[1]};
    my $vcfv = $_[2];
    my $value = $_[3];
    my $ln = $_[4];
    
    # Get the validation level
    my $level = $self->get("validation");
    
    my @req = ("ID","Description");
    
    my @hav = keys(%attr);
    # Invalid cases:
    my %check = %{$helper->_unique(@hav)};
    
    if ($value && $ln) { # From file
        # Validation always
        # 1. hav not in req
        foreach my $r (@req) {
            croak "\nCorrupted VCF SAMPLE header, line ".$ln.": required ".
                "key $r not found!\nOffending line content: $value"
                    if (!defined($check{$r}) || !$check{$r});
        }
        # 2. missing values in hav
        while (my ($k,$v) = each %attr) {
            croak "\nCorrupted VCF SAMPLE header, line ".$ln.": key ".$k.
                " does not have a value!\nOffending line content: $value"
                    if (!defined($v) || length($v)==0);
        }
        
        # Validation in strict only
        if ($level eq "strict") {
            # 3. duplicates in hav
            while (my ($k,$v) = each %check) {
                croak "\nCorrupted VCF SAMPLE header, line ".$ln.": key ".$k.
                    " found multiple times!\nOffending line content: $value"
                        if ($v > 1);
            }
            
            # 4. Description is not a string
            croak "\nCorrupted VCF SAMPLE header, line ".$ln.": key ".
                "Description not a string!\nOffending line content: $value"
                    if (!defined($attr{"Description"})
                        || ref($attr{"Description"}));
        }
    }
    else { # From structure
        # Validation always
        # 1. hav not in req
        foreach my $r (@req) {
            croak "\nCorrupted VCF SAMPLE header : required key $r not found!"
                if (!defined($check{$r}) || !$check{$r});
        }       
        # 2. missing values in hav
        while (my ($k,$v) = each %attr) {
            croak "\nCorrupted VCF SAMPLE header key $k does not have a value!"
                if (!defined($v) || length($v)==0);
        }
        
        # Validation in strict only
        if ($level eq "strict") {
            # 3. duplicates in hav
            while (my ($k,$v) = each %check) {
                croak "\nCorrupted VCF SAMPLE header: key $k found multiple ".
                    "times!" if ($v > 1);
            }
            
            # 4. Description is not a string
            croak "\nCorrupted VCF SAMPLE header: key Description not a string!"
                if (!defined($attr{"Description"}) 
                    || ref($attr{"Description"}));
        }
    }
    
    return(0);
}

# TODO: Validate basic, another function should check the keys, requiring to be
#       present in both samples and meta sample info
sub _validate_header_field_meta {
    my $self = $_[0];
    my %attr = %{$_[0]};
    my $vcfv = $_[1];
    my $value = $_[2];
    my $ln = $_[3];
    
    # Get the validation level
    my $level = $self->get("validation");
    
    my @req = ("ID","Type","Number","Values");
    my @types = ("Integer","Float","Character","String");
    
    my @hav = keys(%attr);
    # Invalid cases:
    my %check = %{$helper->_unique(@hav)};
    
    if ($value && $ln) { # From file
        # Validation always
        # 1. hav not in req
        foreach my $r (@req) {
            croak "\nCorrupted VCF META header, line ".$ln.": required key ".$r.
                " not found!\nOffending line content: $value"
                    if (!defined($check{$r}) || !$check{$r});
        }   
        # 2. missing values in hav
        while (my ($k,$v) = each %attr) {
            croak "\nCorrupted VCF META header, line ".$ln.": key ".$k.
                " does not have a value!\nOffending line content: $value"
                    if (!defined($v) || length($v)==0);
        }
        
        # Validation in relaxed and strict
        if ($level eq "strict" || $level eq "relaxed") {
            # 3. Type not in ("Integer","Float","Character","String")
            croak "\nCorrupted VCF META header, line ".$ln.": key Type not ".
                "one of (Integer,Float,Character,String)!\nOffending line ".
                    "content: $value" if (!(grep {$_ eq $attr{"Type"}} @types));
        }
        
        # Validation in strict only
        if ($level eq "strict") {
            # 4. duplicates in hav
            while (my ($k,$v) = each %check) {
                croak "\nCorrupted VCF META header, line ".$ln.": key ".$k.
                    " found multiple times!\nOffending line content: $value"
                        if ($v > 1);
            }
        
            # 5. Description is not a string
            croak "\nCorrupted VCF META header, line ".$ln.": key Description ".
                "not a string!\nOffending line content: $value"
                    if (!defined($attr{"Description"})
                        || ref($attr{"Description"}));
        }
                   
        # Validate version specific fields
        my $expr = $self->_init_header_expr($vcfv);
        my $number = $expr->{"META"}->{"Number"};
        # 1. Number is not 0 or a positive integer or not in @numbers
        if ($level eq "strict" || $level eq "relaxed") {
            croak "\nCorrupted VCF META header, line ".$ln.": key Number not ".
                "one an integer or a valid character!\nOffending line ".
                "content: $value" if ($attr{"Number"} !~ m/$number/);
        }
    }
    else { # From structure
        # Validation always
        # 1. hav not in req
        foreach my $r (@req) {
            croak "\nCorrupted VCF META header: required key $r not found!"
                if (!defined($check{$r}) || !$check{$r});
        }   
        # 2. missing values in hav
        while (my ($k,$v) = each %attr) {
            croak "\nCorrupted VCF META header: key $k does not have a value!"
                if (!defined($v) || length($v)==0);
        }
        
        # Validation in relaxed and strict
        if ($level eq "strict" || $level eq "relaxed") {
            # 3. Type not in ("Integer","Float","Character","String")
            croak "\nCorrupted VCF META header: key Type not one of (Integer,".
                "Float,Character,String)!" 
                    if (!(grep {$_ eq $attr{"Type"}} @types));
        }
        
        # Validation in strict only
        if ($level eq "strict") {
            # 4. duplicates in hav
            while (my ($k,$v) = each %check) {
                croak "\nCorrupted VCF META header: key $k found multiple ".
                    "times!" if ($v > 1);
            }
        
            # 5. Description is not a string
            croak "\nCorrupted VCF META header: key Description not a string!"
                if (!defined($attr{"Description"}) 
                    || ref($attr{"Description"}));
        }
                    
        # Validate version specific fields
        my $expr = $self->_init_header_expr($vcfv);
        my $number = $expr->{"META"}->{"Number"};
        # 1. Number is not 0 or a positive integer or not in @numbers
        if ($level eq "strict" || $level eq "relaxed") {
            croak "\nCorrupted VCF META header: key Number not an integer or ".
                "a valid character!" if ($attr{"Number"} !~ m/$number/);
        }
    }
    
    return(0);
}

sub _extract_header_ids {
    my $self = $_[0];
    my $header = $self->get_header;
    
    # We suppose that header is validated, so we don't any further check and at
    # least the required fields exist.
    
    # Initialize a hash with the header fields that may have ID attribute
    my $header_ids = {
        "contig" => {},
        "INFO" => {},
        "FORMAT" => {},
        "FILTER" => [],
        "ALT" => [],
        "SAMPLE" => []
    };
    
    # VCF version must always be there
    $header_ids->{"fileformat"} = $header->{"fileformat"};
    
    # Check if contigs exist so we can use them later for validation
    if (scalar @{$header->{"contig"}} > 0) {
        foreach my $contig (@{$header->{"contig"}}) {
            $header_ids->{"contig"}->{$contig->{"ID"}} = 
                int $contig->{"length"};
        }
    }
    
    # Check if INFO exists so we can use them later for validation
    if (scalar @{$header->{"INFO"}} > 0) {
        foreach my $info (@{$header->{"INFO"}}) {
            $header_ids->{"INFO"}->{$info->{"ID"}} = {
                "Number" => $info->{"Number"},
                "Type" => $info->{"Type"}
            }
        }
    }
    
    # Check if FORMAT exists so we can use them later for validation
    if (scalar @{$header->{"FORMAT"}} > 0) {
        foreach my $format (@{$header->{"FORMAT"}}) {
            $header_ids->{"FORMAT"}->{$format->{"ID"}} = {
                "Number" => $format->{"Number"},
                "Type" => $format->{"Type"}
            }
        }
    }
    
    # Check if FILTER exists so we can use them later for validation
    if (scalar @{$header->{"FILTER"}} > 0) {
        foreach my $filter (@{$header->{"FILTER"}}) {
            push(@{$header_ids->{"FILTER"}},$filter->{"ID"});
        }
    }
    
    # Check if ALT exists so we can use them later for validation
    if (scalar @{$header->{"ALT"}} > 0) {
        foreach my $alt (@{$header->{"ALT"}}) {
            push(@{$header_ids->{"ALT"}},$alt->{"ID"});
        }
    }
    
    # Check if SAMPLE exists so we can use them later for validation
    if (scalar @{$header->{"SAMPLE"}} > 0) {
        foreach my $sample (@{$header->{"SAMPLE"}}) {
            push(@{$header_ids->{"SAMPLE"}},$sample->{"ID"});
        }
    }
    
    return($header_ids);
}

sub _extract_sample_names {
    my $self = $_[0];
    my $header = $self->get_header;
    my @a = @{$header->{"columns"}};
    return([]) if ($#a eq 7);
    return([$a[9]]) if ($#a eq 9);
    return(\@a[9..$#a]);
}

sub _prepare_header_lines {
    my $self = $_[0];
    
    my $header = $self->get_header;
    my @header_lines;
    my $line;
    
    # Write fileformat
    $line = $self->_format_header_lines("fileformat",$header);
    push(@header_lines,$line->[0]);
    
    # Write fileDate if present
    if ($header->{"fileDate"}) {
        $line = $self->_format_header_lines("fileDate",$header);
        push(@header_lines,$line->[0]);
    }
    
    # Write source if present
    if ($header->{"source"}) {
        $line = $self->_format_header_lines("source",$header);
        push(@header_lines,$line->[0]);
    }
    
    # Write reference if present
    if ($header->{"reference"}) {
        $line = $self->_format_header_lines("reference",$header);
        push(@header_lines,$line->[0]);
    }
    
    # Write assembly if present
    if ($header->{"assembly"}) {
        $line = $self->_format_header_lines("assembly",$header);
        push(@header_lines,$line->[0]);
    }
    
    # Write phasing if present
    if ($header->{"phasing"}) {
        $line = $self->_format_header_lines("phasing",$header);
        push(@header_lines,$line->[0]);
    }
    
    # Write commandline if present
    if ($header->{"commandline"}) {
        $line = $self->_format_header_lines("commandline",$header);
        push(@header_lines,$line->[0]);
    }
    
    # Write PEDIGREE if present
    if ($header->{"PEDIGREE"}) {
        $line = $self->_format_header_lines("PEDIGREE",$header);
        push(@header_lines,$line->[0]);
    }
    
    # Write pedigreeDB if present
    if ($header->{"pedigreeDB"}) {
        $line = $self->_format_header_lines("pedigreeDB",$header);
        push(@header_lines,$line->[0]);
    }
    
    # Write contig if present
    if ($header->{"contig"}) {
        $line = $self->_format_header_lines("contig",$header);
        push(@header_lines,@{$line});
    }
    
    # Write INFO if present
    if ($header->{"INFO"}) {
        $line = $self->_format_header_lines("INFO",$header);
        push(@header_lines,@{$line});
    }
    
    # Write FORMAT if present
    if ($header->{"FORMAT"}) {
        $line = $self->_format_header_lines("FORMAT",$header);
        push(@header_lines,@{$line});
    }
    
    # Write FILTER if present
    if ($header->{"FILTER"}) {
        $line = $self->_format_header_lines("FILTER",$header);
        push(@header_lines,@{$line});
    }
    
    # Write ALT if present
    if ($header->{"ALT"}) {
        $line = $self->_format_header_lines("ALT",$header);
        push(@header_lines,@{$line});
    }
    
    # Write SAMPLE if present
    if ($header->{"SAMPLE"}) {
        $line = $self->_format_header_lines("SAMPLE",$header);
        push(@header_lines,@{$line});
    }
    
    # Write FILTER if present
    if ($header->{"META"}) {
        $line = $self->_format_header_lines("META",$header);
        push(@header_lines,@{$line});
    }
    
    # Write other if present
    if ($header->{"other"}) {
        foreach my $o (keys(%{$header->{"other"}})) {
            $line = $self->_format_header_lines($o,$header->{"other"});
        }
        push(@header_lines,@{$line});
    }
    
    return(\@header_lines);
}

sub _format_header_lines {
    my ($self,$key,$header) = @_;
    
    my @lines;
    my $data = $header->{$key};
    my $ref = ref($data);
    if ($ref) {
        if ($ref eq "ARRAY") { # Then it's array of hashes
            foreach my $kvh (@{$data}) {
                my $line = $self->_format_header_line($key,$kvh);
                push(@lines,$line);
            }
        }
        elsif ($ref eq "HASH") {
            # Only "other" is HASH... Maybe safe to remove at some point
        }
    }
    else {
        my $line = $self->_format_header_line($key,$data);
        push(@lines,$line);
    }
    
    return(\@lines);
}

sub _format_header_line {
    my ($self,$key,$ref) = @_;
    if (defined($ref)) {
        if (ref($ref) eq "HASH") {
            my @conts;
            while (my ($k,$v) = each (%{$ref})) {
                push(@conts,"$k=$v");
            }
            return("##$key=<".join(",",@conts).">");
        }
        else {
            return("##$key=$ref");
        }
    }
}

sub _get_column_names {
    my $self = $_[0];
    return($self->get_header->{"columns"});
}

sub _prepare_new_header_field_info {
    my ($self,$ra) = @_;
    
    my %attr;
    tie %attr,"Tie::IxHash::Easy";
    
    # Get fields other than the 4 mandatory
    my @mandatory = ("ID","Number","Type","Description");
    my %in_mandatory = map {$_ => 1} @mandatory;
    my @not_mandatory = grep {not $in_mandatory{$_}} keys(%{$ra});
    
    $attr{"ID"} = $ra->{"ID"} if (defined($ra->{"ID"}));
    $attr{"Number"} = $ra->{"Number"} if (defined($ra->{"Number"}));
    $attr{"Type"} = $ra->{"Type"} if (defined($ra->{"Type"}));
    $attr{"Description"} = $ra->{"Description"} 
        if (defined($ra->{"Description"}));
    
    # If there are no mandatory fields (VCFv>=4.3), sort and include. If not
    # allowed by the version, validation routines will take care.   
    if (scalar @not_mandatory > 0) {
        my @snot_mandatory = sort @not_mandatory;
        foreach my $s (@snot_mandatory) {
            $attr{$s} = $ra->{$s};
        }
    }
    
    return(\%attr);
}

sub _prepare_new_header_field_format {
    my ($self,$ra) = @_;
    
    my %attr;
    tie %attr,"Tie::IxHash::Easy";
    
    $attr{"ID"} = $ra->{"ID"} if (defined($ra->{"ID"}));
    $attr{"Number"} = $ra->{"Number"} if (defined($ra->{"Number"}));
    $attr{"Type"} = $ra->{"Type"} if (defined($ra->{"Type"}));
    $attr{"Description"} = $ra->{"Description"} 
        if (defined($ra->{"Description"}));
    
    return(\%attr);
}

sub _prepare_new_header_field_filter {
    my ($self,$ra) = @_;
    
    my %attr;
    tie %attr,"Tie::IxHash::Easy";
    
    $attr{"ID"} = $ra->{"ID"} if (defined($ra->{"ID"}));
    $attr{"Description"} = $ra->{"Description"} 
        if (defined($ra->{"Description"}));
    
    return(\%attr);
}

sub _prepare_new_header_field_alt {
    my ($self,$ra) = @_;
    
    my %attr;
    tie %attr,"Tie::IxHash::Easy";
    
    $attr{"ID"} = $ra->{"ID"} if (defined($ra->{"ID"}));
    $attr{"Description"} = $ra->{"Description"} 
        if (defined($ra->{"Description"}));
    
    return(\%attr);
}

# TODO: Check what SAMPLE fields are mandatory
sub _prepare_new_header_field_sample {
    my ($self,$ra) = @_;
    
    my %attr;
    tie %attr,"Tie::IxHash::Easy";
        
    # Get fields other than the 2 mandatory
    my @mandatory = ("ID","Description");
    my %in_mandatory = map {$_ => 1} @mandatory;
    my @not_mandatory = grep {not $in_mandatory{$_}} keys(%{$ra});
    
    $attr{"ID"} = $ra->{"ID"} if (defined($ra->{"ID"}));
    $attr{"Description"} = $ra->{"Description"} 
        if (defined($ra->{"Description"}));
    
    # If there are no mandatory fields (VCFv>=4.3), sort and include. If not
    # allowed by the version, validation routines will take care.   
    if (scalar @not_mandatory > 0) {
        my @snot_mandatory = sort @not_mandatory;
        foreach my $s (@snot_mandatory) {
            $attr{$s} = $ra->{$s};
        }
    }
    
    return(\%attr);
}

# TODO: Check what META fields are mandatory - not yet clear
sub _prepare_new_header_field_meta {
    my ($self,$ra) = @_;
    
    my %attr;
    tie %attr,"Tie::IxHash::Easy";
        
    # Get fields other than the 2 mandatory
    my @mandatory = ("ID","Type","Number","Values");
    my %in_mandatory = map {$_ => 1} @mandatory;
    my @not_mandatory = grep {not $in_mandatory{$_}} keys(%{$ra});
    
    $attr{"ID"} = $ra->{"ID"} if (defined($ra->{"ID"}));
    $attr{"Type"} = $ra->{"Type"} if (defined($ra->{"Type"}));
    $attr{"Number"} = $ra->{"Number"} if (defined($ra->{"Number"}));
    $attr{"Values"} = $ra->{"Values"} if (defined($ra->{"Values"}));
    
    if (scalar @not_mandatory > 0) {
        my @snot_mandatory = sort @not_mandatory;
        foreach my $s (@snot_mandatory) {
            $attr{$s} = $ra->{$s};
        }
    }
    
    return(\%attr);
}

sub _parse_sample_names_internal {
    my ($self,$line) = @_;
    
    my @fields = split("\t",$line);
    my $n = scalar @fields;
    my @nf = @fields[9..($n-1)];
    return(\@nf) if ($fields[8] && $fields[8] eq "FORMAT");
    return([]);
}

sub _format_column_names {
    my $self = shift @_;
    return(join("\t",@_));
}

sub _has_header {
    my $self = $_[0];
    return(1) if (defined($self->{"_content"}));
    return(0);
}

sub _init_header_scaffold {
    my ($self,$ver) = @_;
    my $scaf = $self->_define_header_scaffold;
    return($scaf->{$ver});
}

sub _init_header_metanames {
    my ($self,$ver) = @_;
    my $mn = $self->_define_header_metanames;
    return($mn->{$ver});
}

sub _init_header_expr {
    my ($self,$ver) = @_;
    my $specs = $self->_define_header_expr;
    return($specs->{$ver});
}

sub _define_header_scaffold {
    my %scaffold;
    tie %scaffold, "Tie::IxHash::Easy";
    
    $scaffold{"VCFv4.1"}{"fileformat"} = undef;
    $scaffold{"VCFv4.1"}{"fileDate"} = undef;
    $scaffold{"VCFv4.1"}{"source"} = undef;
    $scaffold{"VCFv4.1"}{"reference"} = undef;
    $scaffold{"VCFv4.1"}{"assembly"} = undef;
    $scaffold{"VCFv4.1"}{"phasing"} = undef;
    $scaffold{"VCFv4.1"}{"commandline"} = undef;
    $scaffold{"VCFv4.1"}{"pedigreeDB"} = undef;
    $scaffold{"VCFv4.1"}{"PEDIGREE"} = undef;
    $scaffold{"VCFv4.1"}{"contig"} = [];
    $scaffold{"VCFv4.1"}{"INFO"} = [];
    $scaffold{"VCFv4.1"}{"FORMAT"} = [];
    $scaffold{"VCFv4.1"}{"FILTER"} = [];
    $scaffold{"VCFv4.1"}{"ALT"} = [];
    $scaffold{"VCFv4.1"}{"SAMPLE"} = [];
    $scaffold{"VCFv4.1"}{"other"} = {};
    $scaffold{"VCFv4.1"}{"columns"} = [];
    
    $scaffold{"VCFv4.2"}{"fileformat"} = undef;
    $scaffold{"VCFv4.2"}{"fileDate"} = undef;
    $scaffold{"VCFv4.2"}{"source"} = undef;
    $scaffold{"VCFv4.2"}{"reference"} = undef;
    $scaffold{"VCFv4.2"}{"assembly"} = undef;
    $scaffold{"VCFv4.2"}{"phasing"} = undef;
    $scaffold{"VCFv4.2"}{"commandline"} = undef;
    $scaffold{"VCFv4.2"}{"pedigreeDB"} = undef;
    $scaffold{"VCFv4.2"}{"PEDIGREE"} = undef;
    $scaffold{"VCFv4.2"}{"contig"} = [];
    $scaffold{"VCFv4.2"}{"INFO"} = [];
    $scaffold{"VCFv4.2"}{"FORMAT"} = [];
    $scaffold{"VCFv4.2"}{"FILTER"} = [];
    $scaffold{"VCFv4.2"}{"ALT"} = [];
    $scaffold{"VCFv4.2"}{"SAMPLE"} = [];
    $scaffold{"VCFv4.2"}{"other"} = {};
    $scaffold{"VCFv4.2"}{"columns"} = [];
    
    $scaffold{"VCFv4.3"}{"fileformat"} = undef;
    $scaffold{"VCFv4.3"}{"fileDate"} = undef;
    $scaffold{"VCFv4.3"}{"source"} = undef;
    $scaffold{"VCFv4.3"}{"reference"} = undef;
    $scaffold{"VCFv4.3"}{"assembly"} = undef;
    $scaffold{"VCFv4.3"}{"phasing"} = undef;
    $scaffold{"VCFv4.3"}{"commandline"} = undef;
    $scaffold{"VCFv4.3"}{"pedigreeDB"} = undef;
    $scaffold{"VCFv4.3"}{"PEDIGREE"} = undef;
    $scaffold{"VCFv4.3"}{"contig"} = [];
    $scaffold{"VCFv4.3"}{"INFO"} = [];
    $scaffold{"VCFv4.3"}{"FORMAT"} = [];
    $scaffold{"VCFv4.3"}{"FILTER"} = [];
    $scaffold{"VCFv4.3"}{"ALT"} = [];
    $scaffold{"VCFv4.3"}{"SAMPLE"} = [];
    $scaffold{"VCFv4.3"}{"META"} = [];
    $scaffold{"VCFv4.3"}{"other"} = {};
    $scaffold{"VCFv4.3"}{"columns"} = [];
    
    return(\%scaffold);
}

sub _define_header_metanames {
    my %metanames = (
        "VCFv4.1" => [
            "fileformat","fileDate","source","assembly","phasing","commandline",
            "pedigreeDB","contig","reference","PEDIGREE","INFO","FORMAT",
            "FILTER","ALT","SAMPLE"
        ],
        "VCFv4.2" => [
            "fileformat","fileDate","source","assembly","phasing","commandline",
            "pedigreeDB","contig","reference","PEDIGREE","INFO","FORMAT",
            "FILTER","ALT","SAMPLE"
        ],
        "VCFv4.3" => [
            "fileformat","fileDate","source","assembly","phasing","commandline",
            "pedigreeDB","contig","reference","PEDIGREE","INFO","FORMAT",
            "FILTER","ALT","SAMPLE","META"
        ]
    );
    return(\%metanames);
}

sub _define_header_expr {
    my %specs = (
        "VCFv4.1" => {
            "INFO" => {
                "Number" => qr/^\d?$|^[AG]$|^\.$/
            },
            "FORMAT" => {
                "Number" => qr/^\d?$|^[AG]$|^\.$/
            }
        },
        "VCFv4.2" => {
            "INFO" => {
                "Number" => qr/^\d?$|^[AGR]$|^\.$/
            },
            "FORMAT" => {
                "Number" => qr/^\d?$|^[AGR]$|^\.$/
            }
        },
        "VCFv4.3" => {
            "INFO" => {
                "Number" => qr/^\d?$|^[AGR]$|^\.$/,
                "ID" => qr/^([A-Za-z ][0-9A-Za-z .]*|1000G)$/
            },
            "FORMAT" => {
                "Number" => qr/^\d?$|^[AGR]$|^\.$/
            },
            "META" => {
                "Number" => qr/^\d?$|^[AGR]$|^\.$/
            }
        }
    );
    return(\%specs);
}

sub _hack_comma_in_quotes {
    my ($self,$str) = @_;
    
    # Whether Description tag is in the middle or in the end, the description
    # text will always be $vals[1]. If there are no quotes, value remains
    # unchanged, we simply check if $vals[1] exists.
    my @vals = split("\"",$str);
    if ($vals[1]) {
        $vals[1] =~ s/,/___COMMA___/g;
    }
    return(join("",@vals));
}

sub _unhack_comma_in_quotes {
    my ($self,$str) = @_;
    $str =~ s/___COMMA___/,/g;
    return($str);
}

sub _init {
    my ($self,$params) = @_;
    
    $self->{'validation'} = $params->{'validation'};
    $self->{'file'} = $params->{'file'};
    
    if ($params->{'data'}) {
        $self->{'_content'} = {%{$params->{'data'}}}; # Copy
        undef($params->{'data'});
        #$self->validate;
    }
    else {
        $self->{'_content'} = undef;
    }
    
    $self->{'_output'} = undef;
    
    return($self);
}

sub _check_params {
    my ($self,$params) = @_;
    
    # $params can be a hash reference with parameters or an existing filename
    # If the second, create the $params with this file
    if (-f $params) {
        $params = {'file' => $params};
    }
    
    # Check fatal
    my $stop;
    $stop .= "--- Please specify VCF file to read/write or data ---\n" 
        if (!$params->{"file"} && !$params->{"data"});
        
    if ($stop) {
        print STDERR "\n$stop\n";
        croak "Type perldoc $MODNAME for help in usage.\n\n";
        exit;
    }
    
    # Accepted parameters
    my @accept = ("file","data","validation");
    
    # Check and warn for unrecognized parameters
    foreach my $p (keys(%{$params})) {
        carp "\nUnrecognized parameter : $p   --- Ignoring..." 
            if (!$helper->_smatch($p,@accept));
    }
    
    # Warn if both data and file provided
    if ($params->{"file"} && $params->{"data"}) {
        carp "\nBoth file and data provided, file will be used.\n";
        undef($params->{"data"});
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

    perldoc VCF::Header


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

1; # End of VCF::Header
