package VCF::Record;

use 5.10.0;
use strict;
use warnings;
use utf8;

use Carp;
use File::Basename;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Tie::IxHash::Easy;
use List::Util qw(first);

use VCF::_Helper;

use vars qw($helper);

BEGIN {
    $helper = VCF::_Helper->new;
    binmode STDOUT, ":utf8";
    $|=1;
}

use Data::Dumper qw(Dumper);

=head1 NAME

VCF::Record - The great new VCF::Record!

=head1 VERSION

Version 0.01

=cut

our $MODNAME = "VCF::Record";
our $VERSION = '0.01';

# ID checker
our %GLOBAL_IDS;
# Allele uniqueness
our %META_ALLELES;

=head1 SYNOPSIS

VCF::Record provides methods for parsing, validating, manipulating and writing
VCF file records and meta-information lines based on the VCF file 
specifications. A VCF::Record object can be initialized by providing a VCF file
as a single input, or by providing a hash with parameters with the following 
posssible keys: C<file>, C<data>, C<validation>.

C<file> is a valid, existing VCF file, C<data> is a hash with header data. To
see an example, just parse a VCF file and inspect the has returned by the
function C<get_header>. If both C<file> and C<data> are provided, C<file> will
be used. C<validation> sets the header validation strictness.

Usage:
    
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

=head1 SUBROUTINES/METHODS

=head2 new

The VCF::Record object constructor. It accepts an existing VCF file.

    my $record = VCF::Record->new({
        'file' => 'my_variant.vcf'
    });
    my $record = VCF::Record->new('my_variant.vcf');

=cut

sub new {
    my ($class,$params) = @_;
    my $self = {};
    
    # Bless so that we start using methods (including parameter validation)
    bless($self,$class);
    
    # Validate the input parameters
    $params = $self->_check_params($params);
    
    # If all OK, continue with object initialization
    $self->_init($params);
    return($self);
}

=head2 parse_all

VCF::Record object simple OO parsing of all records in a VCF file. This method 
does not perform any essential validations on record corectness.

    $record->parse_all;

=cut

sub parse_all {
    my $self = $_[0];
    my $vcf = $self->get('file');
    
    my @records;
    open(VCF_FILE,$vcf);
    while (my $line = <VCF_FILE>) {
        next if ($line =~ m/^#/);
        $line =~ s/\r|\n$//g;
        print STDERR "  parsed $. entries\r" if ($.%10000 == 0);
        my $fields = $self->parse_one($line,$.);
        #my $u = $self->_check_unique_fields($fields,$.);
        #if ($u) {
        #    close(VCF_FILE);
        #    die "\nCorrupted VCF ID entry, line $.: ID $u has been found ".
        #        "more than once!";
        #}
        push(@records,$fields);
    }
    close(VCF_FILE);
    
    $self->set_records(\@records);
}

=head2 validate_all

VCF::Record object validation of a VCF file records. This method performs 
record field validations according to each supported file format
specifications and according to validation level. It can be used when
parsing and validating at separate stages, e.g. when a record must be added to
the file. A VCF::Header object is required for validation.

    $record->validate_all;

=cut

sub validate_all {
    my ($self,$header) = @_;
    
    croak "\nNo header provided! Cannot validate VCF records. Use parse for ".
        "just parsing.\n" if (!$header);
    croak "\nNo records found in VCF::Record object! have you parsed a VCF ".
        "file first?\n" if (!$self->has_records);
    
    my $records = $self->get_records;
    
    my ($header_ids,$sample_names) = $self->_extract_validation_data($header);
    
    my $n = @{$records};
    for (my $i=0; $i<$n; $i++) {
        print STDERR "  validated $i entries\r" if (($i+1)%10000 == 0);
        $records->[$i] = 
            $self->_validate_one($records->[$i],$header_ids,$sample_names);
    }
    
    $self->set_records($records);

    return($self);
}

=head2 parse_and_validate_all

VCF::Record object parsing and validation of a VCF file records. This method
method performs record validations according to each supported file format
specifications and according to validation level. A VCF::Header object is 
required for validation.

    $header->parse_and_validate($header);

=cut

sub parse_and_validate_all {
    my ($self,$header) = @_;
    my $vcf = $self->get('file');
    
    $self->parse_all;
    $self->validate_all($header);
    
    ## This implementation was 5x slower because _extract_validation_data was
    ## called #records times instead of once! parse_and_validate_one remains as
    ## API but not for bulk usage.
    #my @records;
    #open(VCF_FILE,$vcf);
    #while (my $line = <VCF_FILE>) {
    #    next if ($line =~ m/^#/);
    #    $line =~ s/\r|\n$//g;
    #    print STDERR "  processed $. entries\r" if ($.%10000 == 0);
    #    my $fields = $self->parse_and_validate_one($line,$header,$.);
    #    push(@records,$fields);
    #}
    #close(VCF_FILE);
    #$self->set_records(\@records);
    
    return(0);
}

=head2 parse_one

Create a hash with VCF data for a single record by parsing a VCF record line. 
This method does not perform any essential validations on record corectness.

    $record->parse_one;

=cut

sub parse_one {
    my ($self,$line,$ln) = @_;
    
    $ln = " ".$ln if ($ln);
    croak "\nCorrupted VCF line$ln! Empty lines are not permitted."
        if ($line =~/^\s*$/);
    
    my %fields = (
        "_id" => undef,
        "CHROM" => undef,
        "POS" => undef,
        "ID" => undef,
        "REF" => undef,
        "ALT" => [],
        "QUAL" => undef,
        "FILTER" => [],
        "INFO" => {},
        "FORMAT" => [],
        "SAMPLE" => []
    );
    
    my @cols = split(/\t/,$line);
    
    my $n = @cols;
    croak "\nCorrupted VCF line$ln! $n fields found, at least 8 expected."
        if ($n < 8);
    
    # We have simple parsing without validation. So, the validation parameter
    # must be temporarily switched to "basic"
    my $current_level = $self->get("validation");
    $self->set('validation',"basic");
    
    # A mock hash as argument (e.g. $contig to chromosome validation)
    my %mock;
    
    # A mock VCF version - because of no validation, header is not really
    # required
    my $vcfv = "VCFv4.1";
    
    $fields{"CHROM"} = $self->_validate_record_field_chrom($cols[0],\%mock,
        $ln);
    $fields{"POS"} = $self->_validate_record_field_pos($cols[1],
        $fields{"CHROM"},\%mock,$ln);
    $fields{"ID"} = $self->_validate_record_field_id($cols[2],$ln);
    $fields{"REF"} = $self->_validate_record_field_ref($cols[3],$ln);
    $fields{"ALT"} = $self->_validate_record_field_alt($cols[4],$fields{"REF"},
        [],$vcfv,$ln);
    $fields{"QUAL"} = $self->_validate_record_field_qual($cols[5],$ln);
    $fields{"FILTER"} = $self->_validate_record_field_filter($cols[6],[],$ln);
    $fields{"INFO"} = $self->_validate_record_field_info($cols[7],
        $fields{"ALT"},[],$ln);
    
    if ($n > 8) { # Continue with validation
        croak "\nCorrupted VCF line $ln! $n fields found, at least 10 ".
            "expected when FORMAT column is present." if ($n < 10);
        $fields{"FORMAT"} = $self->_validate_record_field_format($cols[8],
            [],$ln);
        $fields{"SAMPLE"} = $self->_validate_record_field_sample(
            join("\t",@cols[9..($n-1)]),[],$fields{"FORMAT"},$fields{"INFO"},
                $fields{"ALT"},[],$ln);
    }
    
    # Finally construct the internal id
    my $sep = "^^^"; # A strange separator to avoid conflict with data
    $fields{"_id"} = $fields{"CHROM"}.$sep.$fields{"POS"}.$sep.$fields{"REF"}.
        $sep.join(",",@{$fields{"ALT"}});
    
    return(\%fields);
}

=head2 validate_one

Validate a hash with VCF data for a single record parsed by C<parse_record>. In
order to perform validation, a VCF::Header object must be also provided. It can
be the VCF::Header object or the result of C<VCF::Header::get_header>. It's
better to provide the VCF::Header object, as it's usually validated whille in 
the other case, validation will be reperformed according to validation level. 

    $record->validate_one($fields,$header);

=cut

sub validate_one {
    my ($self,$fields,$header,$header_ids,$sample_names) = @_;
    # By also passing $header_ids and $sample_names as optional arguments, we 
    # can use validate_one instead of the internal _validate_one with the same
    # speed.
    
    croak "\nNo header provided! Cannot validate VCF record.\n"if (!$header);
    
    ($header_ids,$sample_names) = $self->_extract_validation_data($header)
        if (!$header_ids && !$sample_names);
    
    return($self->_validate_one($fields,$header_ids,$sample_names));
}

=head2 parse_and_validate_one

Create a hash with VCF data for a single record by parsing a VCF record line,
while validating on record corectness according to the validation level
declared on VCF::Record object initialization. In order to perform validation,
a VCF::Header object must be also provided. It can be the VCF::Header object 
or the result of C<VCF::Header::get_header>. It's better to provide the 
VCF::Header object, as it's usually validated whille in  the other case, 
validation will be reperformed according to validation level.

    $record->parse_and_validate_one($line,$header);

=cut

sub parse_and_validate_one {
    my ($self,$line,$header,$ln) = @_;
    
    croak "\nCorrupted VCF line $ln! Empty lines are not permitted."
        if ($line =~/^\s*$/);
        
    croak "\nNo header provided! Cannot validate VCF record.\n"if (!$header);
        
    my ($header_ids,$sample_names) = $self->_extract_validation_data($header);
    
    my @cols = split(/\t/,$line);
    
    my $n = @cols;
    $ln = " ".$ln if ($ln);
    croak "\nCorrupted VCF line$ln! $n fields found, at least 8 expected."
        if ($n < 8);
    
    # Temporary fields, normal form after validation
    my %fields = (
        "_id" => undef,
        "CHROM" => $cols[0],
        "POS" => $cols[1],
        "ID" => $cols[2],
        "REF" => $cols[3],
        "ALT" => $cols[4],
        "QUAL" => $cols[5],
        "FILTER" => $cols[6],
        "INFO" => $cols[7],
        "FORMAT" => [],
        "SAMPLE" => []
    );
    
    if ($n > 8) { # Continue with validation
        croak "\nCorrupted VCF line$ln! $n fields found, at least 10 ".
            "expected when FORMAT column is present." if ($n < 10);
        $fields{"FORMAT"} = $cols[8];
        $fields{"SAMPLE"} = join("\t",@cols[9..($n-1)]);
    }
    
    # The internal id construction takes place inside _validate_one
    return($self->_validate_one(\%fields,$header_ids,$sample_names));
}

=head2 add_key_value_attr

Adds a new key-value attribute to a record in fields that allow this, that is
the INFO and SAMPLE (FORMAT). The key must have been added to the
VCF header first (using the VCF::Header::add method), otherwise the validation
will fail.

We remind that the SAMPLE attribute is a helper VCF::Record attribute for
storing each VCF samples' genotypes and other data.
    
    my $rec = $record->get_records->[0];
    my %kv = (
        MYANN => ["A|geneA|dis1^dis2^dis3|phen1^phen2^phen3"],
        MYSCORE => 0.458
    );
    $record->add_key_value_attr(\%kv,"INFO",$header,$rec);

=cut

sub add_key_value_attr {
    my ($self,$keyval,$type,$header,$record) = @_;
    
    croak "\nYou did not provide a record hash to add the key-value pair!\n"
        if (!$record);
    croak "\nYou did not provide a header object for validation!\n"
        if (!$header);
    croak "\nPlease provide a valid type to add the key-value pair!\n"
        if (!$type);
    
    # Get validation elements
    my ($header_ids,$sample_names) = $self->_extract_validation_data($header);
    
    # 1. Check that the type is valid AND exists in the header. 
    my @allow_add_fields = ("INFO","SAMPLE");
    croak "\nInvalid or immutable type $type! You are allowed to add ".
        "key-value pair attributes only to INFO, SAMPLE fields."
            if (!grep {$_ eq $type} @allow_add_fields);
    my @header_keys = keys(%{$header_ids});
    croak "\nThe type $type that you are trying to add must exist in the VCF ".
        "header!" if (!grep {$_ eq $type} @header_keys);
    
    # 2. Check that the key is valid (from the header). If we want to add an
    #    attribute to a record, this must have been added to the header first
    #    through add_vcf_header_record
    my %kv = %{$keyval};
    my @keys = keys(%kv);
    my @curr_keys = keys(%{$header_ids->{$type}});
    foreach my $k (@keys) {
        croak "You are trying to add $type attribute $k which has not been ".
            "yet added to the VCF header!" if (!(grep{$_ eq $k} @curr_keys));
    }
    
    # 3. Verify that the keys do not already exist in the record. If existing,
    #    a set-like function must be used.
    my @ex = keys(%{$record->{$type}});
    foreach  my $k (@keys) {
        croak "You are trying to add $type attribute $k which already exists! ".
            "Please use the change function!" if (grep{$_ eq $k} @ex);
    }
    
    # 4. If above checks OK, add data to the VCF record
    foreach  my $k (@keys) {
        my $v = ($helper->_is_array_ref($kv{$k})) ? ($kv{$k}) : ([$kv{$k}]);
        $record->{$type}->{$k} = $v;
    }
    
    # 5. Verify that the key-value pairs respect the Type definition from the
    #    header by re-inspecting the record
    if ($type eq "INFO") {
        $record->{"INFO"} = $self->_validate_record_field_info($record,$record,
            $header_ids->{"INFO"});
    }
    elsif ($type eq "SAMPLE") {
        $record->{"SAMPLE"} = $self->_validate_record_field_sample($record,
            $record,$record->{"FORMAT"},$record->{"INFO"},$record->{"ALT"},
            $header_ids->{"FORMAT"});
    }
    
    return($record);
}

=head2 add_key_attr

Adds a new key attribute to a record in fields that allow this, that is the ALT
FILTER and FORMAT attributes. The key must have been added to the VCF header
first (using the VCF::Header::add method), otherwise the validation
will fail.
    
    my $rec = $record->get_records->[0];
    my @f1 = ("q10","p30");
    $record->add_key_attr(\@f1,"FILTER",$header,$rec);

=cut

sub add_key_attr {
    my ($self,$keyref,$type,$header,$record) = @_;
    
    croak "\nYou did not provide a record hash to add the key!\n"
        if (!$record);
    croak "\nYou did not provide a header object for validation!\n"
        if (!$header);
    croak "\nPlease provide a valid type to add the key!\n"
        if (!$type);
    
    # Get validation elements
    my ($header_ids,$sample_names) = $self->_extract_validation_data($header);
    
    # 1. Check that the type is valid AND exists in the header. 
    my @allow_add_fields = ("ALT","FILTER","FORMAT");
    croak "\nInvalid or immutable type $type! You are allowed to add a value ".
        "attribute only to ALT, FILTER, FORMAT fields." 
            if (!grep {$_ eq $type} @allow_add_fields);
    my @header_keys = keys(%{$header_ids});
    croak "\nThe type $type that you are trying to add must exist in the VCF ".
        "header!" if (!grep {$_ eq $type} @header_keys);
    
    # 2. Check that the key is valid (from the header). If we want to add an
    #    attribute to a record, this must have been added to the header first
    #    through add_vcf_header_record
    my @keys = @{$keyref};
    my @curr_keys = ($helper->_is_array_ref($header_ids->{$type})) ? 
        @{$header_ids->{$type}} : keys(%{$header_ids->{$type}});
    foreach my $k (@keys) {
        croak "You are trying to add $type attribute $k which has not been ".
            "yet added to the VCF header!" if (!(grep{$_ eq $k} @curr_keys));
    }
    
    # 3. Verify that the keys do not already exist in the record. If existing,
    #    a set-like function must be used.
    my @ex = @{$record->{$type}};
    foreach  my $k (@keys) {
        croak "You are trying to add $type attribute $k which already exists!"
            if (grep{$_ eq $k} @ex);
    }
    
    # 4. If above checks OK, add data to the VCF record
    # TODO: Consider removing the FILTER addition here and move to a set 
    #       function
    if ($type eq "FILTER") {
        if ($record->{$type}->[0] eq "." || $record->{$type}->[0] eq "PASS") {
            # We must remove the missing value or just the PASS, so we simply
            # replace.
            $record->{$type} = undef;
            push(@{$record->{$type}},@keys);
        }
        else {
            push(@{$record->{$type}},@keys);
        }
    }
    else {
        push(@{$record->{$type}},@keys);
    }
    
    # 5. Verify that the key-value pairs respect the Type definition from the
    #    header by re-inspecting the record
    if ($type eq "ALT") {
        my @alt_ids = @{$header_ids->{"ALT"}};
        if (@alt_ids > 0) {
            for (my $i=0; $i<@alt_ids; $i++) {
                $alt_ids[$i] = "<".$alt_ids[$i].">";
            }
        }
        $record = $self->_validate_record_field_alt($record,$record,
            \@alt_ids,$header_ids->{"fileformat"});
    }
    elsif ($type eq "FILTER") { 
        $record->{"FILTER"} = $self->_validate_record_field_filter($record,
            $header_ids->{"FILTER"});
    }
    elsif ($type eq "FORMAT") {
        $record->{"FORMAT"} = $self->_validate_record_field_format($record,
            $header_ids->{"FORMAT"});
    }
    
    return($record);
}

=head2 add

Validates a new record and adds it to the VCF records object. The argument is
a hash with the record information to be added. The record will be validated
according to the validation strictness initialized with the VCF::Record object
and then added to the collection
    
    $record->add($r);
    
TODO: Proper example

=cut

sub add {
    my ($self,$new_rec,$header) = @_;
    
    croak "\nYou did not provide a header object for validation!\n"
        if (!$header);
    
    my $records = $self->get_records;
    
    # Validate the input
    my $val_rec = $self->validate_one($new_rec,$header);
    
    # Construct the record to be added
    my %rec = (
        "CHROM" => $val_rec->{"CHROM"},
        "POS" => $val_rec->{"POS"},
        "ID" => $val_rec->{"ID"},
        "REF" => $val_rec->{"REF"},
        "ALT" => $val_rec->{"ALT"},
        "QUAL" => $val_rec->{"QUAL"},
        "FILTER" => $val_rec->{"CHROM"},
        "INFO" => {},
        "FORMAT" => [],
        "SAMPLE" => []
    );
    
    if ($records && $records->[0]) { # Get from a record
        my $arec = $records->[0];
        # For INFO, FORMAT, SAMPLE, use a tied hash to store order of 
        # existing record. We assume that validation has taken care of
        # missing keys etc.
        
        my ($tmp1,$tmp2);
        
        # INFO fields must respect the other records
        my %info = %{$arec->{"INFO"}};
        # Check if order of keys is the same already, by serializing the arrays
        $tmp1 = join("__",keys(%info));
        $tmp2 = join("__",keys(%{$val_rec->{"INFO"}}));
        if ($tmp1 eq $tmp2) { # Then simply add it
            $rec{"INFO"} = $val_rec->{"INFO"};
        }
        else {
            my %fi_rec;
            tie %fi_rec,"Tie::IxHash::Easy";
            foreach my $key (keys(%info)) {
                $fi_rec{$key} = $val_rec->{"INFO"}->{$key};
            }
            $rec{"INFO"} = \%fi_rec;
        }
        
        # FORMAT keys must have the same order
        my @format = @{$arec->{"FORMAT"}};
        # Check if order of keys is the same already, by serializing the arrays
        $tmp1 = join("__",@format);
        $tmp2 = join("__",@{$val_rec->{"FORMAT"}});
        if ($tmp1 eq $tmp2) { # Then simply add it
            $rec{"FORMAT"} = $val_rec->{"FORMAT"};
        }
        else {
            # Just copy the current, as this field is fixed
            $rec{"FORMAT"} = $arec->{"FORMAT"};
        }
        
        # SAMPLE keys etc. must have the same order...
        my @sample = @{$arec->{"SAMPLE"}};
        my @s_add;
        foreach my $s (@{$val_rec->{"SAMPLE"}}) {
            $tmp1 = join("__",keys(%{$sample[0]}));
            $tmp2 = join("__",keys(%{$s}));
            if ($tmp1 eq $tmp2) { # Then simply add it
                push(@s_add,$s);
            }
            else {
                my %fi_rec;
                tie %fi_rec,"Tie::IxHash::Easy";
                foreach my $key (keys(%{$sample[0]})) {
                    $fi_rec{$key} = $s->{$key};
                }
                push(@s_add,$s);
            }
        }
        $rec{"SAMPLE"} = \@s_add;
    }
    else { # 1st record? From header?
        # TODO: Code this part, $records is empty, so order from header
    }
    
    # Finally, add the new record to the existing records
    push(@{$records},\%rec);
    $self->set_records($records);
}

=head2 change

Change the attributes of a given record. This function will B<not> add new keys
and values to an existing record. Use the C<add_key_attr> and 
C<add_key_value_attr> methods for this.

The C<change> method must generally be used with much caution.
    
B<This method is pending implementation>

=cut

# $with must be either a hash with $type, $key, $value elements, e.g.
# $with = {
#     'type' => 'INFO',
#     'key' => 'DP',
#     'value' => 25
# };
# $with = {
#     'type' => 'FORMAT',
#     'value' => "GT:DP:AF"
# };
# $with = {
#     'type' => 'FILTER',
#     'value' => "PASS"
# };

sub change {
    my ($self,$i,$what,$with) = @_;
    # Stub
}

=head2 lookup

A very simple lookup method to retrieve a record based on a simple internal id
consisting of I<CHROM>_I<POS>_I<REF>_I<ALT>. If I<ALT> is multi-allelic, the
different alleles are joind by comma (,). Thus, the I<ALT> to be looked up
must be an array reference (see example below).

Currently B<ALL> the fields of $search must be specified.

my $search = {
    'chrom' => "1",
    'pos' => 123456,
    'ref' => "A",
    'alt' => ["G"]
}
$record->lookup($search);

TODO: A more sophisticated find method looking with additional criteria, like
      a filter value, a range of coordinates etc. The input VCF will have to be
      sorted.

=cut

sub lookup {
    my ($self,$search) = @_;
    
    # This block will be removed in a future more flexible lookup
    if (!exists($search->{"chrom"}) || !$search->{"chrom"}
        || !exists($search->{"pos"}) || !$search->{"pos"}
        || !exists($search->{"ref"}) || !$search->{"ref"}
        || !exists($search->{"alt"}) || !$search->{"alt"}) {
        carp "\nAll search terms must be specified! Returning empty...\n";
        return([]);
    }
    
    my $str = "";
    my $sep = "^^^";
    $str = $search->{"chrom"}.$sep.$search->{"pos"}.$sep.$search->{"ref"}.$sep.
        join(",",@{$search->{"alt"}});
    
    # This block will be activated in a future more flexible lookup
    #$str .= $search->{"chrom"}.$sep
    #   if (exists($search->{"chrom"}) && $search->{"chrom"});
    #$str .= $search->{"pos"}.$sep
    #   if (exists($search->{"pos"}) && $search->{"pos"});
    #$str .= $search->{"ref"}.$sep
    #   if (exists($search->{"ref"}) && $search->{"ref"});
    #$str .= join(",",@{$search->{"alt"}})
    #   if (exists($search->{"alt"}) && $search->{"alt"});
    
    if (!$str) {
        carp "\nNo search string given! Returning empty array...\n" ;
        return([]);
    }
    
    my $records = $self->get_records;
    my $r = first { $str == $_->{'_id'} } @{$records};
    
    return($r) if ($r);
    return([]);
}

=head2 output

Opens a file for writing the parsed/manipulated VCF records. Essentially a 
wrapper over Perl's C<open>. The second argument determines the writing mode, 
which is one of C<'new'> (default) or C<'append'>.

    $record->output('records.vcf');

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

=head2 write

VCF::Record writing to file. The C<output> method must be called first or a
valid filehandle must be provided.

    $record->write;

=cut

sub write {
    my ($self,$fh) = @_;
    
    # Even if $self->output has been called with a file handle and a 
    # file handle is provided, the second precedes. If $self->output has not
    # been called, and $fh is not a valid file handle output will be written
    # to STDOUT
    $self->{'_output'} = $fh if ($fh && fileno($fh) != -1);
    
    my $records = $self->get_records;
    
    # If file handle given, write to file otherwise to STDOUT
    my $fout = $self->{'_output'};
    
    # Format and write the records
    my $n = $self->length;
    for (my $i=0; $i<$n; $i++) {
        my $r = $records->[$i];
        my $line = $self->_format_record_line($r);
        ($fout) ? print $fout "$line\n" :
            print STDOUT "$line\n";
    }
}

=head2 done

Closes the file where VCF records are written. Essentially a  wraper over Perl's
C<close>.

    $record->done;

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

=head2 length

Returns the number of VCF records.

    my $size = $record->length;
    
=cut

sub length {
    my $self = $_[0];
    my $n = @{$self->get_records};
    return($n);
}

=head2 has_records

Check if a VCF::Record object has any records recorded. It will return a falsey
value if nothing has been parsed yet.

    my $has = $record->has_records;
    
=cut

sub has_records {
    my $self = $_[0];
    my $records = $self->get_records;
    return(1) if ($records && @{$records} > 0);
    return(0);
}

=head2 get_records

Returns a simple array of hashes with the VCF records and their fields.

    my $records = $record->get_records;
    
TODO: Add documentation (the field names)

=cut

sub get_records {
    my $self = $_[0];
    return([]) if (!$self->{"_content"});
    return($self->{"_content"});
}

=head2 get_record

Returns a hash with the $i-th VCF record, where $i can be a record we know or
the result of the C<find> method.

    my $rec = $record->get_record($i);
    
=cut

sub get_record {
    my ($self,$i) = @_;
    croak "\nNegative record index: $i!\n" if ($i < 0);
    return({}) if (!$self->{"_content"});
    return($self->{"_content"}->[$i]);
}

=head2 get

VCF::Header object getter

    my $param_value = $header->get("param_name");
    
=cut

sub get {
    my ($self,$name) = @_;
    return($self->{$name});
}

=head2 set_record

Replaces the $i-th record with the given record hash. Requires a VCF::Header 
object or the result of VCF::Header::get_header.
    
    $record->set_record($i,$rec,$header);
    
=cut

sub set_record {
    my ($self,$i,$rec,$header) = @_;
    my $len = $self->length;
    croak "\nNegative record index: $i!\n" if ($i < 0);
    croak "\nRecord index ($i) exceeds the number of records ($len)!\n" 
        if ($i > $len);
    croak "\nYou must provide the header for validation!\n" if ($header);
    croak "\nYou must provide a record!\n" if (!$rec);
    my $newrec = $self->validate_one($rec,$header);
    $self->{"_content"}->[$i] = $newrec;
    return($self);
}

=head2 set_records

Sets the actual parsed records content of the VCF file (an array of hashes)
either after parsing or constructed manually.
    
    my $file = 'my_variants.vcf';
    my $my_records = &parse_vcf_records($file);
    $record->set_records($my_records);
    
=cut

sub set_records {
    my ($self,$records) = @_;
    $self->{"_content"} = $records;
    return($self);
}

=head2 set

VCF::Header object setter

    $header->set("param_name","param_value");
    
=cut

sub set {
    my ($self,$name,$value) = @_;
    $self->{$name} = $value;
    return($self);
}

sub _init {
    my ($self,$params) = @_;
    
    $self->{'validation'} = $params->{'validation'};
    $self->{'file'} = $params->{'file'};
    
    if ($params->{'data'}) {
        $self->{'_content'} = [@{$params->{'data'}}]; # Copy
        undef($params->{'data'});
        #$self->validate_all;
    }
    else {
        $self->{'_content'} = undef;
    }
    
    $self->{'_output'} = undef;
    
    return($self);
}

sub _validate_one {
    my ($self,$fields,$header_ids,$sample_names) = @_;
    
    $fields->{"CHROM"} = $self->_validate_record_field_chrom($fields,
        $header_ids->{"contig"});
        
    $fields->{"POS"} = $self->_validate_record_field_pos($fields,$fields,
        $header_ids->{"contig"});
        
    $fields->{"ID"} = $self->_validate_record_field_id($fields);
    
    $fields->{"REF"} = $self->_validate_record_field_ref($fields);
    
    # If ALT present in header, we must put brackets around it as the ALT field
    # in VCF records has this format
    my @alt_ids = @{$header_ids->{"ALT"}};
    if (@alt_ids > 0) {
        for (my $i=0; $i<@alt_ids; $i++) {
            $alt_ids[$i] = "<".$alt_ids[$i].">";
        }
    }
    $fields->{"ALT"} = $self->_validate_record_field_alt($fields,$fields,
        \@alt_ids,$header_ids->{"fileformat"});
        
    $fields->{"QUAL"} = $self->_validate_record_field_qual($fields);
    
    $fields->{"FILTER"} = $self->_validate_record_field_filter($fields,
        $header_ids->{"FILTER"});
    
    $fields->{"INFO"} = $self->_validate_record_field_info($fields,$fields,
        $header_ids->{"INFO"});
            
    $fields->{"FORMAT"} = $self->_validate_record_field_format($fields,
        $header_ids->{"FORMAT"}) if (scalar @{$fields->{"FORMAT"}} > 0);
    
    $fields->{"SAMPLE"} = $self->_validate_record_field_sample($fields,
        $sample_names,$fields->{"FORMAT"},$fields->{"INFO"},
        $fields->{"ALT"},$header_ids->{"FORMAT"})
            if (scalar @{$fields->{"SAMPLE"}} > 0);
    
    # If there is no internal id (e.g. coming from manually constructed data,
    # an internal id should be defined)
    if (!exists($fields->{"_id"}) || !defined($fields->{"_id"})) {
        my $sep = "^^^"; # A strange separator to avoid conflict with data
        $fields->{"_id"} = $fields->{"CHROM"}.$sep.$fields->{"POS"}.$sep.
            $fields->{"REF"}.$sep.join(",",@{$fields->{"ALT"}});
    }
    
    return($fields);
}

sub _validate_record_field_chrom {
    my ($self,$chr,$contig,$ln) = @_;
    
    # Input to check from file or structure?
    $chr = $chr->{"CHROM"} if ($helper->_is_hash_ref($chr));
    
    # More informative error string if coming from file or structure
    my $minf = ($ln) ? ", line $ln:" : ":";
    
    # Validation level
    my $level = $self->get('validation');
    
    # Check if chromosome passed, will fail in case of consecutive tabs
    croak "\nCorrupted VCF #CHROM entry$minf #CHROM not found! Are there ".
        "consecutive tabs in the line?" if (!$chr);
    
    # Check the naming specs: no colon (:) or white space allowed
    if ($level eq "relaxed") {
        croak "\nCorrupted VCF #CHROM entry$minf #CHROM cannot contain colon ".
        "(:) or white space! Offending value: $chr" if ($chr =~ m/:|\s+/);
    }
    
    if ($level eq "strict") {
        # Check if contig validator is present in the header so we can validate 
        # the chromosome name, otherwise we don't check
        if (scalar keys(%{$contig}) > 0) {
            croak "\nCorrupted VCF #CHROM entry$minf: #CHROM not in allowed ".
                "contig names defined in the header! Offending value: $chr"
                    if (!$contig->{$chr});
        }
    }
    
    return($chr);
}

sub _validate_record_field_pos {
    my ($self,$pos,$chrom,$contig,$ln) = @_;
    
    # Input to check from file or structure?
    $pos = $pos->{"POS"} if ($helper->_is_hash_ref($pos));
    $chrom = $chrom->{"CHROM"} if ($helper->_is_hash_ref($chrom));
    
    # More informative error string if coming from file or structure
    my $minf = ($ln) ? ", line $ln:" : ":";
    
    # Validation level
    my $level = $self->get('validation');
    
    # Check if position passed, will fail in case of consecutive tabs
    croak "\nCorrupted VCF POS entry$minf POS not found! Are there ".
        "consecutive tabs in the line?" if (!$pos);
    
    # Check that it's an integer
    croak "\nCorrupted VCF POS entry$minf POS must be an integer! ". 
        "Offending value: $pos" if ($pos !~ /^\d+$/);
    
    croak "\nCorrupted VCF POS entry$minf POS must be 0 or positive! ".
        "Offending value: $pos" if ($pos < 0);
        
    # If pass above, it is
    $pos = int $pos;
    
    # Check if contig validator is present in the header so we can validate 
    # that position does not exceed length + 1, otherwise we don't check
    if ($level eq "strict") {
        if (scalar keys(%{$contig}) > 0) {
            croak "\nCorrupted VCF #CHROM entry$minf $pos exceeds the contig ".
                "length!" if ($pos > $contig->{$chrom} + 1);
        }
    }
    
    return($pos);
}

sub _validate_record_field_id {
    my ($self,$id,$ln) = @_;
    
    # Input to check from file or structure?
    $id = $id->{"ID"} if ($helper->_is_hash_ref($id));
    
    # More informative error string if coming from file or structure
    my $minf = ($ln) ? ", line $ln:" : ":";
    
    # Validation level
    my $level = $self->get('validation');
    
    # Check that it is not a null value in the sense of continuous tabs
    croak "\nCorrupted VCF ID entry$minf ID string or missing value not ".
        "found!" if (!$id);
    
    # Depending on hash or text input, $id may be a string or an array reference
    my @ids = ($helper->_is_array_ref($id)) ? (@{$id}) : split(";",$id);
    for (my $i=0; $i<@ids; $i++) {
        $ids[$i] =~ s/^\s+|\s+$//g;
        croak "\nCorrupted VCF ID entry$minf ID cannot contain white ".
            "space! Offending value: $ids[$i]" if ($ids[$i] =~ m/\s+/);
        if ($level eq "strict") {
            croak "\nCorrupted VCF ID entry$minf ID must be alphanumeric or . ".
                "only! Offending value: $ids[$i]" 
                    if ($ids[$i] !~ m/^[a-zA-Z0-9_\-.]*$/);
        }
    }
    
    return(\@ids);
}

sub _validate_record_field_ref {
    my ($self,$ref,$ln) = @_;
    
    # Input to check from file or structure?
    $ref = $ref->{"REF"} if ($helper->_is_hash_ref($ref));
    
    # More informative error string if coming from file or structure
    my $minf = ($ln) ? ", line $ln:" : ":";
    
    # Validation level
    my $level = $self->get('validation');
    
    # Check that it is not a null value in the sense of continuous tabs
    croak "\nCorrupted VCF REF entry$minf REF base or missing value not ".
        "found!" if (!$ref);
        
    # Check ref in A, C, G, T, N
    if ($level eq "relaxed" || $level eq "strict") {
        croak "\nCorrupted VCF REF entry$minf REF is not a valid nucleotide ".
            "representation! Offending value: $ref" if ($ref !~ m/^[ACGTN]+/i);
    }
    
    return($ref);
}

sub _validate_record_field_alt {
    my ($self,$alt,$ref,$alts,$vcfv,$ln) = @_;
    
    # Input to check from file or structure?
    $alt = $alt->{"ALT"} if ($helper->_is_hash_ref($alt));
    $ref = $ref->{"REF"} if ($helper->_is_hash_ref($ref));
    
    # More informative error string if coming from file or structure
    my $minf = ($ln) ? ", line $ln:" : ":";
    
    # Validation level
    my $level = $self->get('validation');
    
    # Check that it is not a null value in the sense of continuous tabs
    croak "\nCorrupted VCF ALT entry$minf ALT base(s) or missing value ".
        "not found!" if (!$alt);
    
    # Depending on hash or text input, $alt may be a string or an array
    # reference
    my @altes = ($helper->_is_array_ref($alt)) ? (@{$alt}) : split(",",$alt);
    
    # Get validation tools
    my $expr = $self->_init_vcf_field_expr($vcfv);
    my $altex = $expr->{"ALT"};
    
    # If the ALT field is present in the VCF header, it must be taken into
    # account
    if (@{$alts} > 0) {
        foreach my $a (@altes) {
            croak "\nCorrupted VCF ALT entry$minf ALT is not a valid ".
                "nucleotide or other valid representation! Offending value: $a" 
                    if ($a !~ m/$altex/ && !(grep {$_ eq $a} @{$alts}));
        }
    }
    else { # Easier
        foreach my $a (@altes) {
            croak "\nCorrupted VCF ALT entry$minf ALT is not a valid ".
                "nucleotide or other valid representation! Offending value: $a" 
                    if ($a !~ m/$altex/);
        }
    }
    
    # FIXME: The problem below
    # Check that first base of each allele must match the reference if their 
    # lengths are different. This may cause a problem with freebayes complex
    # alleles... Will issue a warning upon definition of validation strictness
    #if ($level eq "strict") {
    #   foreach my $a (@altes) {
    #       if (length($ref) != length($a)) {
    #           my $fr = substr($ref,0,1);
    #           my $fa = substr($a,0,1);
    #           die "\nCorrupted VCF ALT entry$minf ALT must start with the ".
    #               "reference base when ALT and REF not of the same length! ".
    #               "Offending value: $a" if ($fr ne $fa);
    #       }
    #    }
    #}
    
    # Check that the alleles are not the same as the reference
    foreach my $a (@altes) {
        croak "\nCorrupted VCF ALT entry$minf ALT cannot be the same as ".
            "REF! Offending value: $a" if ($a eq $ref);
    }
    
    return(\@altes);
}

sub _validate_record_field_qual {
    my ($self,$qual,$ln) = @_;
    
    # Input to check from file or structure?
    $qual = $qual->{"QUAL"} if ($helper->_is_hash_ref($qual));
    
    # More informative error string if coming from file or structure
    my $minf = ($ln) ? ", line $ln:" : ":";
    
    # Validation level
    my $level = $self->get('validation');
    
    # Check if qual passed, will fail in case of consecutive tabs
    croak "\nCorrupted VCF QUAL entry$minf QUAL not found! Are there ".
        "consecutive tabs in the line?" if (!$qual);
    
    # Check that it's an integer or float or "."
    if ($level eq "relaxed" || $level eq "strict") {
        croak "\nCorrupted VCF QUAL entry$minf QUAL must be a positive ".
            "integer or float or \".\" ! Offending value: $qual" 
                if ($qual !~ /^\d+$/ 
                    && $qual !~ /^(?=\d|\.\d)\d*(\.\d*)?([Ee][+-]?(\d+))?$/
                    && $qual ne ".");
    }
    
    return($qual);
}

sub _validate_record_field_filter {
    my ($self,$filter,$filts,$ln) = @_;
    
    # Input to check from file or structure?
    $filter = $filter->{"FILTER"} if ($helper->_is_hash_ref($filter));
    
    # More informative error string if coming from file or structure
    my $minf = ($ln) ? ", line $ln:" : ":";
    
    # Validation level
    my $level = $self->get('validation');
    
    # Check if filter passed, will fail in case of consecutive tabs
    croak "\nCorrupted VCF FILTER entry$minf FILTER not found! Are there ".
        "consecutive tabs in the line?" if (!$filter);
    
    # Check that it's PASS or ".", otherwise it should be a string separated
    # by ";" and the components must be in the filters specfied by the header
    if ($filter eq "PASS" || $filter eq ".") {
        return([$filter]);
    }
    elsif ($filter->[0] eq ".") {
        return($filter);
    }
    else {
        # Depending on hash or text input, $alt may be a string or an array
        # reference
        my @app_filts = ($helper->_is_array_ref($filter)) ? (@{$filter}) : 
            (split(";",$filter));
        foreach my $f (@app_filts) {
            if ($level eq "relaxed" || $level eq "strict") {
                croak "\nCorrupted VCF FILTER entry$minf FILTER not found in ".
                    "the filters specified in the header! Offending value: $f"
                        if (!(grep {$_ eq $f} @{$filts}));
            }
        }
        return(\@app_filts);
    }
}

sub _validate_record_field_info {
    my ($self,$info,$alt,$infos,$ln) = @_;
    
    # Input to check from file or structure?
    $info = $info->{"INFO"} if ($helper->_is_hash_ref($info));
    $alt = $alt->{"ALT"} if ($helper->_is_hash_ref($alt));
    
    # More informative error string if coming from file or structure
    my $minf = ($ln) ? ", line $ln:" : ":";
    
    # Validation level
    my $level = $self->get('validation');
    
    # Check if info passed, will fail in case of consecutive tabs
    croak "\nCorrupted VCF INFO entry$minf INFO not found! Are there ".
        "consecutive tabs in the line?" if (!$info);
    
    my @header_keys;
    if ($helper->_is_hash_ref($info)) {
        #tie %info_hash,"Tie::IxHash::Easy";
        #%info_hash = %{$info};
        if ($level eq "relaxed" || $level eq "strict") {
            @header_keys = keys(%{$infos});
            foreach my $key (keys(%{$info})) {
                my $value = $info->{$key};
                # Check if the key is defined in the header
                croak "\nCorrupted VCF INFO entry$minf key $key not declared ".
                    "in the VCF header!" if (!(grep {$_ eq $key} @header_keys));
                if (defined($value)) { # Can be a Flag
                    my @vals = @{$value};
                    if ($infos->{$key}->{"Number"} eq "A") {
                        my $card = @vals;
                        my $expc = @{$alt};
                        croak "\nCorrupted VCF INFO entry$minf $expc values ".
                            "expected for key $key, found $card! Offending ".
                            "value: ".$helper->_print_array_ref($value)
                                if ($card != $expc);
                                
                        # and check also the type
                        my $exp_type = $infos->{$key}->{"Type"};
                        foreach my $r (@vals) {
                            croak "\nCorrupted VCF INFO entry$minf expected ".
                                "$exp_type for key $key, found $r! Offending ".
                                "value: ".$helper->_print_array_ref($value)
                                    if (!$self->_check_allowed_type(
                                        $r,$exp_type));
                        }
                    }
                    # If G
                    elsif ($infos->{$key}->{"Number"} eq "G") {
                        # TODO: INFO validation - Stub for genotypes
                    }
                    # If R, the rights length must be equal to the number of all 
                    # alleles
                    elsif ($infos->{$key}->{"Number"} eq "R") {
                        my $card = @vals;
                        my $expc = @{$alt} + 1;                    
                        croak "\nCorrupted VCF INFO entry$minf $expc values ".
                            "expected for key $key, found $card! Offending ".
                            "value: ".$helper->_print_array_ref($value)
                                if ($card != $expc);
                        # and check also the type
                        my $exp_type = $infos->{$key}->{"Type"};
                        foreach my $r (@vals) {
                            croak "\nCorrupted VCF INFO entry$minf expected ".
                                "$exp_type for key $key, found $r! Offending ".
                                "value: ".$helper->_print_array_ref($value)
                                    if (!$self->_check_allowed_type(
                                        $r,$exp_type));
                        }
                    }
                    # If . then freedom, only Type checking
                    elsif ($infos->{$key}->{"Number"} eq ".") {
                        my $exp_type = $infos->{$key}->{"Type"};
                        foreach my $r (@vals) {
                            croak "\nCorrupted VCF INFO entry$minf expected ".
                                "$exp_type for key $key, found $r! Offending ".
                                "value: ".$helper->_print_array_ref($value) 
                                    if (!$self->_check_allowed_type(
                                        $r,$exp_type));
                        }
                    }
                    # Cardinality checking
                    else {
                        my $info_number = int($infos->{$key}->{"Number"});
                        my $card = @vals;                    
                        croak "\nCorrupted VCF INFO entry$minf $info_number ".
                            "values expected for key $key, found $card! ".
                            "Offending value: ".
                            $helper->_print_array_ref($value) 
                                if ($card != $info_number);
                        
                        # and check also the type
                        my $exp_type = $infos->{$key}->{"Type"};
                        foreach my $r (@vals) {
                            croak "\nCorrupted VCF INFO entry$minf expected ".
                                "$exp_type for key $key, found $r! Offending ".
                                "value: ".$helper->_print_array_ref($value) 
                                    if (!$self->_check_allowed_type(
                                        $r,$exp_type));
                        }
                    }
                }
                else { # Its header INFO Type must be Flag and Number 0
                    croak "\nCorrupted VCF INFO entry$minf no values found ".
                        "for key $key, but its header declared type is not ".
                        "Flag!" if ($infos->{$key}->{"Number"} ne 0 ||
                            $infos->{$key}->{"Type"} ne "Flag");
                }
            }
        }
        return($info);
    }
    else {
        # For each key-value pair, split per "=" to construct the hash of INFO 
        # data
        my %info_hash;
        tie %info_hash,"Tie::IxHash::Easy";
        my @info_data = split(";",$info);
        @header_keys = keys(%{$infos}) 
            if ($level eq "relaxed" || $level eq "strict");
        foreach my $kv (@info_data) {
            $kv =~ s/^\s+|\s+$//g;
            my ($key,$value) = split("=",$kv);
            # Check if the left side of "=" belongs to the header declared 
            # keys
            if ($level eq "relaxed" || $level eq "strict") {
                croak "\nCorrupted VCF INFO entry$minf key $key not declared ".
                    "in the VCF header!" if (!(grep {$_ eq $key} @header_keys));
            }
            # Validate cardinality of right side of "="
            my @vals = ();
            if (defined($value)) { # Can be a Flag
                @vals = split(",",$value);
                # If A, the rights length must be equal to the number of alt 
                # alleles
                if ($level eq "relaxed" || $level eq "strict") {
                    if ($infos->{$key}->{"Number"} eq "A") {
                        my $card = @vals;
                        my $expc = @{$alt};
                        croak "\nCorrupted VCF INFO entry$minf $expc values ".
                            "expected for key $key, found $card! Offending ".
                            "value: $value" if ($card != $expc);
                            
                        # and check also the type
                        my $exp_type = $infos->{$key}->{"Type"};
                        foreach my $r (@vals) {
                            croak "\nCorrupted VCF INFO entry$minf expected ".
                                "$exp_type for key $key, found $r! Offending ".
                                "value: $value"
                                    if (!$self->_check_allowed_type(
                                        $r,$exp_type));
                        }
                    }
                    # If G
                    elsif ($infos->{$key}->{"Number"} eq "G") {
                        # TODO: INFO validation - Stub for genotypes
                    }
                    # If R, the rights length must be equal to the number of all 
                    # alleles
                    elsif ($infos->{$key}->{"Number"} eq "R") {                     
                        my $card = @vals;
                        my $expc = @{$alt} + 1;
                        croak "\nCorrupted VCF INFO entry$minf $expc values ".
                            "expected for key $key, found $card! Offending ".
                            "value: $value" if ($card != $expc);
                        
                        # and check also the type
                        my $exp_type = $infos->{$key}->{"Type"};
                        foreach my $r (@vals) {
                            croak "\nCorrupted VCF INFO entry$minf expected ".
                                "$exp_type for key $key, found $r! Offending ".
                                "value: $value" 
                                    if (!$self->_check_allowed_type(
                                        $r,$exp_type));
                        }
                    }
                    # If . then freedom, only Type checking
                    elsif ($infos->{$key}->{"Number"} eq ".") {                     
                        my $exp_type = $infos->{$key}->{"Type"};
                        foreach my $r (@vals) {
                            croak "\nCorrupted VCF INFO entry$minf expected ".
                                "$exp_type for key $key, found $r! Offending ".
                                "value: $value" 
                                    if (!$self->_check_allowed_type(
                                        $r,$exp_type));
                        }
                    }
                    # Cardinality checking
                    else {
                        my $info_number = int($infos->{$key}->{"Number"});
                        my $card = @vals;
                        croak "\nCorrupted VCF INFO entry$minf $info_number ".
                            "values expected for key $key, found $card! ".
                            "Offending value: $value" 
                                if ($card != $info_number);
                        
                        # and check also the type
                        my $exp_type = $infos->{$key}->{"Type"};
                        foreach my $r (@vals) {
                            croak "\nCorrupted VCF INFO entry$minf expected ".
                                "$exp_type for key $key, found $r! Offending ".
                                "value: $value" 
                                    if (!$self->_check_allowed_type(
                                        $r,$exp_type));
                        }
                    }
                }
                else { # Its header INFO Type must be Flag and Number 0
                    if ($level eq "relaxed" || $level eq "strict") {
                        croak "\nCorrupted VCF INFO entry$minf no values ".
                            "found for key $key, but its header declared type ".
                            "is not Flag!" 
                                if ($infos->{$key}->{"Number"} ne 0 
                                    || $infos->{$key}->{"Type"} ne "Flag");
                    }
                }
            }
            # If not dead, add to the info hash
            @{$info_hash{$key}} = @vals;
        }
        return(\%info_hash);
    }
    
    #return(\%info_hash);
    
    # 1. Split per ";"
    # 2. For each member, split per "=" to construct the hash of INFO data
    # 3. Gather the lefts from (2) and check if they are in the header
    # 4. Based on $infos, check if the number of rights from (2)
    # 4.1. If 1, 2 etc, the rights length must be that number
    # 4.2. If A, the rights length must be equal to the number of alleles
    # 4.3. If G, not directly applicable here... Maybe another function that
    #      jointly checks INFO and FORMAT
    # 4.4. Version specific: if R, the rights length must be equal to the number
    #      of alleles + 1
    # 5. Based on $infos, check the type of each of the rights from (2) 
    #    (Integer, Float, String etc.). Maybe this can run in parallel with 4.
    #
    # The return value is an array of hashes or a simple hash? The 1st one is
    # more JSON and MongoDB like, the 2nd is more practical for later parsing
    # and annotation... We go with the 2nd since the INFO key names are
    # carefully designed to be JSON friendly.
}

sub _validate_record_field_format {
    my ($self,$format,$formats,$ln) = @_;
    
    # Input to check from file or structure?
    $format = $format->{"FORMAT"} if ($helper->_is_hash_ref($format));
    
    # More informative error string if coming from file or structure
    my $minf = ($ln) ? ", line $ln:" : ":";
    
    # Validation level
    my $level = $self->get('validation');
    
    # Check if info passed, will fail in case of consecutive tabs
    croak "\nCorrupted VCF FORMAT entry$minf FORMAT not found! Are there ".
        "consecutive tabs in the line?" if (!$format);
        
    # Depending on hash or text input, $id may be a string or an array reference
    my @format_keys = ($helper->_is_array_ref($format)) ? (@{$format}) : 
        split(":",$format);
    
    # Get validation keys
    my @header_keys;
    
    # The first key must be the GT field if it is present
    if ($level eq "relaxed" || $level eq "strict") {
        @header_keys = keys(%{$formats});
        croak "\nCorrupted VCF FORMAT entry$minf first key $format_keys[0] ".
            "is not GT and GT is present in the VCF header!" 
                if ($format_keys[0] ne "GT");
                
        # Check if all other format keys are present in the header
        foreach my $k (@format_keys) {
            $k =~ s/^\s+|\s+$//g;
            croak "\nCorrupted VCF FORMAT entry$minf key $k not declared ".
                "in the VCF header!" if (!(grep {$_ eq $k} @header_keys));
        }
    }
    
    # If not dead, return
    return(\@format_keys);
}

sub _validate_record_field_sample {
    my ($self,$sample_string,$sample_names,$format,$info,$alt,$formats,$ln) = 
        @_;
    
    # $sample_string: joined VCF fields corresponding to samples (string)
    # $sample_names: parsed SAMPLE names (array)
    # $format: the parsed FORMAT fields (array)
    # $info: INFO fields for genotype parsing
    # $alt: ALT fields for number of genotypes validation
    # $formats: FORMAT data types from the header
    # $ln: current line for verbosity
    
    # If no sample names parsed, return
    #return([]) if (!$sample_string || scalar @{$sample_names} == 0);
    return([]) if (!$sample_string);
    # If no sample names passed, then give an illegal one for 2-step parsing
    # and validation, otherwise, will fail
    my $foo_name = "Not provided";
    my $samples_passed = scalar @{$sample_names} != 0 ? 1 : 0;
    
    $sample_string = $sample_string->{"SAMPLE"}
        if ($helper->_is_hash_ref($sample_string));
    
    # Initiate the array of hashes
    my @sample_array;
    
    # Have the sample names handy
    my @snames = @{$sample_names};
    
    # More informative error string if coming from file or structure
    my $minf = ($ln) ? ", line $ln:" : ":";
    
    # Validation level
    my $level = $self->get('validation');
    
    # Have the header formats handy, sample format keys to be validated must
    # belong there
    my @fmts = @{$format};
    
    # Init the hash format
    my %format_hash;
    
    # In the structure case, if we are adding a key-value pair to a format, it  
    # must be given for ALL samples, thus it must be an array of hashes and each
    # hash MUST have a valid sample name in the "name" field. The order will be 
    # automatically added.
    if ($helper->_is_array_ref($sample_string)) { # From structure
        @sample_array = @{$sample_string};
        
        # Die if not the sample array of values to add to be validated does not
        # have all the samples
        if ($level eq "relaxed" || $level eq "strict") {
            croak "\nThe sample FORMAT key-value pairs provided does not ".
                "equal to the number of sample names provided!"
                    if (scalar @sample_array != scalar @{$sample_names});
        }
        
        # Foreach key, validate the format in each sample, using the header 
        # data, similarly to the INFO case
        if ($level eq "relaxed" || $level eq "strict") {
            for (my $i=0; $i<scalar @sample_array; $i++) {
                my $s = $sample_array[$i];
                
                # Die if there is no "name" field or if the name does not belong
                # to the sample names
                croak "No sample name provided for FORMAT key-value pairs $i!"
                    if (!defined($s->{"name"}) || !$s->{"name"});
                
                # Case of parsing without validation
                $s->{"name"} = $sample_names->[$s->{"order"}]
                    if ($s->{"name"} eq $foo_name);
                
                my $tmp = $s->{"name"};
                croak "Sample name $tmp provided with key-value pairs not ".
                    "found in provided sample names!"
                        if (!(grep {$_ eq $s->{"name"}} @{$sample_names}));
                    
                # Die if any key not already in the FORMAT column 
                # Allow our "name" and "order" custom attributes
                while (my($key,$value) = each(%{$s})) {                
                    croak "\nFormat key $key not found in the already ".
                        "declared formats!" if (!(grep {$_ eq $key} @fmts) 
                        && $key ne "name" && $key ne "order");
                    
                    next if ($key eq "name" || $key eq "order");
                    my @vals = @{$value};
                    
                    if ($value ne "") {
                        # Provided genotype is properly formed
                        if ($key eq "GT") {
                            croak "\nCorrupted VCF SAMPLE entry$minf ".
                                "malformed genotype! Offending value: ".
                                $helper->_print_array_ref($value) 
                                    if ($value->[0] !~ m/^(\.|\d+)([\|\/])?/);
                        }
                        
                        # If A
                        if ($formats->{$key}->{"Number"} eq "A") {
                            my $card = @vals;
                            my $expc = @{$alt};                        
                            croak "\nCorrupted VCF SAMPLE entry$minf $expc ".
                                "values expected for key $key, found $card! ".
                                "Offending value: ".
                                $helper->_print_array_ref($value)
                                    if ($card != $expc);
                                    
                            # and check also the type
                            my $exp_type = $formats->{$key}->{"Type"};
                            foreach my $r (@vals) {
                                croak "\nCorrupted VCF SAMPLE entry$minf ".
                                    "expected $exp_type for key $key, found $r!"
                                        if (!$self->_check_allowed_type(
                                            $r,$exp_type));
                            }
                        }
                        # If G
                        elsif ($formats->{$key}->{"Number"} eq "G") {
                            # TODO: INFO validation - Stub for genotypes
                        }
                        # If R, the vals length must be equal to number of all 
                        # alleles
                        elsif ($formats->{$key}->{"Number"} eq "R") {
                            my $card = @vals;
                            my $expc = @{$alt} + 1;
                            croak "\nCorrupted VCF SAMPLE entry$minf $expc ".
                                "values expected for key $key, found $card! ".
                                "Offending value: ".
                                $helper->_print_array_ref($value)
                                    if ($card != $expc);
                            
                            # and check also the type
                            my $exp_type = $formats->{$key}->{"Type"};
                            foreach my $r (@vals) {
                                croak "\nCorrupted VCF SAMPLE entry$minf ".
                                    "expected $exp_type for key $key, found $r!"
                                        if (!$self->_check_allowed_type(
                                            $r,$exp_type));
                            }
                        }
                        # If . then freedom, only Type checking
                        elsif ($formats->{$key}->{"Number"} eq ".") {
                            my $exp_type = $formats->{$key}->{"Type"};
                            foreach my $r (@vals) {
                                croak "\nCorrupted VCF SAMPLE entry$minf ".
                                    "expected $exp_type for key $key, found $r!"
                                        if (!$self->_check_allowed_type($r,
                                            $exp_type));
                            }
                        }
                        # Cardinality checking
                        else {
                            my $format_num = int($formats->{$key}->{"Number"});
                            my $card = @vals;
                            croak "\nCorrupted VCF SAMPLE entry$minf ".
                                "$format_num values expected for key $key, ".
                                "found $card! Offending value: ".
                                $helper->_print_array_ref($value)
                                    if ($card != $format_num);
                            
                            # and check also the type
                            my $exp_type = $formats->{$key}->{"Type"};
                            foreach my $r (@vals) {
                                croak "\nCorrupted VCF SAMPLE entry$minf ".
                                    "expected $exp_type for key $key, found $r!"
                                        if (!$self->_check_allowed_type($r,
                                            $exp_type));
                            }
                        }
                    }            
                    else { # Its header SAMPLE Type must be Flag and Number 0
                        croak "\nCorrupted VCF SAMPLE entry$minf no values ".
                            "found for key $key, but its header declared type ".
                            "is not Flag!" 
                                if ($formats->{$key}->{"Number"} ne 0 ||
                                    $formats->{$key}->{"Type"} ne "Flag");
                    }
                }
                # If not dead yet, add the order
                $s->{"order"} = $i;
            }
        }
    }
    else { # From file
        # Tie hash
        tie %format_hash,"Tie::IxHash::Easy";
        
        # Split samples
        my @samples = split("\t",$sample_string);
        # Foreach key, validate the format in each sample, using the header 
        # data, similarly to the INFO case
        for (my $i=0; $i<@samples; $i++) {
            my $s = $samples[$i];
                        
            # Add the sample name to the sample hash
            $format_hash{"name"} = $samples_passed ? $snames[$i] : $foo_name;
            $format_hash{"order"} = $i;
            
            # Start working with genotype data
            my @genotypes = split(":",$s);
            
            # Kill if genotype fields number is not the same as the format keys
            croak "\nCorrupted VCF SAMPLE entry$minf the number of genotype ".
                "fields is not the same as the declared FORMAT keys!" 
                    if (scalar @genotypes != scalar @fmts);
            
            # Validate sample genotype fields according to their header data 
            # types. No way to check if same order, so we go with the FORMAT 
            # order. Keep in mind that each Number type may contain the missing 
            # value "."
            for (my $j=0; $j<@fmts; $j++) {
                my $key = $fmts[$j];
                my $g = $genotypes[$j];
                my @vals = split(",",$g);
                
                if ($level eq "relaxed" || $level eq "strict") {
                    if ($g ne "") {
                        # Genotype is properly formed
                        if ($key eq "GT") {
                            croak "\nCorrupted VCF SAMPLE entry$minf ".
                                "malformed genotype! Offending value: $g" 
                                    if ($g !~ m/^(\.|\d+)([\|\/])?/);
                        }
                        
                        # If A
                        if ($formats->{$key}->{"Number"} eq "A") {
                            my $card = @vals;
                            my $expc = @{$alt};                        
                            croak "\nCorrupted VCF SAMPLE entry$minf $expc ".
                                "values expected for key $key, found $card! ".
                                "Offending value: $g" if ($card != $expc);
                                    
                            # and check also the type
                            my $exp_type = $formats->{$key}->{"Type"};
                            foreach my $r (@vals) {
                                croak "\nCorrupted VCF SAMPLE entry$minf ".
                                    "expected $exp_type for key $key, found ".
                                    "$r! Offending value: $g"
                                        if (!$self->_check_allowed_type($r,
                                        $exp_type));
                            }
                        }
                        # If G
                        elsif ($formats->{$key}->{"Number"} eq "G") {
                            # TODO: INFO validation - Stub for genotypes
                        }
                        # If R, the vals length must be equal to number of all
                        # alleles
                        elsif ($formats->{$key}->{"Number"} eq "R") {
                            my $card = @vals;
                            my $expc = @{$alt} + 1;
                            croak "\nCorrupted VCF SAMPLE entry$minf $expc ".
                                "values expected for key $key, found $card! ".
                                "Offending value: $g" if ($card != $expc);
                            
                            # and check also the type
                            my $exp_type = $formats->{$key}->{"Type"};
                            foreach my $r (@vals) {
                                croak "\nCorrupted VCF SAMPLE entry$minf ".
                                    "expected $exp_type for key $key, found ".
                                    "$r! Offending value: $g"
                                        if (!$self->_check_allowed_type($r,
                                            $exp_type));
                            }
                        }
                        # If . then freedom, only Type checking
                        elsif ($formats->{$key}->{"Number"} eq ".") {
                            my $exp_type = $formats->{$key}->{"Type"};
                            foreach my $r (@vals) {
                                croak "\nCorrupted VCF SAMPLE entry$minf ".
                                    "expected $exp_type for key $key, found ".
                                    "$r! Offending value: $g"
                                        if (!$self->_check_allowed_type($r,
                                            $exp_type));
                            }
                        }
                        # Cardinality checking
                        else {
                            my $format_num = int($formats->{$key}->{"Number"});
                            my $card = @vals;
                            croak "\nCorrupted VCF SAMPLE entry$minf ".
                                "$format_num values expected for key $key, ".
                                "found $card! Offending value: $g" 
                                    if ($card != $format_num);
                            
                            # and check also the type
                            my $exp_type = $formats->{$key}->{"Type"};
                            foreach my $r (@vals) {
                                croak "\nCorrupted VCF SAMPLE entry$minf ".
                                    "expected $exp_type for key $key, found ".
                                    "$r! Offending value: $g"
                                        if (!$self->_check_allowed_type(
                                            $r,$exp_type));
                            }
                        }
                    }        
                    else { # Its header SAMPLE Type must be Flag and Number 0
                        croak "\nCorrupted VCF SAMPLE entry$minf no values ".
                            "found for key $key, but its header declared type ".
                            "is not Flag!" if ($formats->{$key}->{"Number"} ne 0 
                                || $formats->{$key}->{"Type"} ne "Flag");
                    }
                }
                # If not dead, add to the info hash
                @{$format_hash{$key}} = @vals;
            }
            push(@sample_array,\%format_hash);
        }
    }
    return(\@sample_array);
}

sub _format_record_line {
    my ($self,$r) = @_;
    
    # Format #CHROM
    my $chrom = $r->{"CHROM"};
    # Format POS
    my $pos = $r->{"POS"};
    # Format ID
    my $id = join(",",@{$r->{"ID"}});
    # Format REF
    my $ref = $r->{"REF"};
    # Format ALT
    my $alt = join(",",@{$r->{"ALT"}});
    # Format QUAL
    my $qual = $r->{"QUAL"};
    # Format FILTER
    my $filter = join(",",@{$r->{"FILTER"}});
    # Format INFO
    my @info_contents;
    my $info_tmp = $r->{"INFO"};
    foreach my $key (keys(%{$info_tmp})) {
        my $val = join(",",@{$info_tmp->{$key}});
        push(@info_contents,$key."=".$val);
    }
    my $info = join(";",@info_contents);
    # Format FORMAT
    my $format = join(":",@{$r->{"FORMAT"}});
    # Format SAMPLE...
    my @samples;
    my @sample_tmp = @{$r->{"SAMPLE"}};
    my @sample_sorted =  sort { $a->{"order"} <=> $b->{"order"} } @sample_tmp;
    foreach my $s (@sample_sorted) {
        my @genotypes;
        foreach my $f (@{$r->{"FORMAT"}}) {
            push(@genotypes,join(",",@{$s->{$f}}));
        }
        push (@samples,join(":",@genotypes));
    }
    my $jsamples = join("\t",@samples);
    
    # Final line
    my $line = "$chrom\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format".
        "\t$jsamples";
    return($line);
}

sub _is_record_object {
    my $self = $_[0];
    return(1) if (defined($self->{"_content"}));
    return(0);
}

sub _extract_validation_data {
    my ($self,$header) = @_;
    
    my ($header_ids,$sample_names);
    if (ref($header) eq "VCF::Header") {
        $header_ids = $header->_extract_header_ids;
        $sample_names = $header->_extract_sample_names;
    }
    else {
        my $ho = VCF::Header->new({'data' => $header});
        $ho->set('validation',$self->get('validation'));
        $ho->validate;
        $header_ids = $ho->_extract_header_ids;
        $sample_names = $ho->_extract_sample_names;
    }
    
    return($header_ids,$sample_names);
}

sub _get_sample_names_from_record {
    my ($self,$rec) = @_;
    my @names;
    foreach my $r (@{$rec->{"SAMPLE"}}) {
       push(@names,$r->{"name"});
    }
    return(\@names);
}

sub _check_unique_fields {
    my ($self,$fields,$ln) = @_;
    
    my @ids = @{$fields->{"ID"}};
    my $allele = $fields->{"CHROM"}."_".$fields->{"POS"}.$fields->{"REF"}.
        join(",",@{$fields->{"ALT"}});
    
    foreach my $id (@ids) {
        $GLOBAL_IDS{$id}++;
        return($id) if ($id ne "." && $GLOBAL_IDS{$id} > 1);
    }
    #FIXME: The unique allele check does not apply to structural variants
    #$meta_alleles{$allele}++;
    #die "\nCorrupted VCF ID entry, line $ln: ID $allele has been found more ".
    #    "than once!" if ($meta_alleles{$allele} > 1);
    
    return(0);
}

sub _init_vcf_field_expr {
    my ($self,$ver) = @_;
    my $specs = $self->_define_vcf_field_expr;
    return($specs->{$ver});
}

# Enriched with dbVar structural variation types!
sub _define_vcf_field_expr {
    my %specs = (
        "VCFv4.1" => {
            "ALT" => qr/^(\.)?((\]{1}(?![:\s]).+:\d+\]{1})|(\[{1}(?![:\s]).+:\d+\[{1}))?[acgtnACGTN]+((\]{1}(?![:\s]).+:\d+\]{1})|(\[{1}(?![:\s]).+:\d+\[{1}))?(\.)?$|^<(ALT|INS|DEL|DUP|INV|CNV|DUP:TANDEM|INS:NOVEL|INS:ME(:(ALU|L1|LINE1|SVA|HERV))?|DEL:ME(:(ALU|L1|LINE1|SVA|HERV))?)>$/
        },
        "VCFv4.2" => {
            "ALT" => qr/^(\.)?((\]{1}(?![:\s]).+:\d+\]{1})|(\[{1}(?![:\s]).+:\d+\[{1}))?[acgtnACGTN]+((\]{1}(?![:\s]).+:\d+\]{1})|(\[{1}(?![:\s]).+:\d+\[{1}))?$|^(\*)?(\.)?$|^<(ALT|INS|DEL|DUP|INV|CNV|DUP:TANDEM|INS:NOVEL|INS:ME(:(ALU|L1|LINE1|SVA|HERV))?|DEL:ME(:(ALU|L1|LINE1|SVA|HERV))?)>$/
        },
        "VCFv4.3" => {
            "ALT" =>  qr/^(\.)?((\]{1}(?![:\s]).+:\d+\]{1})|(\[{1}(?![:\s]).+:\d+\[{1}))?[acgtnryswkmbdhvnACGTNRYSWKMBDHV\.-]+(\.)?((\]{1}(?![:\s]).+:\d+\]{1})|(\[{1}(?![:\s]).+:\d+\[{1}))?$|^\*?$|^<(INS|DEL|DUP|INV|CNV|DUP:TANDEM|INS:NOVEL|INS:ME(:(ALU|L1|LINE1|SVA|HERV))?|DEL:ME(:(ALU|L1|LINE1|SVA|HERV))?)>$/,
            "INFO" => qr/^([A-Za-z ][0-9A-Za-z .]*|1000G)$/
        }
    );
    return(\%specs);
}

sub _check_allowed_type {
    my ($self,$val,$type) = @_;
    if ($type eq "Integer") {
        return(1) if ($val =~ m/^(\.)?|([+-]?\d+)$/);
    }
    elsif ($type eq "Float") {
        return(1) 
            if ($val =~ m/^(\.)?|([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee][+-]?(\d+))?$/);
    }
    elsif ($type eq "String") {
        return(1) if !ref($val);
    }
    elsif ($type eq "Character") {
        return(1) if ($val =~ m/^[A-Za-z\.]{1}$/);
    }
    return(0);
}

sub _check_params {
    my ($self,$params) = @_;
    
    # $params can be a hash reference with parameters or an existing filename
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
        undef($params->{'data'});
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

    perldoc VCF::Record


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

1; # End of VCF::Record
