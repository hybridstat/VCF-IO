package VCF::_Helper;

use 5.10.0;
use strict;
use warnings;
use utf8;

binmode STDOUT, ":utf8";

$|=1;

use Carp;

use Data::Dumper qw(Dumper);

=head1 NAME

VCF::_Helper - Helper class for the VCF::IO package!

=head1 VERSION

Version 0.01

=cut

our $MODNAME = "VCF::_Helper";
our $VERSION = '0.01';

=head1 SYNOPSIS

VCF::_Helper contains certain common function used throught the VCF::IO package.
They are not meant to be used outside the scope of the package.

    use VCF::_Helper;

    my $helper = VCF::_Helper->new;
    $helper->_smatch('foo',('foo','bar'));
    
=head1 SUBROUTINES/METHODS

=head2 new

The VCF::_Helper object constructor.

    my $helper = VCF::_Helper->new;

=cut

sub new {
    my $class = shift @_;
    my $self = {};
    bless($self,$class);
    return($self);
}

sub _disp {
    my $self = shift @_;
    my @msg = @_;
    print STDERR "\n@msg";
}

sub _smatch {
    my ($self,$s,@a) = @_;
    if (grep(/^$s$/,@a)) {
        return(1);
    }
    return(0);
}

sub _unique {
    my ($self,@list) = @_;
    my (%seen,$item);
    foreach $item (@list) {
        $seen{$item}++;
    }
    return(\%seen);
}

sub _is_hash_ref {
    my ($self,$h) = @_;
    return(1) if (ref($h) eq "HASH");
    return(0);
}

sub _is_array_ref {
    my ($self,$a) = @_;
    return(1) if (ref($a) eq "ARRAY");
    return(0);
}

sub _print_array_ref {
    my ($self,$a) = @_;
    return(join(",",@{$a}));
}

sub _print_array {
    my ($self,@a) = @_;
    return(join(",",@a));
}

=head1 AUTHOR

Panagiotis Moulos, C<< <pmoulos at hybridstat.gr> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-vcf-io at rt.cpan.org>, or through
the web interface at L<https://rt.cpan.org/NoAuth/ReportBug.html?Queue=VCF-IO>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc VCF::_Helper


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

1; # End of VCF::_Helper
