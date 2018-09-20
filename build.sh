#!/bin/bash

if [ "$(id -u)" != "0" ]; then
	echo "You must execute this script as root!"
else
	# Try to install required Perl packages in case they are missing
	perl -MCPAN -e 'install qw(IO::Uncompress::Gunzip List::Util 
        Tie::IxHash::Easy Test::Exception Test::File)'

	perl Makefile.PL
	make
	make test
	make install
fi
