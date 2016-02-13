#! /usr/bin/perl

# Wed Apr  3 13:16:31 EDT 2013
# Author: sbhatia 
# You need wget for this script to work.

# Usage: ./pigget.pl <pigeon file name>

# Get the file from the command prompt.
$pigeonFile = $ARGV[0];

# Get it's contents into a variable.
$desc = `cat $pigeonFile`;

# Set up Pigeon's real URL in a variable for convenience.
$url = "http://synbiotools.bu.edu:5801";

# Post Pigeon program to the server and get the php file it generates.
$a = `wget -q --post-data="desc=$desc" $url/pigeon1.php $url/pigeon.php`;

# Load up the php file it generates into a variable.
$a = `cat pigeon.php`;

# Look for the name of the image file.
$a = `grep "Weyekin" pigeon.php`;
($fileName) = ($a =~ /src\s*=([0-9a-zA-Z\/\\\.\-]+.png)\w*/);

# Download the image file.
$a = `wget -q $url/$fileName -O $ARGV[1]`;

# Tell the user the name of the file downloaded.
print STDOUT "Image $fileName retrieved. Bye.\n";

# Clean previous files, if any.
$a = `rm pigeon.php* pigeon1.php*`
