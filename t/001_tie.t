#-*-perl-*-
#$Id$
use strict;
use warnings;
use lib 'lib';
use Test::More tests => 9;
use File::Temp qw(tempfile);
use_ok('IO::Seekable');
use_ok('PerlIO::via::SeqIO');

use PerlIO::via::SeqIO qw(O open);

# secret tied handle class...
ok(my ($fh,$fn) = tempfile(), "get tempfile");
open(my $dup, "+<&", $fh); # an untied version 
ok( tie(*$fh,'_TH'), "tie handle");
(tied *$fh)->fh($dup); # hack
ok( print $fh "hey dude\n", "print");
ok( $fh->flush, "flush");
ok( $fh->seek(0,0), "seek to top" );
ok( <$fh> eq "hey dude\n", "what we just wrote");
ok( $fh->close, "close");

