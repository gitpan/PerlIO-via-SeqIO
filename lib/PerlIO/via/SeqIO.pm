#$Id: SeqIO.pm 513 2009-10-22 18:09:34Z maj $
# PerlIO layer for sequence format with BioPerl guts
# Enjoy!
package PerlIO::via::SeqIO;
use strict;
use warnings;
no warnings qw(once);
use Bio::Seq;
use Carp;
use IO::String;
use IO::Seekable;
use File::Temp qw(tempfile);
use Exporter;
our $VERSION = '0.021';
our @ISA = qw(Exporter);
our @EXPORT = qw(open);
our @EXPORT_OK = qw(O T);
our %OBJS;
our %ITERS;
our $INSTANCE = 128; # big "fileno"
BEGIN {
    # crazy setup machinations...
    use Bio::SeqIO;
    our @SUPPORTED_FORMATS = qw(embl fasta genbank gcg pir);
    foreach (@SUPPORTED_FORMATS) {
	no strict qw(refs);
	my $pkg = "PerlIO::via::SeqIO::$_";
	my $io_class = "Bio::SeqIO::$_";
	Bio::SeqIO->_load_format_module($_);
	$ITERS{$io_class} = \&{$io_class."::next_seq"};
	my $code =<<END;
package ${pkg};
use strict;
use warnings;
use base qw(PerlIO::via::SeqIO);
1;
END
        eval $code;
	croak($@) if $@;
	{ # hook
	    no warnings qw(redefine);
	    *{$io_class."::next_seq"} = sub {
		my $seq = $PerlIO::via::SeqIO::ITERS{$io_class}->(@_);
		return unless $seq;
		$PerlIO::via::SeqIO::OBJS{sprintf("%s\n",$seq)} = $seq;
		return $seq;
	    };
	}
    }
    return;
}

our %MODE_SYM = (
    'r' => '<',
    'w' => '>',
    'a' => '>>',
    'r+' => '+<',
    'w+' => '+>',
    'a+' => '+>>'
    );

# init the layer:

sub PUSHED {
    no strict qw(refs);
    ${'PUSHED'.$INSTANCE} = { format => (split m/::/, $_[0])[-1], mode => $MODE_SYM{$_[1]}, instance => $INSTANCE };
    bless ( \*{'PUSHED'.$INSTANCE}, $_[0] );
}

# define our own open

sub open(*;@) {
    no strict qw(refs);
    my ($out_handle, $mode, $file) = @_;
    my ($source, $engine, $in_handle, $fh, $dumh, $fmt);
    unless (defined $file or $mode =~ /via/) {
	($mode, $file) = ($mode =~ /^(\+?(?:<?|>>?|\|)?)(.+)$/);
	$mode ||= '<';
    }
#   passthrough to real open here
    if ($mode !~ m/:via/) {
	my $ret;
	if ( $file =~ s/^&// ) {
	    $file = (caller)[0]."::$file" unless $file =~ /::/;
	    $ret = CORE::open($fh, $mode.'&'.$file);
	}
	else {
	    $ret = CORE::open($fh, $mode, $file);
	}
	croak("$!") unless $ret;
	if ( $_[0] and $_[0] =~ /^[A-Z]+$/ ) {
	    # bareword handle
	    *{(caller)[0]."::$_[0]"} = *$fh;
	}
	$_[0] ||= $fh;
    }
    else {
	# deal with first-arg file handles
	if ( $out_handle and $out_handle =~ /^[A-Z]+$|::/ ) {
	    no strict qw(refs);
	    $out_handle = (caller)[0]."::$out_handle" unless $out_handle =~ /::/;
	    $out_handle = *{$out_handle};
	    ($dumh, $file) = tempfile unless $file;
	}
	# deal with third-arg file handles
	if ($file =~ /^&/ or
	    (! -e $file and $file =~ /^[A-Z]+$|::/)) {
	    no strict qw(refs);
	    $file =~ s/^&//;
	    # use the caller namespace unless $file appears
	    # to be fully qualified...
	    $file = (caller)[0]."::$file" unless $file =~ /::/;
	    # create the real glob
	    $in_handle = \*{$file}; 
	    # to kick the PUSHED routine with a CORE::open,
	    # need to provide a real (but dummy) file
	    ($dumh, $file) = tempfile;
	}
	elsif (ref($file) eq 'GLOB') {
	    $in_handle = $file;
	    # to kick the PUSHED routine with a CORE::open,
	    # need to provide a real (but dummy) file
	    ($dumh, $file) = tempfile;
	}
	$fmt = (split m/::/, __PACKAGE__)[-1];
	$fmt eq 'SeqIO' and $fmt = ($mode =~ /.*?SeqIO::([a-zA-Z]+)/)[0];
	croak("Can't guess sequence format on write") 
	    if (!$fmt and $mode =~ />/);
	tie($fh = *{'TFH'.$INSTANCE}, '_TH');
	# to kick the PUSHED method
	if ( $out_handle ) {
	    if ($out_handle =~ /(DATA|STD(?:IN|OUT))/) {
		# special...don't reopen, but pass output thru
		# $fh first...
		my $bareword = $1;
		my $redirect = ($bareword =~ /OUT/ ? ">&" : "<&");
		# get a pristine copy of DATA/STDIN/STDOUT
		CORE::open(my $dup, $redirect, (caller)[0]."::$bareword");
		if ($bareword =~ /OUT/) {
		    $| && $dup->autoflush(1);
		}
		else { #STDIN
                    # clean copy for Bio::SeqIO, unless already attached
		    # to a file
		    $in_handle = $dup unless $in_handle; 
		}
		CORE::open( $fh, $mode, $file) or croak( $! );
		# retie
		$fh = $out_handle;
		tie $fh, '_TH';
		(tied $fh)->fh($dup);
	    }
	    else {
		$fh = *$out_handle;
		# retie
		tie($fh, '_TH');
		CORE::open( $fh, $mode, $file) or croak( $! );
	    }
	}
	else {
	    CORE::open( $fh, $mode, $file) or croak($!);
	}
	my $o = O(\*{'PUSHED'.$INSTANCE});
	(tied $fh)->o($o);
	$$o{mode} = $mode;
	# reads and writes use different machines:
	if ($mode !~ />/) { #read
	    $$o{source} = Bio::SeqIO->new( ( $in_handle ? (-fh => $in_handle) : (-file => $file)), -format=>$fmt); #RO
	    $$o{format} eq 'SeqIO' and 
		$$o{format} = (split m/::/,ref($$o{source}))[-1];
	}
	else { #write
	    (tied $fh)->set_write_format($$o{format});
	    if ($in_handle) {
		# tell the tied handle 
		# where we really want output to go
		(tied $fh)->fh($in_handle);
	    }
	    
	}

	$$o{fh}=$fh;
	$$o{fname}=$file;
	$$o{instance}=$INSTANCE;
	$_[0] ||= $fh;
	$INSTANCE++;
	return 1;
    }
}

### the via class provides the hook into the PerIO mechanism

### the filehandle tie takes care of the reading and writing
### in particular, the tie provides the hook into readline, 
### allows control over what constitutes a "line", independent of
### the value of $/ or anything else.

### so the FILL is just a dummy, because its definition is 
### required by the via API.

sub FILL {
    my ($obj, $fh) = @_;
    return;
}

# seq object converter

sub T {
    my @objs = @_;
    my @ret;
    foreach my $s (@objs) {
	unless (defined $s and ref($s) and
		($s->isa('Bio::SeqI') || $s->isa('Bio::Seq') || 
		 $s->isa('Bio::PrimarySeq'))) {
	    carp "Item undefined or not a sequence object; returning an undef";
	    push @ret, undef;
	    next;
	}
	$s->isa('Bio::PrimarySeq') and $s = _pseq_to_seq($s);
	push @ret, sprintf("%s\n", $s);
	$OBJS{$ret[-1]} = $s;
    }
    return wantarray ? @ret : $ret[0];
}

# object getter ...

sub O {
    my $sym = shift;
    no strict qw(refs);
    if (ref($sym) =~ /via|GLOB/) {
	return ${*$sym{SCALAR}};
    }    
    elsif (!ref($sym)) {
	for ($sym) {
	    m/Bio/ && do {
		return $OBJS{$sym}; last;
	    };
	    m/via/ && do {
		return (tied $sym)->o; last;
	    };
	    m/^[A-Z]+$/ && do {
		$sym = (caller)[0]."::$sym";
		return *$sym{SCALAR};
	    };
	}
    }
    else {
	croak("Don't understand the arg");
    }
}

# wrap Bio::PrimarySeqs (incl. Bio::LocatableSeqs) in a Bio::Seq
# for Bio::SeqIO use

sub _pseq_to_seq {
    my @pseqs = @_;
    my @ret;
    foreach (@pseqs) {
	unless (defined $_ && ref($_) && $_->isa('Bio::PrimarySeq')) {
	    push @ret, undef;
	    next;
	}
	my $seq = Bio::Seq->new();
	$seq->id( $_->display_id || $_->id );
	$seq->primary_seq( $_ );
	push @ret, $seq;
    }
    return wantarray ? @ret : $ret[0];
}

1;

package _TH;
use strict;
use warnings;
no strict qw(refs);
use IO::Seekable;
push @IO::Handle::ISA, qw(IO::Seekable);
our $AUTOLOAD;

sub TIEHANDLE { bless({}, shift) }

sub READLINE {
    my ($self, @args) = @_;
    my $fh = $self->fh;
    my $o = $self->o;
    if ($o) { # do SeqIO 
	my $seq = $$o{source}->next_seq;
	return unless $seq;
	if (!wantarray) {
	    return sprintf("%s$/",$seq);
	} 
	else {
	    local $_;
	    my @ret = ($seq);
	    while ( $_ = $$o{source}->next_seq ) {
		push @ret, $_;
	    }
	    return map { sprintf("%s$/", $_) } @ret;
	}
    }
    else { # passthrough
	return (wantarray ? $fh->getlines(@args) : $fh->getline(@args));
    }
}

sub OPEN {
    my ($self, $mode,$file) = @_;
    my ($fh, $lfh);
    my $ret = CORE::open($fh, $mode,$file);
    if ($mode =~ /via/) {
	# get a lower level handle
	$mode =~ s/:via.*//;
	CORE::open($lfh, $mode, $file);
	$self->fh($lfh);
    }
    else {
	$self->fh($fh);
    }
    return $ret;
}

sub PRINT{
    my ($self, @args) = @_;
    my $o = $self->o;
    my $fh = $self->fh;
    my $ios = $$o{io_string};
    my $ret = 0;
    foreach (@args) {
	my @input = split( /(Bio::.*?=HASH\(0x[0-9a-f]+\)\s)/);
	foreach my $item (@input) {
	    if ($item =~ /HASH/) { # string rep of object
		$item = $OBJS{$item};
		$$o{engine}->write_seq($item);
		$ios->pos(0); # seek to top
		my $line = join('', <$ios>);
		$ios->pos(0); ${$ios->string_ref}='';
		$self->write($line, length $line);
		undef $OBJS{$item}; # clean up
	    }
	    else { # write the raw buffer
		$self->write($item, length);
	    }
	}
    }
    1;
}


sub FLUSH { shift->fh->flush }

sub CLOSE {
    my $self = shift;
    my $ret = $self->fh->close;
    return $ret;
}
# autoload the other filehandle methods

sub AUTOLOAD {
    my $func = lc((split(m/::/,$AUTOLOAD))[-1]);
    my $self = shift;
    $func =~ s/destroy/DESTROY/;
    $self->fh && $self->fh->$func(@_);
}

# on the fly format changes

sub set_write_format {
    my $self = shift;
    my $format = shift;
    my $o = $self->o;
    unless (grep( /^$format$/, @PerlIO::via::SeqIO::SUPPORTED_FORMATS)) {
	carp("The format '$format' isn't supported; current format unchanged");
	return;
    }
    unless ($$o{mode} && $$o{mode} =~ />|\+/) {
	carp("Can't set format; handle not open for writing");
	return;
    }
    $self->o->{format} = $format;
    $$o{io_string} = IO::String->new();
    $$o{engine} = Bio::SeqIO->new(-fh=>$$o{io_string},
				  -format=>$$o{format});
    $$o{engine}->_flush_on_write(1);
    return 1;
}

# provide the connection to the via instance

sub o {
    my ($self, $o) = @_;
    return $self->{o} = $o if (defined $o);
    return $self->{o};
} 

# store the "real" filehandle
   
sub fh {
    my ($self, $fh) = @_;
    return $self->{fh} = $fh if ($fh);
    return $self->{fh};
}

1;	

__END__

=pod

=head1 NAME

PerlIO::via::SeqIO - PerlIO layer for biological sequence formats

=head1 SYNOPSIS
 
 use PerlIO::via::SeqIO;

 # open a FASTA file for reading:
 open( my $f, "<:via(SeqIO)", 'my.fas');

 # open an EMBL file for writing
 open( my $e, ">:via(SeqIO::embl)", 'my.embl');

 # convert
 print $e $_ while (<$f>);

 # add comments (this really works)
 while (<$f>) {
   # get the real sequence object
   my $seq = O($_);
   if ($seq->desc =~ /Pongo/) {
     print $e "# this one is almost human...";
   }
   print $e $_; 
 }

 # a one-liner, sort of
 $ alias scvt="perl -Ilib -MPerlIO::via::SeqIO -e \"open(STDIN, '<:via(SeqIO)'); open(STDOUT, '>:via(SeqIO::'.shift().')'); while (<STDIN>) { print }\""
 $ cat my.fas | scvt gcg > my.gcg

=head1 DESCRIPTION

C<PerlIO::via::SeqIO> attempts to provide an easy option for
harnessing the magic sequence format I/O of the BioPerl
(L<http://bioperl.org>) toolkit. Opening a biological sequence file
under C<via(SeqIO)> yields a filehandle that can be used to read and
write L<Bio::Seq> objects sequentially with an absolute minimum of
setup code.

C<via(SeqIO)> also allows the user to mix plain text and sequence formats
on a single filehandle transparently. Different sequence formats
can be written to a single file by a simple filehandle tweak.

=head1 DETAILS

=over

=item Basics

Here's the basic idea, in code converting FASTA to EMBL format:

 open($in, '<:via(SeqIO)', 'my.fas');
 open($out, '>:via(SeqIO::embl)', 'my.embl');
 while (<$in>) {
   print $out $_;
 }

Scalar and bareword filehandles both are understood by C<via(SeqIO)>,
as well as STDIN, STDOUT, and DATA. For example:

 open(STDIN, '<:via(SeqIO)');
 ...

allows 

 cat my.gcg | perl your.pl > out

where C<your.pl> can read STDIN and acquire the sequence objects by
using the object getter L</UTILITIES/O()>. The format
of the input in this case will be guessed by the C<Bio::SeqIO>
machinery.

On reading, you can rely on L<Bio::SeqIO>'s format guesser by invoking
an unqualifed

 open($in, '<:via(SeqIO)', 'mystery.txt');

or you can specify the format, like so:

 open($in, '<:via(SeqIO::embl)', 'mystery.txt');

On writing, a qualified invocation is required;

 open($out, '>:via(SeqIO)', 'my.fas');        # throws
 open($out, '>:via(SeqIO::fasta)', 'my.fas'); # that's better

=item Retrieving the sequence object itself

This does what you mean:

 open($in, '<:via(SeqIO)', 'my.fas');
 open($out, '>:via(SeqIO::embl)', 'my.embl');
 while (<$in>) {
   print $out $_;
 }

However, C<$_> here is not the sequence object itself. To get that use 
the all-purpose object getter L</UTILITIES/O()>:

 while (<$in>) {
   print join("\t", O($_)->id, O($_)->desc), "\n";
 }

=item Writing a I<de novo> sequence object

Use the L</UTILITIES/T()> mapper to convert a Bio::Seq object into a thing that can be formatted by C<via(SeqIO)>:

 open($seqfh, ">:via(SeqIO::embl)", "my.embl");
 my $result = Bio::SearchIO->new( -file=>'my.blast' )->next_result;
 while(my $hit = $result->next_hit()){
   while(my $hsp = $hit->next_hsp()){
     my $aln = $hsp->get_aln;
       print $seqfh T($_) for ($aln->each_seq);
     }
   }

=item Writing plain text

Interspersing plain text among your sequences is easy; just print the
desired text to the handle. See the L</SYNOPSIS>.

Even the following works:

 open($in, "<:via(SeqIO)", 'my.fas')
 open($out, ">:via(SeqIO::embl)", 'annotated.txt');

 $seq = <$in>;
 print $out "In EMBL format, the sequence would be rendered:", $s;

=item Switching write formats

You can also easily switch write formats. (Why? Because...who knows?)
Use C<set_write_format> off the tied handle object:

 open($in, "<:via(SeqIO)", 'my.fas')
 open($out, ">:via(SeqIO::embl)", 'multi.txt');

 $seq1 = <$in>;
 print "This is sequence 1 in embl format:\n";
 print $out $seq1;
 (tied $out)->set_write_format(gcg);
 print $out "while this is sequence 1 in GCG format:\n"
 print $out $seq1;

=item Supported Formats

The supported formats are contained in
C<@PerlIO::via::SeqIO::SUPPORTED_FORMATS>. Currently they are

 fasta, embl, gcg, genbank

=back

=head1 IMPLEMENTATION

This is essentially a hack, but one that attempts to behave fairly 
well. The handles are highly overloaded, with one foot in PerlIO::via
and the other in C<tie>. Things to keep in mind:

=over

=item C<PerlIO::via::SeqIO> exports C<open()>

Neither L<PerlIO::via> nor C<tie> provided a low enough hook. When the
mode does not contain a C<:via()> call, your opens are passed through
to C<CORE::open>. If you run into problems, please ping me. "Why didn't
you do it like ..." comments are also most welcome.

=item Peeking at the guts

The filehandle takes notes in a hash under the hood. To look at it,
use the object getter:

 $o = O($fh);
 print join("\n", keys %$o), "\n";

The "public" interface (see L</UTILITIES>) is available thru the tied
object; that is

 (tied $fh)

and not 

 O($fh).

=back

=head1 UTILITIES

In the C<PerlIO::via::SeqIO> namespace. To use, do

 use PerlIO::via::SeqIO qw(open O T);

(The C<open> hook needs to be available for the package to
function. It is a member of C<@EXPORT>. See
L</IMPLEMENTATION/C<PerlIO::via::SeqIO> exports C<open()>>.)

=head2 O()

 Title   : O
 Usage   : $o = O($sym) # export it; not an object method
 Function: get the object "represented" by the argument
 Returns : the right object
 Args    : PerlIO::via::SeqIO GLOB, or 
           *PerlIO::via::SeqIO::TFH (tied fh) or
           scalar string (sprintf-rendered Bio::SeqI object)
 Example : $seqobj = O($s = <$seqfh>);

=head2 T()

 Title   : T
 Usage   : T($seqobj) # export it; not an object method
 Function: Transform a real Bio::Seq object to a
           via(SeqIO)-writeable thing
 Returns : A thing writeable as a formatted sequence
           by a via(SeqIO) filehandle
 Args    : a[n array of] Bio::Seq or related object[s]
 Example : print $seqfh T($seqobj);

=head2 set_write_format()

 Title   : set_write_format
 Usage   : (tied $fh)->set_write_format($format)
 Function: Set a write handle to write a specified 
           sequence format
 Returns : true on success
 Args    : scalar string; a supported format 
           (see @PerlIO::via::SeqIO::SUPPORTED_FORMATS)

=head1 SEE ALSO

L<perlio>, L<PerlIO::via>, L<Bio::SeqIO>, L<Bio::Seq>, 
L<http://bioperl.org>

=head1 AUTHOR - Mark A. Jensen

 Email maj -at- fortinbras -dot- us
 http://fortinbras.us
 http://bioperl.org/wiki/Mark_Jensen

=cut
