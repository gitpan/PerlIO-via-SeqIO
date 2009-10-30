#-*-perl-*-
#$Id: 002_seqio.t 529 2009-10-29 18:48:14Z maj $

use strict;
use warnings;
use Module::Build;
our $home;
BEGIN {
    use File::Spec;
    our $home = '.';
    push @INC, File::Spec->catfile($home,'lib');
}
use Test::More tests => 32;
use File::Temp qw(tempfile);
use_ok('IO::Seekable');
use_ok('PerlIO::via::SeqIO');
use_ok('Bio::SeqIO');
use_ok('Bio::SearchIO');
use PerlIO::via::SeqIO qw(open O T);
use subs qw(O);

diag("Using exported open...");
# seqio tests...
my $testf = File::Spec->catfile($home, 't', 'test.fas');

ok open( my $seqfh, "<:via(SeqIO)", $testf ), "SeqIO open for reading";
ok my $seq = <$seqfh>, "Readline an object";
isa_ok(O($seq), 'Bio::Seq');
is(O($seq)->id, '183.m01790', "correct seq");
my @rest = <$seqfh>;
is( @rest, 13, "got the rest");
close($seqfh);
ok open( 'DATA', "<:via(SeqIO)" ), "open DATA through SeqIO layer";
my $dtapos = tell(DATA);
ok $seq = <DATA>, "Readline an object";
isa_ok( O($seq), 'Bio::Seq');
is(O($seq)->id, '183.m01790', "correct seq");
{ local $_ = $seq;
  is( O->id, '183.m01790', "correct seq (O with \$_)");
}
@rest=<DATA>;
is( @rest, 13, "got the rest");


seek(DATA, $dtapos, 0);
ok open( 'DATA', "<:via(SeqIO)" ), "open DATA through SeqIO layer";

my ($tmph, $tmpf) = tempfile(DIR=>'.');
$tmph->close;
unlink $tmpf or diag("tempfile unlink issue: $!");

ok open( $tmph, ">:via(SeqIO::embl)", $tmpf) , "open tempfile for Embl write";
while(<DATA>) {
    print $tmph $_;
}
close($tmph);

ok open(my $emblh, "<",  $tmpf), "open embl outfile as text file";
#checksum on id numbers...
my $a = 0; $a += (split /\s+/)[2] for grep (/^SQ/, <$emblh>);
is( $a, 14262 , "fasta converted to embl" );
close($emblh);

seek(DATA, $dtapos, 0);
ok open ($tmph, ">:via(SeqIO::gcg)", $tmpf), "open tempfile for GCG write";
while(<DATA>) { print $tmph $_; }
close($tmph);

ok open(my $gcgh, "<", $tmpf), "open gcg outfile as text file";
# check sum on lhs coordinates
my @a = grep(/^.+$/, <$gcgh>);
$a = 0; $a += (split /\s+/)[1] for grep (/^\s+[0-9]/, @a);
is ($a, 142639, "fasta converted to gcg");
close($gcgh);

# change formats on fly
seek(DATA, $dtapos, 0);
ok open(my $formh, ">:via(SeqIO::embl)", $tmpf), "open tempfile for writing Embl";
$seq = <DATA>;
eval { print $formh $seq };
ok !$@, "print as embl";
ok ($formh->set_write_format('gcg'), "set gcg format");
eval { print $formh $seq };
ok !$@, "print as gcg";
close($formh);
undef $formh;

ok open($formh, "<", $tmpf), "open multi-format dump";
@rest = <$formh>;
@rest = grep /1080/, @rest;
ok $rest[1] =~ /^SQ.*?\s1080/, "embl there first";
ok $rest[3] =~ /^183.*?Length: 1080/, "gcg there second";
close($formh);

# render de novo seq objects
undef $tmph;
ok open($tmph, ">:via(SeqIO::embl)", $tmpf), "open tempfile for writing Embl";
diag("redefine warnings may appear; do not be alarmed...");
$testf = File::Spec->catfile($home,'t', 'test.bls');
my $result = Bio::SearchIO->new( -file=>$testf, -format=>'blast' )->next_result;
my $hsp_ct = 0;
while(my $hit = $result->next_hit()){
    while(my $hsp = $hit->next_hsp()){
	$hsp_ct++;
	my $aln = $hsp->get_aln;
	for ($aln->each_seq) {
	    print $tmph T($_);
	}
	1;

    }
}
close($tmph);
undef $tmph;

open($tmph, "<", $tmpf);
@rest = <$tmph>;
is ( scalar grep(/^SQ/, @rest), 2*$hsp_ct, "wrote converted hsps as embl");
close($tmph);
unlink $tmpf or diag("tempfile unlink issue: $!");
undef $tmph;

# STDIN/STDOUT checks 
SKIP : { 
    my $cat;
    for ($^O) {
	/unix/i and $cat = 'cat';
	/cygwin/i and $cat = 'cat';
	/ms/i and $cat = 'type';
    }
    skip "Don't know your OS; complain to author", 1 unless $cat;
    my $test = File::Spec->catfile($home, 't', 'test.fas');
    my $lib = File::Spec->catfile($home, 'lib');
    my $perl = Module::Build->current->perl;
    my @gcg = `$cat $test | $perl -I$lib "-MPerlIO::via::SeqIO qw(open)" -e "open(STDIN, '<:via(SeqIO::fasta)'); open(STDOUT, '>:via(SeqIO::gcg)'); while (<STDIN>) { print }"`;
    my @a  = grep(/^.+$/, @gcg);
    $a = 0; $a += (split /\s+/)[1] for grep (/^\s+[0-9]/, @gcg);
    is ($a, 142639, "STDIN/STDOUT work in shell");
}
1;
__END__
>183.m01790 |||similar to unknown protein||chr_13|chr13|183
ATGGACGACAAAGAACTCGAAATACCGGTAGAACATTCCACGGCTTTCGGTCAGCTCGTG
ACGGGTCCCCCGGGAGCGGGTAAATCGACCTATTGTCATGGCTTACATCAGTTCCTTACA
GCCATCGGTAGACCAGTGCATATCATCAACCTCGATCCTGCAGTCCCAAACCCTCCGTAT
CCATGCTCTATAAACATCACGGAACTCATCACACTCGAAAGTGTTATGGAGGAATACAAT
CTAGGACCGAATGGGGCGATGCTTTATTGTATAGAATTCTTAGAGGCCAATTTTGACTGG
CTAGTGGAGAGGCTGGATGAGGTCTTGGCTGAAGAGGGGGGGAATGGATATGTGGTGTTT
GATACGCCGGGTCAAGCAGAGTTATGGACGAACCATGATAGTTTGAAGAACGTGGTCGAA
AAGTTGGTCAAGATGGACTATAGACTAGCGGCTGTGCATCTCAGCGACGCGCACTACATA
ACAGATGCCTCAAAATTCATCTCTGTAGTTTTGCTAGCTCTTCGGGCGATGCTGCAAATG
GAAATGCCGCATTTGAATGTGCTCAGCAAAATAGATTTGATATCAACTTATGGAGAGCTC
CCGTTCGACTTGAGCTATTACACAGAAGTCCAAGATCTGTCATACTTACTGGGCAGTCTG
GATTCAGACCCTCGAACAGCAAAGTACCACAAGTTAAATAAAGCGTTGGTAGAGCTTATA
GAAGGCTTTTCATTAGTCGGATTTCAAACCCTCGCTGTTGAGGACAAAGAATCAATGCTT
AATATCGTCCGTCTTGTCGATAAGATGACGGGCTACATATTTATTCCGTCTGGCGACCTC
GAAGGAACCAACGCCATCAATACCCAAGCTCTGTTTGGTAGTGCCATGTCGTCGGCGAAG
CTTACAGGAAGAGCAGGCGGGGACGTAAGAGATGTTCAGGAGAGATGGATGGATAACAAG
GAGGCTTGGGATGAATGGGAGAAGAAAGAATGGAAGAGAGAAGCAGAGATAAGAGCCCAG
ATGGGCACTGGAATACCAGAAGGGATGAAAGGCGGTGAAGATGCGGAAAGTACAGGTATA
>AAL117C location=AgChr1:complement(140329..141372)
ATGGCGTATGGACAGATTGTGATAGGTCCACCGGGGTCTGGGAAGTCGACATACTGTAAT
GGGTGCAGCCAGTTCTTTAATGCCATCGGCAGACACGCTCGGATCGTGAACATGGACCCT
GCAAACGACTCGCTGCCCTACCAATGCGATGTAGACATTCGAGACTTTATTACTCTGGAG
GAAATCATGAACGAGCAGCACCTGGGGCCCAACGGAGGGCTGGTGTATGCGTTTGAGTCG
GTGGAGCACTCACTGTCGCTGTTTGCGCTGCAGATCAAGACGCTGGTCAAGGATGAGAAC
GCATATCTCGTCTTTGACTGCCCCGGTCAGGTGGAGCTGTTCACGCATCACTCGGCGCTC
TCCAAGATATTCCAGCAGCTGGTGCGCGACTTGGACCTACGAGTGTGCGTGGTGAACTTG
ATGGACAGCATCTACATTACATCGCCGTCGCAGTATGTCTCGGTACTGCTGCTGGCGCTG
CGCTCAATGTTGATGATGGACCTGCCCCATATTAACGTTCTCTCTAAGATCGATATGCTG
AGCTCGTACGGCGACCTGCCGTTCCGGCTCGACTACTATACCGAGGTGCAAGACTTGGAG
TATCTGCAACCGCATATTGAACGCGAACACAAGGGAGCCAAGGCGTTGAGGTACCGCCGA
CTAACGGAGGCCATAGGAGAGGTGGTTTCGGACTTCAACCTGGTCGCCTTCGAGGTGCTT
TGCGTCGATGACAAACAGAGCATGATCAACTTGCAAAGCGCAATCGACAAGGCCAATGGT
TATATTTTTGGTGCCTCCGAAGTTGGTGGCGATACTGTGTGGGCGGAGGCAACCCGCCAG
GGCACTGCTGCAATTGAATATGACATTCAGGACAGATGGATCGACAACAAGGACTTTTAT
GACAAGGAGGAAGAGGCTAGGCGCAAGAAGTTACTTGAGGAGCATGAGCTTCTGGAGAAA
GAAGTTGATGTCAACCAGGATGATGAATGGGAACGCGCAGTGAAGGAATGGGAGTCCCAG
CACTCTGTGAACTTCGTTAAA
>AN2438.1 hypothetical protein (53856 - 52862)
ATGAGTGAGGATCAATTGGGTCCGAACGGCGGTGTTTTGTATGCGTTGGAAGAGCTAGAG
GAGAACTTTGACTTCTTGGAGGAAGGGTTGAAAGAGCTCGGAGAGGACTATATTATCTTC
GATTGTCCCGGCCAGGTAGAAATTTTCACTCACCATTCGTCCTTACGGAATATCTTCTTC
AAGATCCAGAAGATGGGCTATAGACTAATAGTACTACACCTAATCGACTCCTACAACCTC
ACCCTGCCATCGATGTACATCTCCTCTCTTATTCTATGCTTGCGTGCCATGCTCCAAATG
GACCTTCCACATCTCAACGTCCTAACAAAAATCGATAATTTGTCCAATTATACTTCGCTG
CCTTTCAACCTAGATTTCTACACCGAGGTTCAGGACCTTACATACCTCCTCCCCCACTTA
GAGGCAGAGTCCTCCCGGCTATCGCACGAGAAGTTCGGAGCACTGAACAACGCCATCATC
ACACTGATTGAGGAGTTTGGACTCGTGGGCTTCGAAACACTGGCTGTAGAAGATAAAAAG
AGCATGATGAATTTGCTCCGGGCCATTGACCGCGCAAGTGGATACGTGTTTGGGCCTGCA
GAAGGCGCAAATGACTCCGTTTGGCAAGTGGCTGTTCGGGAAGGAATGGGGTCCATGGAT
ATCCGTGATATTCAAGAGCGTTGGATAGATGCCAAAGACGAGTACGATGAGTTGGAACGA
CGGCAGCGAGAGGAGGAGATAAAAAATCACCAGCAAGCTGCAACCTACCAGGCAGGGAAC
GAGGACGACGACGATGATAACGATTACGAATTCGGGCGCAGGATGCCTGTACCAGACAGT
GGAGTGAAAGTGATGCGGAAG
>FG05298.1 hypothetical protein (258181 - 259340)
ATGCCTTTCGCGCAACTCGTTCTCGGTAGTCCGGGCTGCGGAAAGAGTACATACTGTGAT
GGCATACAGCTGACCGGTCAAGTGCATCAGTTCCTAGGCGCCATCGGGCGAGCCTGTTCA
GTCGTCAATCTCGATCCTGCCAACGATCATACCAACTACCCTGCAGCTCTCGACATTCGC
AGTTTGATTAAGCTCGAGGAGATTATGAAAGATGATAAATTAGGACCTAATGGCGGCATC
CTGTATGCCCTCGAAGAGTTGGAACACAATTTCGAGTGGTTGGAAGAAGGACTGAAAGAA
TTCAGCGAAGACTATATTCTTTTCGACTGTCCGGGACAAGTGGAACTATATACACACCAC
AACTCCTTGCGAAACATATTCTACAAGCTCCAGAAGATTGGATTCAGGCTTGTTTCCGTC
CACCTCTCCGACTCCTTCTGCCTCACGCAACCGTCGTTATACGTATCGAACGTCCTCCTC
TCCCTTCGTGCGATGATCCAGATGGATATGCCACACATAAATATTCTCTCCAAGATCGAC
AAAGTTGCCGACTACGACGAACTCCCTTTCAACCTCGATTACTACACAGACGTGGACGAC
CTTACATATTTGACACCCCATCTTGAGACAGAGTCGCCCGCTCTGAGGAGTGAGAAATTC
GGCAAGCTCAACGAGGCGATTGCGAATCTGATCGAGAGCTACGGTCTGGTGCGCTATGAA
GTCCTGGCTGTCGAGAACAAGAAAAGCATGATGCATATCCTCCGTGTCATTGACCGTGCT
GGTGGATACGTCTTTGGTAGTGCTGAAGGAGCCAATGATACAGTCTGGTCAGTTGCCATG
AGGAACGAGTCGTCCATGTTGGGGGTGCAGGACATCCAAGAGCGTTGGATCGACCAAAAG
GTGGAATATGATCAAATGGAGCGTGAGGCCGAAGAAGAACAGGCGCGCATCCAAGAAGAA
CAAGCCATGGAGATGGAACAATCACAGCCACCTCCTGCGCCGACAGGTGGCATGGATCCT
GATTTTGGTGACATGACGGTGCCCAAAGATAGTGGGATCAAAGTAGTTAGAAAG
>MG06110.4 hypothetical protein similar to (NCU09745.1) hypothetical protein (25629 - 24026)
ATGGGATTTCTAGGCGCAATAGGGAGAGCATGTTCCGTAGTAAACCTTGACCCGGCCAAT
GACCATACGAGCTATCCATGTGCCCTCGACATACGAAATCTTGTCACGCTGGAGGAAATC
ATGGGAGACGACAATTTGGGGCCAAACGGTGGCATCCTCTACGCTATTGAAGAGCTGGAG
CATAACTTTGAGTGGTTGGAAGATGGTCTGAAAGAGCTTGGGGACGACTACATACTATTC
GACTGCCCGGGCCAGGTCGAGCTGTACACACATCACAATTCATTGCGCAATATCTTCTTC
AAGTTACAAAAGCTCGGCTACAGACTTGTGGTTGTTCACCTCTCGGACAGCATTTGCCTC
ACTCAACCATCGTTGTACATCTCGAATCTCCTCCTCGCTTTGCGCGCCATGCTCCAGATG
GATCTTTCCCATGTCAATGTCCTCACCAAAATCGACAAGGTGTCTTCATATGACAGACTA
GCCTTCAACCTCGACTTTTATACCGAGGTCCACGATCTTTCGTACCTCCTCCCCGAGCTC
GAAGCCGAGAATCCGTCGCTACGCAGCGAAAAGTTCGCCAAGCTAAACCGAGCCGTCGCA
AACTTGATTGAAGACTTTGGGCTCGTCCGGTTCGAAGTCTTGGCTGTCGAGAATAAGAAA
AGTATGATGCATTTGCTCCGGGTCCTCGATCGTGCCAACGGGTACGTTTTTGGTGGGGCC
GAGGGAGCCAACGACACCGTTTGGCAAGTAGCCATGCGCAACGAGGGCTCCCTGATGGGG
GTCCAAGATATCCAGGAGCGCTGGATCGATAACAAAGAGGCTTATGACGAGATGGAGCAG
CGTGAATGGGAGGAACAGGTCAAGGCACAAGAAGCCATGGCCGAAGCCGATGCAGCAGCT
GCTGAAGAGGGCGACGATGACTTGATGGGAGGCCCAGGTGCTCGA
>NCU09745.1 (NCU09745.1) hypothetical protein (81475 - 83184)
ATGACCTCCCCACTGCCAGTGCAGCAGTTTATGGGCGCCATCGGGCGACAATGCTCGGTA
GTCAACCTCGACCCTGCGAACGACCACACCAACTACCCATGCGCGCTCGACATTCGCGAC
CTTGTCACTTTGGAGGAGATTATGGCAGACGACAAATTGGGTCCCAATGGCGGTATTCTG
TACGCACTTGAAGAGCTGGAAAATAACATGGAATGGCTCGAGAACGGCCTCAAGGAGCTT
GGAGAAGACTATGTGCTTTTTGACTGCCCTGGTCAAGTCGAGCTCTACACCCACCACAAC
TCGTTACGCAACATCTTTTACCGGTTACAGAAGCTGGGCTACAGGCTGGTAGTTGTCCAC
CTTTCCGACTGCTTCTGCCTCACACAACCATCGCTCTACATTTCCAACGTCCTCCTCTCT
TTGCGCGCCATGTTGCAAATGGACCTTCCCCACATCAACGTCCTGACCAAGATTGACAAG
ATCTCGTCCTACGATCCTCTTCCATTCAACCTCGACTATTACACCGAAGTACAAGACCTA
CGGTACCTCATGCCGTCCCTCGACGCGGAATCGCCTGCCCTGAAGAAAGGCAAGTTCACC
AAGCTTAACGAGGCCGTTGCGAACATGGTTGAGCAGTTCGGCCTTGTCAGCTTCGAGGTG
CTGGCAGTCGAGAACAAGAAGAGTATGATGCATCTGTTGCGCGTGATTGACCGTGCAAGT
GGGTACGTCTTTGGCGGCGCTGAGGGAACGAACGACACCGTCTGGCAGGTTGCCATGCGC
AACGAGTCATCATTGCCCGATGCTCTTGATATTCAAGAGAGGTGGATCGATAGCAAAGAA
GAGTATGACGAGATGGAGCGGAAGGAGGAGGAAGAACAAGAAAAACTGCGGGCGGAGCAG
GCACGGGCCGCTGAAGAAGCAGGTCTCGGTGACGGCTCGGTCCCTGGAGTGGCGCCACAG
TTCACCAGTGGCTCGGGAATCCGTGTGACGCTTAGCCTAGTGGCCGCTTTTACCAAATAT
AGCGATCTT
>SPAC144.07c SPAC144.07c conserved eukaryotic protein; ATP-binding protein; similar to S. cerevisiae YOR262W
ATGCCATTTTGTCAAGTGGTCGTTGGACCTCCGGGTTCTGGGAAATCAACTTACTGTTTC
GGAATGTACCAATTATTATCTGCCATAGGAAGGAGTAGTATTATCGTCAATCTTGACCCA
GCAAATGACTTTATCAAATACCCATGCGCAATTGATATTCGTAAAGTTCTCGATGTTGAG
ATGATCCAAAAAGACTATGATTTAGGACCAAATGGAGCACTTATTTATGCTATGGAAGCA
ATTGAATATCACGTTGAATGGTTGCTTAAGGAGCTAAAAAAGCATCGAGATTCATATGTG
ATATTTGATTGCCCTGGTCAAGTTGAGTTATTTACAAACCATAATTCCTTACAAAAAATA
ATCAAAACTTTGGAAAAGGAACTGGATTATAGACCTGTGTCCGTACAACTTGTAGATGCA
TATTGCTGCACGAATCCTTCTGCATATGTTAGTGCACTGCTTGTTTGCCTAAAGGGGATG
CTTCAGCTGGACATGCCACATGTAAATATTTTGTCGAAGGCTGATTTGCTTTGTACGTAT
GGAACTTTACCAATGAAACTAGATTTTTTTACCGAAGTACAAGACCTTTCATATTTGGCG
CCTTTGCTTGATAGAGATAAACGTCTTCAGCGCTATAGTGATTTAAACAAAGCTATTTGT
GAACTTGTTGAAGATTTTAATCTTGTTTCTTTTGAAGTTGTTGCAGTAGAAAATAAAGCC
AGTATGTTACGTGTTCTTCGAAAAATCGATCAAGCAGGTGGATATGCATATGGATCTACA
GAAATTGGTGGTGATGCCGTTTGGGTGAATGCCGTTCGTCAAGGTGGAGACCCTCTTCAA
GGTATTTCGCCTCAGGAAAGATGGATTGACAAGAAAGAGGAATATGACAAATATGAATGG
GAATTAGAGCAAAAATCGACCATGGACGAAGATGAAAATGAAGGG
>Sbay_Contig635.43 YOR262W, Contig c635 67551-68594
ATGCCTTTTGCTCAGATTGTTATTGGACCCCCGGGTTCAGGGAAGTCTACGTATTGTAAC
GGATGTTCACAATTTTTTAATGCTATTGGGAGACATTCTCAGGTGGTAAATATGGATCCC
GCCAATGATGCCTTACCTTATCCGTGTGCTGTGGATATCAGAGATTTTATAACTTTGGAA
GAGATCATGAAAGAGCAACACTTGGGCCCTAATGGTGGTTTGATGTATGCCGTTGAATCT
CTAGATAAGTCCATTGATTTATTTATACTACAGATCAAATCACTTGTAGAAGAAGAGAAG
GCATATGTTGTGTTTGACTGCCCGGGACAAGTTGAGCTGTTTACGCATCATTCTTCATTA
TTCAGCATTTTCAAGAAATTAGAAAAAGAACTAGATATGAGATTCTGTGTGGTGAATTTG
ATTGATTGTTTTTACATGACATCTCCTTCACAATATGTCTCGATTTTGCTCCTGGCATTA
AGGTCTATGCTGATGATGGACCTGCCCCATATCAACGTCTTTTCGAAGATAGATAAGTTG
AAATCATATGGAGAATTGCCATTTAGATTAGATTATTATACAGAAGTTCAAGATTTGGAT
TATTTGGAGCCGTATATTGAAAAAGAAGGTTCTGGTGCACTGGGAAAAAGATATAGCAAA
TTGACTGAAACGATTAGTGAGCTGGTTTCTGATTTTAACCTGGTTTCCTTTGAAGTTTTG
GCTGTGGATGACAAAGAAAGTATGATAAATCTCCAGGGTGTTATTGATAAAGCCAATGGT
TACATATTTGGTGCATCTGAAGTGGGCGGCGACACGGTATGGGCCGAGGCCTCGAGAGAA
GGTGCATTGCTAGCAAGCTATGATATTCAAGATAGGTGGATAGATAATAAAGAAAAATAT
GATAAAGAAGAACAAGAGAAACGGGCTGCAATGGTGAAAGAGCAGGAACTGCAAAATAAA
GAGGTTAATGTAGACGAAGAAGACGAGTGGGAAAATGCACTAAACGACTGGGAAGAAAAA
CAAGGCACAGATTTTGTCAGG
>Scas_Contig692.20 YOR262W, Contig c692 40768-41811
ATGCCATTTGCCCAAATTGTTATCGGACCCCCCGGTTCAGGAAAATCAACATACTGTAAC
GGGTGTTCTCAATTTTTCAACGCCATCGGCAGGCATGGCCAAATAGTGAACATGGATCCA
GCTAATGATGCTCTACCATATCCATGTGCAGTAGACATTCGAGATTTTGTGACTCTGGAG
GAGATTATGCAAGAGCAACAACTGGGCCCCAATGGAGGGTTGATGTATGCTGTGGAATCG
TTAGATGAATCCATCGATCTTTTCATACTACAAATAAAATCTCTAGTTCAAGAGGAGAAG
GCATATTTAGTCTTTGATTGTCCTGGACAAGTAGAGTTGTTTACTCATCATTCATCTCTG
TTCAAAATCTTCAAAAAATTGGAAAAGGAACTAGATATGCGATTTTGTGTGGTGAATTTG
ATTGATTCTTTCTATATTACCTCCCCATCACAGTATGTTTCCATTTTGCTGTTGGCTTTG
AGATCTATGTTAATGATGGACCTACCGCAAATCAATGTTTTCTCCAAGATTGATATGCTG
AAATCCTATGGAGAACTACCTTTTAGATTGGATTATTACACAGAAGTGCAAGATTTAGAT
TATTTACAGCCATTTATTGAGAAGGAGAGTTCCAGTGTTTTGGGTAGAAGATATAGCAAG
TTAACAGAAACGATTAGTGAATTGGTTTCCGATTTTAATTTGGTCTCATTTGAAGTCTTA
GCTGTAGATGATAAACAAAGCATGATTAATTTACAAAGTGTAGTAGACAAGGCTAATGGA
TATATATTTGGAGCATCTGAAGTAGGTGGTGATACTGTTTGGGCAGAAGCCACGCGAGAA
GGTGCAATGATGGTAAATTATGATATACAGGACAGATGGATAGATAACAAAGAAAAGTAC
GATGAAGAGGAGAGAAAAAGACAAGAGGAACAAGCCAAAGAGCAGAACATGCAAGAAAAG
GAGGTAGACGTGGATAATGAGGACGAATGGGAAAAGGCATTGAAGGATTGGGAAGAAAAA
CAAGGAACAGGCTATGTAAGG
>Sklu_Contig2277.4 YOR262W, Contig c2277 4093-5136
ATGCCCTTTGGTCAGATTGTTATCGGCCCTCCTGGTTCAGGAAAGTCTACCTATTGTAAT
GGTTGCTCCCAGTTTTTTAATGCTGTCGGTAGACATGCCCAAGTAATCAACATGGATCCA
GCAAATGATTCGTTACCTTACCCATGTGCCGTTGACATTCGAGATTTCATCACCTTAGAG
GAAATTATGACAGAACAGCAGCTGGGGCCTAATGGTGGATTGATGTACGCCCTAGAATCT
TTGGATAAATCAATCGACTTATTTGTTTTGCAGATCAAATCACTAGTTCAGGATGAACAT
GCTTACGTAGTATTTGATTGTCCGGGGCAAGTGGAGCTTTTTACGCACCATTCGTCCTTG
TTCCGCATATTCAAGAAGTTGGAAAGAGAACTAGATATGAGGTTATGCGTGGTTAATTTA
ATCGATTGTTTTTACATCACCTCTCCTTCACAGTATGTCTCTATTCTTTTGCTAGCTTTG
AGGTCGATGCTGATGATGGACTTACCACACATTAATGTCTTTTCTAAAATTGATTTGTTG
AAATCCTACGGTGAGCTGCCATTCCGACTAGATTATTATACCGAAGTTCAAGAGCTAGAT
TACTTGAAGCCACATATTGACAAGGAAGGGAGCAGCGTCCTTGGAAGGAAATATAGTAGG
TTGACAGAAACCATTAGTGAACTGGTTTCTGACTTTAATCTGGTTTCCTTTGAAGTTTTG
TGTGTTGATGATAAGCAGAGCATGATCAATTTGCAAAGTATTGTGGATAAAGCAAATGGT
TACATATTTGGTGTTTCTGAGATCGGTGGAGATACGGTATGGGCAGAGGCAACGCGACAA
GGCAGTGCAATTGCTAATTACGACATTCAAGAGAGATGGATAGATAATAAAGATATGTAC
GACAGAGAGGAACAGGAAAAACGTGAACAGTTGCTCAAAGAAGAAGAGCTACAGAATAAA
GAAGTAGACGTGGATAAAGGTGATGAGTGGGAAAATGCTTTAAAAGAATGGGAAGAAAAG
CAAGGCATGAGTTATGTAAAA
>Skud_Contig1703.7 YOR262W, Contig c1703 9292-10335 reverse complement
ATGCCATTTGCTCAAATTGTTATCGGCCCACCAGGCTCGGGAAAGTCAACGTATTGTAAC
GGGTGTTCGCAGTTCTTCAACGCCATTGGAAGACATTCTCAAGTGGTGAATATGGATCCC
GCTAATGATGCTTTGCCTTATCCGTGTGCTGTAGATATTAGAGATTTTATAACTTTGGAA
GAGGTTATGCAGGAGCAACAGTTGGGTCCTAATGGTGGTTTAATGTATGCCGTTGAATCC
CTAGATAACTCCATTGATCTATTCATATTACAGATCAAGTCACTTGTAGAAGAAGAAAAG
GCCTACCTTGTGTTTGACTGTCCTGGACAAGTTGAGCTATTCACGCACCATTCATCTTTA
TTTAGCATTTTCAAGAAAATGGAGAAAGAATTGGATATGAGATTCTGTGTCGTAAACTTG
ATTGATTGCTTTTATATGACATCTCCTTCTCAGTATGTTTCAATTTTGCTACTGGCATTA
AGGTCCATGCTAATGATGGATTTGCCTCACATAAACGTTTTTTCCAAAATAGATATGTTA
AAATCATATGGGGAATTACCCTTCAGATTGGATTATTATACAGAGGTCCAGGAGCTAGAT
CATTTGGAGCCATATATTGAAAAGGAAGGCTCTAGCGTTCTAGGAAAAAAATATAGTAAG
TTGACTGAAACGATCAAAGAATTAGTCTCCGATTTTAACTTAGTTTCTTTTGAGGTTCTG
TCCGTGGATGACAAAGAAAGTATGATAAATCTCCAGGGTGTTATTGATAAAGCGAATGGC
TACATATTCGGAGCATCCGAAGTTGGAGGTGATACAGTGTGGGCCGAAGCTTCGAGAGAA
GGTGCATTGTTAGAAAACTACGACATACAGGATAGGTGGATAGATAATAAAGAAACGTAT
GATAAAGAAGAACAAGAGAAGCGTGCATCGCTGTTAAAAGAACAAGAACTGCAGAATAAA
ACGGTTGATGTGAAAGAAGAAGATGAATGGGAAAATGCATTAAAGGAGTGGGAAGAAAAG
CAAGATACGGAGTTTGTCAGA
>Smik_Contig1103.1 YOR262W, Contig c1103 447-1490 reverse complement
ATGCCGTTTGCTCAGATTGTTATTGGCCCACCGGGTTCAGGCAAGTCCACTTATTGTAAC
GGCTGCTCACAGTTCTTCAATGCCATTGGGAGACATTCTCAGGTGGTGAACATGGATCCC
GCTAATGATGCTTTGCCTTATCCTTGTGCTGTGGATATCAGAGATTTTATAACGTTGGAA
GAGATTATGCAAGAGCAACAGTTAGGCCCCAATGGTGGTTTAATGTATGCAGTCGAATCC
TTGGATAAGTCTATTGATTTGTTTTTATTACAGATCAAATCGCTTGTAGAAGAAGAAAAA
GCCTATCTTGTATTCGACTGTCCAGGCCAGGTCGAGTTATTTACTCATCACTCATCCTTA
TTCAATATATTTAAGAAAATGGAGAAAGAATTGGACATGAGGTTCTGTGTAATAAACTTG
ATTGACTGTTTTTACATGACGTCACCCTCACAATATGTCTCAATTTTACTGCTTGCACTA
AGATCCATGTTGATGATGGATCTGCCCCACATAAATGTTTTTTCTAAGATAGATATGTTG
AAATCATATGGAGAACTACCATTTAGACTAGATTATTATACAGAGGTACAGGATCTAGAT
TATTTGGAACCGTATATTGAAAAAGAAGGCTCTAGTGTATTAGGAAAGAAATACAATAAG
TTGACCGACGCAATCAAAGAGCTTGTTTCTGATTTTAACTTGGTTTCCTTTGAGGTTTTG
TCCGTGGATGACAAAGAAAGTATGATAAATCTCCAGGGTGTGATTGATAAAGCAAATGGC
TACATATTTGGTGCGTCTGAGGTTGGTGGTGATACAGTGTGGGCAGAGGCTTCTAGGGAA
GGTGCTCTTTTAACAAGTTACGATATTCAAGATAGGTGGATAGATAATAAGGAAAAGTAT
GACAAAGAAGAAGAAGAGAAACGTGTAATCTTGTTAAAAGAGCAAGAGCTGCAAAATAAA
GCAGTTGACGTGAATGAAGACGATGAGTGGGAAAGTGCGCTCAAGGAATGGGAAGAAAAA
CAAGGTATGGATTTTGTTAGA
>Spar_21273 YOR262W, Contig c261 8817-9860
ATGCCCTTTGCTCAAATTGTTATTGGCCCACCGGGTTCAGGAAAATCAACCTATTGCAAC
GGCTGTTCACAGTTTTTCAATGCCATTGGAAGACATTCTCAGGTAGTAAATATGGACCCT
GCTAATGATGCGTTACCTTACCCATGTGCTGTGGATATTCGAGATTTTATAACTTTGGAG
GAGATTATGCAAGAGCAACAGTTAGGCCCCAATGGTGGTTTGGTGTATGCTGTTGAATCC
TTGGATAAGTCCATTGACTTGTTCATATTACAAATCAAGTCGCTTGTAGAAGAAGAAAAG
GCATATCTCGTATTTGACTGTCCCGGACAAGTGGAGTTATTTACTCATCACTCATCTTTA
TTCAGCATTTTTAAGAAAATGGAAAAAGAATTGGACATGAGATTCTGTGTAGTAAATTTG
ATAGACTGTTTTTACATGACTTCTCCTTCACAATACATCTCCATTTTGCTACTCGCATTA
AGGTCTATGTTAATGATGGATCTACCCCACATTAACGTTTTTTCTAAGATAGATATGTTG
AAATCCTACGGGGAATTACCCTTTAGATTAGATTATTATACAGAGGTTCAGGATCTAGAT
TATTTGGAGCCATATATCGAAAAGGAAGGCTCTAGTGTACTGGGAAAGAAATATAGCAAG
TTAACTGAGACAATCAAAGAGCTTGTTTCAGATTTCAATCTGGTTTCATTTGAGGTCCTG
TCTGTGGATGATAAAGAAAGTATGATAAATCTTCAAGGTGTTATAGATAAAGCAAATGGC
TACATATTCGGCGCATCTGAAGTTGGCGGTGATACAGTTTGGGCTGAGGCCTCTAGAGAA
GGTGCATTACTAGCAAATTACGACATTCAGGACAGATGGATAGACAATAAAGAGAAGTAC
GATAAAGAGGAAGAAGAGAAACGCGCGGCGTTGCTAAAAGAACAAGAGTTGCAAAATAAA
GCCGTTGATGTGAATGAAGAGGATGAGTGGGAAAATGCGCTGAAGGAATGGGAAGAAAAA
CAGGGTACGGATTTCGTTAGA
>YOR262W YOR262W SGDID:S0005788, Chr XV from 817289-818332
ATGCCCTTCGCTCAGATTGTTATTGGTCCACCAGGTTCAGGGAAGTCAACCTATTGCAAC
GGCTGCTCACAGTTCTTCAATGCCATCGGAAGACATTCCCAGGTAGTGAATATGGATCCT
GCTAATGATGCCTTACCTTACCCATGCGCTGTGGATATTCGTGATTTTATAACATTAGAG
GAGATCATGCAAGAGCAACAGTTAGGCCCTAATGGAGGTTTGATGTATGCTGTTGAATCA
TTGGATAATTCTATTGATTTGTTCATTTTACAGATCAAGTCACTTGTAGAAGAAGAAAAA
GCATATCTTGTATTCGACTGTCCGGGCCAAGTGGAGCTATTTACTCATCACTCATCTTTG
TTCAACATCTTTAAAAAAATGGAAAAGGAATTGGACATTAGGTTTTGTGTTGTAAATTTG
ATTGACTGTTTTTACATGACATCCCCTTCACAATATATCTCGATTTTGTTACTTGCATTG
AGGTCTATGTTAATGATGGATCTCCCTCACATCAACGTTTTTTCTAAAATAGATATGCTG
AAATCATACGGAGAATTACCCTTTAGATTAGACTATTATACAGAGGTCCAGGATCTGGAT
TATTTGGAGCCATATATTGAAAAGGAAGGCTCTAGTGTACTGGGAAAGAAATATAGCAAG
TTAACTGAAACAATCAAAGAGCTAGTCTCAGATTTCAACTTAGTATCATTTGAGGTTTTG
TCCGTGGATGACAAAGAAAGTATGATAAATCTTCAAGGTGTTATAGATAAAGCAAATGGC
TACATATTCGGCGCATCCGAAGTTGGTGGTGATACCGTGTGGGCTGAGGCTTCGCGAGAA
GGTGCATTAATAGCGAATTACGACATTCAAGACAGGTGGATAGACAATAAAGAGAAGTAT
GATAAAGAAGAAGAAGAAAAACGTACGGCGTTGTTAAAAGAACAAGAATTGCAAAATAAA
GCTGTTGATGTGAATGAAGAAGATGAGTGGGAAAATGCGCTGAAGGAGTGGGAAGAGAAA
CAAGGAATGGATTTTGTTAGG
