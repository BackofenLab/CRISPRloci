#!/usr/bin/perl
#use feature ':5.10';
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use Cwd qw(getcwd abs_path);
use StructureLibrary::Sequence;
use StructureLibrary::Tools;
use StructureLibrary::Structure;
use StructureLibrary::RNAtoolsConfig; 
use vars qw(%CONFIG);
*CONFIG = \%StructureLibrary::RNAtoolsConfig::CONFIG;
#use FindBin;
#use FindBin qw($Bin);
#use lib $FindBin::Bin();
#use FindBin::Real qw(Bin);

=head1 NAME
compare_accuracies.pl 

=head2 SYNOPSIS

Analyses the structure of CRIPSR direct repeats and predicts the most likely
structure, w.r.t. the context sequence.

COMPULSORY

  --DR=s    CRISPR direct repeat sequence in fasta format 
        OR
  --CRT=s   (AND -msfas) 
            Part of the CRT CRISPR prediction tool output that contains
            the repeats and spacers (this is necessary if there are repeats
            that contain mutations). 
            NOTE: candicates have to be given here with --msfas, because they
            can't be calculated easily.
  
	-a=s      pre-crRNA: CRISPR array (known transcript?)	
	NOTE: 		the headers for the direct repeat and the array need to be unique and
            not identical for one CRISPR locus.
	
OPTIONAL

  -msfas=s  AND -str-f=s		
            motif sequence + structure multiple fasta file
            AND structure format [DB,C] 
  -dot=s    Dotplot of CRISPR array, output from folding algorithm
  -o        Output directory
  -rev      Calculate the reverse complimentary direction of input
  -str2D    Generate 2D structures of candidates for -t many structures
  -nout     creates no output files --> only the summary file to STDOUT
  -t        Number of top results 
  -c        The minimum candidate energy
  --help    Present this message		
					
=head3 DESCRIPTION

You have a CRISPR loci and want information about the putative functional structure.
A possible approach to predict the functional structure is to fold the entire
CRISPR array and see which structure has the highest accuracy accross crRNA units.
This is usually NOT the MFE structure of the repeat!

FOLDING
The direct repeat is folded with RNAfold to see the global dotplot of just
the repeat sequence. Then the entire array is folded with a local folding 
approach (unless a previously calculated dotplot is given).
The default local folding of the array is done with RNAplfold with the following 
options:
RNAplfold  -W 200 -L 100 -noLP

CANDIDATE STRUTURES
To calculate the possible substructures, we use RNAsubopt and then filter the
long list of suboptimal structures according to the following rules and keep
them as a list of candidate strutures.
(1) It is the MFE structure
(2) The substructure has a negative free energy
(3) A stem of at least 3bps exist
(4) No bulges greater than 1nt (1-bulges)

If you do not want the candidates to be automatically generated, then it is 
also give a list of your own candidates (see -msfas). 

BASE-PAIR ACCURACIES
The base-pair accuracy (we do not consider unpaired bases as this was unreliable) 
is calculated for each DR in the array and for each structure candidate. A plot
is created to show the accuracy profile accross the array.

The base-pair accuracy is defined here as the average base-pair probability of a given
structure. 

REVERSE COMPLIMENT 
If you are not sure of the orientation, you can set the -rev flag and then the
calculations are done for the reverse compliment. If you give a dotplot, make
sure that this fits with the reverse compliment sequence!

PREDICTED FUNCTIONAL CRISPR STRUCTURE
For each structure candidate, the mean, median, and standard deviation of the
accuracies across the entire array is calculated (only for perfect repeat matches).
These values indicate which is the functional CRISPR structure, under the assumption
that the functional structure occurs most often in the folded array.

EXTRA INFORMATION
In literature you will often see the DR 3' cleavage site to have the following 
properties
(1) It cuts the repeat so the last 8nt (very conserved) 
	remain attached to the following spacer
(2) These 8nt often begin with the sequence (A)TTG and have the typical
 	ending of the CRISPR, i.e. AAAC, or slight mutations of this.


NOTE1: No mismatches in the repeat within the single array are supported (yet).
NOTE2: use absolute paths to your files if using the cluster to calculate your results!



AUTHOR(S): It was writing by Sita Lange and modified by Omer Alkhnbashi. 

=head4 OPTIONS

        -help   This message.
        
        -DR	<CRISPR_DR_FASTA> e.g. repeat_seq.fasta (either --DR OR --CRT COMPULSORY)
        		CRISPR direct repeat sequence in fasta format	
       			If the orientation is unkown it doesn't matter 
       			as both orientations of the repeat will be 
       			analysed (unless -norev is given).
       			NOTE: only one sequence is allowed here use --CRT instead if
       			there are mutations in the repeat.
        --CRT Part of the CRT CRISPR prediction tool output that contains
            the repeats and spacers (this is necessary if there are repeats
            that contain mutations). 
            NOTE: candicates have to be given here with --msfas, because they
            can't be calculated easily.
       
       	-a=s	
       			pre-crRNA: CRISPR array (COMPULSORY)
       			This is the entire CRISPR array as found in the genome.
       			If the exact transript is known, use this. Otherwise
       			include some flanking sequence to the left of the 1st 
       			repeat and to the right of the last repeat.
       			NOTE: only one sequence is allowed here.
        
        -msfas	<MOTIF_STRUCTURE_FASTA> e.g. motif-structure.fasta 
        		This file contains a list of different structural motifs that
        		should be compared for the same reference sequence. The format
        		should be a multiple fasta file, where each structure motif is
        		given a name in the header and then the next line contains the
        		sequence of the motif and the final (3rd) line contains the
        		structure, given by the format -str-f. See -str-f for more
        		information on the format.
        		(All motifs must be in the correct orientation!)
        		
       	-str-f	<STRUCTURE_FORMAT> e.g. DB (COMPULSORY with -msfas) 
       			This gives the format of your structure in -msfas. The options are:
       			
       			(1) DB = DOT-BRACKET = "().-"
       				matching brackets ( ): base i pairs with base j 
       				. : base i is unpaired
       				- : structure of base i is unkown
       				
       			(2) C = CONSTRAINT = "()x.|<>", which are the symbols for RNAfold
       				matching brackets ( ): base i pairs with base j
       			 	| : paired with another base
 					. : no constraint at all
 					x : base must not pair
 					< : base i is paired with a base j<i
 					> : base i is paired with a base j>i
 					
        -dot   <DOTPLOT> e.g. "seq_dp.ps"
        		The dotplot of the CRISPR array in ps format 
        		(e.g. output of RNAplfold).
        		

        -o		<DIR> e.g. CRISPR_structure
        		The directory name in which to save the results, otherwise
        		the current directory will be used
        		
        -rev	Calculate the reverse compliment of the given array and direct 
        		repeat sequences. If you do not know the orientation of your
        		CRISPR, or you suspect it could operate in both orientations,
        		then you can do the calculations for the reverse compliment.
        
		-str2D	Generate 2D structures of the top -t candidates. 
				This will create nice 2D structures of your repeat sequence 
				for as the -t best candidates, according to the -s score.
		
		-nout	Creates no output files. Only gives the summary file of
				the candidate structures to SDTOUT. This may be useful
				when using the script for predicting the most likely
				functional structure to use in follow-up scripts or
				of large CRISPR datasets. In this case you don't
				want to save all the single files.
				The structure with the highest accuracy score (-s) is
				predicted to be the functional one. 
		
		-t	<NUMBER> e.g. 3
				Number of top candidate results that shall be printed
				to STDOUT and for which special files are created, i.e.
				accuracy profile, 2D structures. 
				Default = 3.
				
		-s 	<MEAN|MEDIAN|STD_DEV|ENERGY> e.g "MEDIAN"
				This is the value by which to sort the top candidate structures	
				Default = MEAN	
				
		-c	<FLOAT> e.g. -0.05
				This is the minimum energy a candidate structure may have. Default
				is zero.

        		
=head5 EXAMPLES



=cut


###############################################################################
# PARSE COMMAND LINE OPTIONS
###############################################################################


# command line options
my ($i_help, $i_man, $i_DR, $i_a, $i_motif_str_fasta, $i_str_format, $i_dot,
	$i_out, $i_rev, $i_str2D, $i_nout, $i_topX, $i_sortby, $i_candE, $i_crt);

my $options = GetOptions ("help"	=> \$i_help,
                         "man"  	=> \$i_man,
                         "DR=s"		=> \$i_DR,
                         "a=s"		=> \$i_a,
                         "msfas=s"	=> \$i_motif_str_fasta,
                         "str-f=s"	=> \$i_str_format,
                         "dot=s" 	=> \$i_dot,
                         "o=s"		=> \$i_out,
                         "rev"	=> \$i_rev,
                         "str2D"	=> \$i_str2D,
                         "nout"		=> \$i_nout,
                         "t=s"		=> \$i_topX,
                         "s=s"		=> \$i_sortby,
                         "c=s"		=> \$i_candE,
                         "CRT=s"	=> \$i_crt
							);
                      
                      
# check input
pod2usage(-exitstatus => 1, -verbose => 1) if $i_help;
pod2usage(-exitstatus => 0, -verbose => 2) if $i_man;
($options) or pod2usage(2);

# has the CRIPSR information been given
if($i_DR){
	(-e $i_DR) or pod2usage("INPUT ERROR: $i_DR does not exist!\n");
	$i_DR = abs_path($i_DR);
} elsif ($i_crt){
	(-e $i_crt) or pod2usage("INPUT ERROR: $i_crt does not exist!\n");
	$i_crt = abs_path($i_crt);
} else {
	pod2usage("INPUT ERROR: --DR or --CRT is compulsory!\n")
}

($i_a) or pod2usage("INPUT ERROR: -a is compulsory!\n");
(-e $i_a) or pod2usage("INPUT ERROR: $i_a does not exist!\n");
$i_a = abs_path($i_a);

#if -msfas is given, check whether the structure format has been given
if($i_motif_str_fasta){
	($i_str_format) or 
		pod2usage("INPUT ERROR: -str-f is compulsory if -msfas is given!\n");
}

$i_str_format = "DB" unless ($i_str_format);
if($i_str_format eq "DB"){
	$i_str_format = $CONFIG{STR_FORMATS}->{DOTBRACKET};
} elsif($i_str_format eq "C"){
	$i_str_format = $CONFIG{STR_FORMATS}->{CONSTRAINT};
} elsif($i_str_format eq "()."){
	$i_str_format = $CONFIG{STR_FORMATS}->{DOTBRACKET};
} elsif($i_str_format eq ".()"){
	$i_str_format = $CONFIG{STR_FORMATS}->{DOTBRACKET};
}


#check the output directory if it is given, otherwise set it to the current dir
if($i_out){
	unless(-e $i_out){
		system("mkdir -R $i_out");
	} else {
		print STDERR "WARNING: writing results to an existing directory: $i_out\n".
					 "This may write over existing files!\n";
	}
} else {
	$i_out = getcwd;
	$i_out .= "/CRISPR_Structure";
	(-e $i_out) or system("mkdir $i_out");
}
## get absolute path of output directory
$i_out = abs_path($i_out);

print STDOUT "Results will be saved to $i_out!\n";
	
unless($i_topX){
	$i_topX = 3;
}

if($i_sortby){
	$i_sortby = uc($i_sortby);
	unless ($i_sortby eq "MEAN" || $i_sortby eq "MEDIAN" || $i_sortby eq "STD_DEV" ||
		$i_sortby eq "ENERGY"){
			pod2usage( "INPUT ERROR: wrong value for -s. ".
				"Choose from [MEAN,MEDIAN,STD_DEV,ENERGY]\n");
	}
}else {
	$i_sortby = "MEAN";
}

$i_candE = 0 unless ($i_candE);
		
###############################################################################
# GLOBAL VARIABLES
###############################################################################
my $CURRUSER = getlogin;
my $CURRDIR = getcwd;
print "-- current directory\t  $CURRDIR\n";
my $TMPDIR = "${CURRDIR}/TMP/";
print "-- TMP ->   $TMPDIR\n";
(-e $TMPDIR) or system("mkdir $TMPDIR");
my $SCRIPTNAME = "CRISPR_structure_accuracies";

my $RNAPLFOLD = $CONFIG{RNAPLFOLD}." ".$CONFIG{RNAPLFOLDPARAMS};
my $RNASUBOPT = $CONFIG{RNASUBOPT}." ".$CONFIG{RNASUBOPTPARAMS}." -e 10";
my $RNAFOLD = $CONFIG{RNAFOLD}." ".$CONFIG{RNAFOLDPARAMS};
my $DOTVIEW = $CONFIG{DOTVIEW}." -strict -pdf ";

#%sequences has the keys ARRAY and DR
# each of the four elements have the keys SEQUENCE, FILE and NAME!
my %SEQUENCES;

my %CANDIDATES;

# candidate names, sorted by the -sortby value
my @SORTED_CANDIDATES;

	
#only one direct repeat sequence has been given
if($i_DR){	
		
set_sequences($i_a, $i_DR, $i_rev);
#sequences_dir($i_out);

fold_DRs_with_RNAfold($i_out, $i_rev) unless ($i_nout);

$i_dot = fold_arrays($i_out, $i_nout, $i_dot) 
	unless ($i_dot);
	
if($i_motif_str_fasta){
	read_given_candidates($i_motif_str_fasta, $i_str_format);

## create candidates automatically with RNAsubopt
} else {
	calculate_candidates($i_candE, $i_rev);
	
	## create candidates directory
	candidates_dir() unless($i_nout);
	
}


my $bp_prob_href = calculate_accuracies($i_dot, $i_str_format);

my ($dr_acc_file, $all_candidates_file, $topX_candidates_aref, $acc_folder) = 
	accuracies_output($i_topX, $i_sortby);
	
## if $i_motif_str_fasta is given, then the "top" candidates are all candidates
## in the given file
if($i_motif_str_fasta){
	my @candidates = sort {$CANDIDATES{$b}->{$i_sortby} <=> 
		$CANDIDATES{$a}->{$i_sortby}} (keys (%CANDIDATES));
	$topX_candidates_aref = \@candidates;
}
	
plot_accuracies($topX_candidates_aref, $dr_acc_file, $i_sortby, $i_motif_str_fasta) 
	unless ($i_nout);


my $ave_dotplot = averaged_DR_dotplot_output($bp_prob_href, $i_sortby, $i_str_format) 
	unless ($i_nout);

# the CRT output was given.
} else{
	
	
}

system("rm -R $TMPDIR");


sub set_sequences{
	my ($array_fasta, $dr_fasta, $isrev) = @_;
	my $FUNCTION = $SCRIPTNAME."::set_sequences";
	
	$SEQUENCES{ARRAY}->{FILE} = $array_fasta;
	$SEQUENCES{DR}->{FILE} = $dr_fasta;
	
	## read given array sequence
	my @array_contents = StructureLibrary::Sequence::read_fasta_file($array_fasta);
	my $num = @{$array_contents[1]};
	die "INPUT ERROR in $FUNCTION: Give only one array sequence!\n" 
		if ($num != 1);

		
	## read given direct repeat sequence
	my @dr_contents = StructureLibrary::Sequence::read_fasta_file($dr_fasta);
	$num = @{$dr_contents[1]};
	die "INPUT ERROR in $FUNCTION: give only one DR sequence!\n" if ($num != 1);

	## process and set sequences and sequence names
	$SEQUENCES{ARRAY}->{NAME} = $array_contents[1]->[0];
	$SEQUENCES{DR}->{NAME} = $dr_contents[1]->[0];
	my $array = $array_contents[0]->{$SEQUENCES{ARRAY}->{NAME}};
	my $dr = $dr_contents[0]->{$SEQUENCES{DR}->{NAME}};
	$array = clean_RNA_sequence($array);
	$dr = clean_RNA_sequence($dr);
			
	if($isrev){
		$SEQUENCES{ARRAY}->{NAME} = $array_contents[1]->[0]."_rev";
		$SEQUENCES{ARRAY}->{SEQUENCE} = 
			StructureLibrary::Sequence::reverse_compliment($array, "RNA");
		$SEQUENCES{DR}->{NAME} = $dr_contents[1]->[0]."_rev";
		$SEQUENCES{DR}->{SEQUENCE} = 
			StructureLibrary::Sequence::reverse_compliment($dr, "RNA");
		
	} else {
		$SEQUENCES{ARRAY}->{SEQUENCE} = $array;
		$SEQUENCES{DR}->{SEQUENCE} = $dr;
	}
}

## folds direct repeat sequences with RNAfold and creates pdfs
sub fold_DRs_with_RNAfold{
	my ($dir, $isrev) = @_;
	my $FUNCTION = $SCRIPTNAME."::fold_DRs_with_RNAfold";
	
	## output folder
	my $res_folder = $dir."/FOLDED_DRs/";
	(-e $res_folder) or system("mkdir $res_folder");
	
	## tmp folder for RNAfold results
	my $tmp_dir = $res_folder."tmp";
	(-e $tmp_dir) or system("mkdir $tmp_dir");
	
	chdir($tmp_dir);
	## go to that directory and fold DRs
	print "running $RNAFOLD < $SEQUENCES{DR}->{FILE}...\n";
	print "*------$RNAFOLD-----------\n";
	system($RNAFOLD." < ".$SEQUENCES{DR}->{FILE});
	
	# Open the directory.
    opendir (DIR, $tmp_dir)
        or die "ERROR in $FUNCTION: Unable to open $tmp_dir: $!";

    # Read in the files.
    my @files = readdir (DIR);
	
	## generate pdf files 
	foreach (@files){
		system("ps2pdf $_ ../$SEQUENCES{DR}->{NAME}_dp.pdf") if($_ =~ /_dp.ps/);
	}
    
    # Close the directory.
    closedir (DIR);	

	chdir($CURRDIR);
#	rmdir($tmp_dir);
	system("rm -R $tmp_dir");
}

############################################################################
## fold the arrays, if no dot is given. save dotplots in the folder and
## highlight the crispr repeats.
sub fold_arrays{
	my ($res_folder, $nout, $dot) = @_;
	my $FUNCTION = $SCRIPTNAME."::fold_arrays";
	
	chdir($TMPDIR);
	system("rm $TMPDIR/*.*"); # first remove any file that may be in here
	
	## fold forward array with RNAplfold
	unless($dot){
		## fold with default folding
		print "running $RNAPLFOLD < $SEQUENCES{ARRAY}->{FILE}...\n";
		system("$RNAPLFOLD < $SEQUENCES{ARRAY}->{FILE}");
	}
	
	
	unless($dot){
		# Open the tmp directory.
		opendir (DIR, $TMPDIR)
	    or die "ERROR in $FUNCTION: Unable to open $TMPDIR: $!";
	    # Read in the files.
	    my @files = readdir (DIR);
		## generate pdf files 
		foreach my $f (@files){
			
			if($f =~ /(.+)_dp\.ps/){
				my $pref = $1;
				if($SEQUENCES{ARRAY}->{NAME} =~ /$pref/){
					$dot = $TMPDIR.$f;
				} else {
					print "file reject: $f\n";
					print STDERR "WARNING in $FUNCTION: this is recognised as not being a".
						" dotplot file: $f!\n"
				}
			} 
		}
	    # Close the directory.
	    closedir (DIR);	
	}
    
    ## check that we now have dot
    die "ERROR in $FUNCTION: There is no dotplot file for the array!" unless (-e $dot);
	
	## create folding directory
	unless($nout){
		## create folder
		my $res_folder = $i_out."/FOLDED_pre-crRNA/";
		(-e $res_folder) or system("mkdir $res_folder");
		chdir($res_folder);
		
		## copy dotplots into this folder and remove the old ones
		system("cp $dot $SEQUENCES{ARRAY}->{NAME}_dp.ps ; rm $dot");
		$dot = $res_folder."/$SEQUENCES{ARRAY}->{NAME}_dp.ps";
		
		## highlight the direct repeat sequence
		print "running ".$DOTVIEW." -dot $dot -seq \"$SEQUENCES{DR}->{SEQUENCE}\"".
			" -o $res_folder/$SEQUENCES{ARRAY}->{NAME}_DRmarked.pdf...\n";
		#system("pwd");
		chdir($CURRDIR);
		#system("pwd");
		system($DOTVIEW." -dot $dot -seq \"$SEQUENCES{DR}->{SEQUENCE}\"".
			" -o $res_folder/$SEQUENCES{ARRAY}->{NAME}_DRmarked.pdf");
	}
	
	chdir($CURRDIR);
	return $dot;
}	 

## this calculates the candidates for a direct repeat in the given direction
## if both directions are needed, call this sub twice
## dr_file is the fasta file of the DR repeat sequence
sub calculate_candidates{
	my ($candidate_cutoff_energy, $isrev) = @_;
	my $FUNCTION = $SCRIPTNAME."::calculate_candidates";
	
	print "running ".$RNASUBOPT." < $SEQUENCES{DR}->{FILE}...\n";
	my @RNAsuboptres = readpipe($RNASUBOPT." < $SEQUENCES{DR}->{FILE}");
	
	my $mfe_energy = 0;
	my $mfe = 0;
	my $subopt = 1;
	
	foreach my $line (@RNAsuboptres){
		
		chomp($line);
		my $n = $line =~ s/\(/\(/g;
		
		if($n >= 3 && $n <= 10 && $line =~ /^([\.\(\)]+)\s+(-?\d+\.\d+)\s*$/){
			my $structure = $1;
			my $energy = $2;
			if($energy < $candidate_cutoff_energy){
				## a suboptimal structure
				if($energy > $mfe_energy){
					## filter out suboptimal structures
					if($structure =~ /^\.*\(*\.?\(+\(\(\.+\)\)\)+\.?\)*\.*$/){
						if($isrev){
							$CANDIDATES{"SUBOPT".$subopt."_R"}->{STRUCTURE} = $structure;
							$CANDIDATES{"SUBOPT".$subopt."_R"}->{ENERGY} = $energy;
						} else {
							$CANDIDATES{"SUBOPT$subopt"}->{STRUCTURE} = $structure;
							$CANDIDATES{"SUBOPT$subopt"}->{ENERGY} = $energy;							
						}
						++$subopt;	
					} 
					
				## mfe structure with a structure energy below 0
				} else {
					if($mfe){
						++$mfe;
						if($isrev){
							$CANDIDATES{"MFE".$mfe."_R"}->{STRUCTURE} = $structure;
							$CANDIDATES{"MFE".$mfe."_R"}->{ENERGY} = $energy;							
						} else {
							$CANDIDATES{"MFE$mfe"}->{STRUCTURE} = $structure;
							$CANDIDATES{"MFE$mfe"}->{ENERGY} = $energy;
						}
					} else {
						$mfe_energy = $energy;
						$mfe = 1;
						if($isrev){
							$CANDIDATES{MFE_rev}->{STRUCTURE} = $1;
							$CANDIDATES{MFE_rev}->{ENERGY} = $2;	
						} else {
							$CANDIDATES{MFE}->{STRUCTURE} = $1;
							$CANDIDATES{MFE}->{ENERGY} = $2;
						}						
					}
				}
			}
		}
	}
	chdir($CURRDIR);
}

## creates candidates diretory
sub candidates_dir{
	my $FUNCTION = $SCRIPTNAME."::candidates_dir";
	
	die "ERROR in $FUNCTION: There are no candidates" unless (%CANDIDATES);
	
	
	## make folder for candidates
	my $res_folder = $i_out."/CANDIDATES/";
	(-e $res_folder) or system("mkdir $res_folder");
	
	my $dr_candies_file = $res_folder.$SEQUENCES{DR}->{NAME}."_candidates.fasta";
	my %contents = ();
	my @sorted_candies = sort {candidate_number($a) <=> candidate_number($b)} 
		keys (%CANDIDATES);
	foreach (@sorted_candies){
		if(exists($CANDIDATES{$_}->{ENERGY})){
			$contents{"$_\t$CANDIDATES{$_}->{ENERGY}"} = 
				"$SEQUENCES{DR}->{SEQUENCE}\n$CANDIDATES{$_}->{STRUCTURE}";
		} else {
			$contents{"$_\t$CANDIDATES{$_}"} = 
				"$SEQUENCES{DR}->{SEQUENCE}\n$CANDIDATES{$_}->{STRUCTURE}";
		}

	}
	StructureLibrary::Sequence::write_hash_to_fasta($dr_candies_file, \%contents);
}



sub read_given_candidates{
	my ($ms_fasta_file, $str_format) = @_;
	my $FUNCTION = $SCRIPTNAME."::read_given_candidates";
	
	## read candidates
	my ($ms_contents_href, $ms_order_aref) = 
		read_seq_structure_fasta_file($ms_fasta_file,$str_format);
	
	die "INPUT ERROR in $FUNCTION: There are no candidates in $ms_fasta_file\n" 
		unless (@{$ms_order_aref} > 0);
	
	foreach my $c (@{$ms_order_aref}){
		## check that the sequence is the same as the dr sequence
		$ms_contents_href->{$c}->{SEQUENCE} = 
			clean_RNA_sequence($ms_contents_href->{$c}->{SEQUENCE});
		
		## check that the structure fits to the sequence
		StructureLibrary::Structure::check_bps_wrt_seq(
			$ms_contents_href->{$c}->{SEQUENCE}, 
			StructureLibrary::Structure::structure2basepairs(
			$ms_contents_href->{$c}->{STRUCTURE}, 1, $str_format, 1));
			
		## set candidate
		$CANDIDATES{$c}->{STRUCTURE} = $ms_contents_href->{$c}->{STRUCTURE};
		
		
	}
}


sub calculate_accuracies{
	my ($dot, $str_format) = @_;
	my $FUNCTION = $SCRIPTNAME."::calculate_accuracies";
	
	my $array_seq = $SEQUENCES{ARRAY}->{SEQUENCE};
	my $dr_seq = $SEQUENCES{DR}->{SEQUENCE};
	
	die "INPUT ERROR in $FUNCTION: The dotplot is not a valid file: $dot\n" 
		unless(-e $dot);
	
	## contents for output file
	my @out_contents = ();
	
	# $MOTIF_STRUCTURE_HREF = $candidates_href;
	
	## check number of candidates
	my @candidate_names = keys (%CANDIDATES);
	my $num_candidates = @candidate_names;
	unless ($num_candidates > 0){
		die "INPUT ERROR in $FUNCTION: There are no candidate struture motifs!\n";
	}
	
	# get base-pair probabilities
	unless(-e $dot){
		die "INPUT ERROR in $FUNCTION: The given dotplot does not exist: $dot\n";
	}
	my ($bp_prob_href,$tmp_seq) = parse_dotplot_return_bp_probs_and_sequence($dot);
	$tmp_seq = clean_RNA_sequence($tmp_seq);
	unless ($tmp_seq eq $array_seq){
		die "INPUT ERROR in $FUNCTION: The given dotplot $dot does not ".
			"fit to the array sequence!\n";
	}
	
	## calculate lengths
	my $dr_length = length($dr_seq);
	my $array_length = length($array_seq);
	
	## get repeat positions (only perfect matches will be found!)
	my @positions = StructureLibrary::Sequence::find_motif($array_seq, $dr_seq);
	(@positions) or die "ERROR in $FUNCTION: the repeat sequence could not be ".
		"found in the array sequence!\n"; 
	
	# process each motif one-by-one
	foreach my $mID (@candidate_names){
		
		
		my @accuracies = ();
		my $sum_acc = 0;
		
		# each occurence is calculated separately
		foreach my $occ (@positions){
			
#			
			
#			
			
			## calculate base-pair accuracy score
			my @acc_all = accuracy($CANDIDATES{$mID}->{STRUCTURE}, 
						$str_format, ($occ+1), $bp_prob_href,0,0,0,
						$array_seq, $array_length,"normalise",1);
			
			## only use the base-pair accuracy!!
			push(@accuracies, StructureLibrary::Tools::round_Xdp($acc_all[1], 4));
			$sum_acc += $acc_all[1];
		}
		## save the bp accuracies and the occurrence positions for each str motif
		$CANDIDATES{$mID}->{BPACCURACIES} = \@accuracies;
		$CANDIDATES{$mID}->{OCCURRENCES} = \@positions;	
	}
	return $bp_prob_href;
}


## prints the accuracies to STDOUT and creates accuracy files
sub accuracies_output{
	my ($topX, $sortby, $nout) = @_;
	my $FUNCTION = $SCRIPTNAME."::accuracies_output";
	
	## make folder for accuracies
	my $res_folder = $i_out."/ACCURACIES/";
	unless ($nout){
		(-e $res_folder) or system("mkdir $res_folder");
	}
	
	my @dr_acc_contents = ();
#	my %all_candidate_summary = ();
	my @candidates = sort {candidate_number($a) <=> 
			candidate_number($b)} keys (%CANDIDATES);
#	my @sorted_candidates = ();
	my @candidate_summary_contents = ();
	
	## calculate mean, median, std_dev 
	push(@dr_acc_contents, "Candidate\tDRnumber\tDRposition\tBPaccuracy\tSequence_Structure");
	my $dr_acc_file = $res_folder."ARRAY_".$SEQUENCES{ARRAY}->{NAME}."_".
		$SEQUENCES{DR}->{NAME}."_BPaccuracies.tab";
	foreach my $c (@candidates){
		
		## get the energy, if it exists
		unless(exists($CANDIDATES{$c}->{ENERGY})){
			$CANDIDATES{$c}->{ENERGY} = "NA";
		}
		
		## calculate mean, median, and standard deviation
		my ($std_dev, $mean) = StructureLibrary::Tools::standard_deviation(
			@{$CANDIDATES{$c}->{BPACCURACIES}});
		my $median = StructureLibrary::Tools::median(
			@{$CANDIDATES{$c}->{BPACCURACIES}});
		$CANDIDATES{$c}->{STD_DEV} = StructureLibrary::Tools::round_Xdp($std_dev, 4);
		$CANDIDATES{$c}->{MEAN} = StructureLibrary::Tools::round_Xdp($mean, 4);
		$CANDIDATES{$c}->{MEDIAN} = StructureLibrary::Tools::round_Xdp($median, 4);
	}
	
	## sort candidates according to the sortby value
	@SORTED_CANDIDATES = sort {$CANDIDATES{$b}->{$sortby} <=> 
		$CANDIDATES{$a}->{$sortby}} (@candidates);
	
	## create the accuracies file for each DR position and candidate
	for (my $x = 1 ; $x <= @SORTED_CANDIDATES ; $x++){
		my $c = $SORTED_CANDIDATES[$x-1];
	
		## go through all positions
		my $num_pos = @{$CANDIDATES{$c}->{OCCURRENCES}};
		for (my $i = 0 ; $i < $num_pos ; $i++){
			push(@dr_acc_contents, "$c\t".($i+1)."\t".
				($CANDIDATES{$c}->{OCCURRENCES}->[$i]+1)."\t".
				$CANDIDATES{$c}->{BPACCURACIES}->[$i]."\t".$x.
				get_structure_sequence_notation($CANDIDATES{$c}->{STRUCTURE}, 
				$SEQUENCES{DR}->{SEQUENCE}));
		}
	}
	
	## write file for accuracies
	StructureLibrary::Tools::write_file_from_array($dr_acc_file, \@dr_acc_contents)
		unless ($nout);
	
	## write summary file
	my $all_candidates_file = $res_folder."ARRAY_".$SEQUENCES{ARRAY}->{NAME}."_".
		$SEQUENCES{DR}->{NAME}."_candidates_summary.tab";

	push(@candidate_summary_contents, "Rank\tCandidate\tSequence_Structure\tEnergy\t".
		"Mean_BPaccuracy\tMedian_BPaccuracy\tStandard_deviation_BPaccuracy");
	print STDOUT "\nRank\tCandidate\tSequence_Structure\tStem\tEnergy\t".
		"Mean_BPaccuracy\tMedian_BPaccuracy\tStandard_deviation_BPaccuracy\n";
#	print STDOUT "\nTop candidate structures\n";
	for (my $i=0 ; $i < @SORTED_CANDIDATES ; $i++){
		my $c = $SORTED_CANDIDATES[$i];
		my $line = ($i+1)."\t$c\t".
			get_structure_sequence_notation($CANDIDATES{$c}->{STRUCTURE}, 
			$SEQUENCES{DR}->{SEQUENCE})."\t";
		$line .= $CANDIDATES{$c}->{ENERGY}."\t";
		$line .= $CANDIDATES{$c}->{MEAN}."\t";
		$line .= $CANDIDATES{$c}->{MEDIAN}."\t";
		$line .= $CANDIDATES{$c}->{STD_DEV};
		push(@candidate_summary_contents, $line);
		
		if($i  < $topX){
			my $t = $i+1;
			print STDOUT "$line\n";
		}
	}
		
	StructureLibrary::Tools::write_file_from_array($all_candidates_file, 
		\@candidate_summary_contents) unless ($nout);
	
	## get top candidate names
	my $n = @SORTED_CANDIDATES;
	$topX = $n if($n < $topX);
	my @top_candidates = @SORTED_CANDIDATES[0..($topX-1)];
	

	print STDOUT "\n";
	return ($dr_acc_file, $all_candidates_file, \@top_candidates, $res_folder);
}


sub get_structure_sequence_notation{
	my ($structure, $seq) = @_;
	my $FUNCTION = $SCRIPTNAME."::get_structure_sequence_notation";
	
	my $len = length($structure);
	unless ($len eq length($seq)){
		die "INPUT ERROR in $FUNCTION: ".
			"The structure and sequence are of unequal lengths!\n";
	}
	
	my $descr = "";
	my $ispaired = 0;
	my $isunpaired = 0;
	my @str = split("", $structure);
	for (my $i = 0 ; $i < $len ; $i++){
		my $s = substr($structure, $i, 1);
		my $c = substr($seq, $i, 1);
		if($s eq "(" || $s eq ")" ){
			$descr .= "-" unless ($ispaired);
			$descr .= $c;
			$isunpaired = 0;
			$ispaired = 1;
		} else {
			$descr .= "-" unless ($isunpaired);
			$descr .= lc($c);
			$isunpaired = 1;
			$ispaired = 0;
		}
	}
	$descr .= "-";
	return $descr;
}

## create bar graph of accuracy profile
sub plot_accuracies{
	my ($topX_candidates_aref, $acc_tab, $sortby, $ms_fas) = @_;
	my $FUNCTION = $SCRIPTNAME."::plot_accuracies";
	
	## make folder for plots
	my $res_folder = $i_out."/RPLOTS/";
	(-e $res_folder) or system("mkdir $res_folder");
	
	my $topX = @{$topX_candidates_aref};
	
	
	my $topXfile = "";
	if($ms_fas){
		$topXfile = $res_folder.$SEQUENCES{ARRAY}->{NAME}."_".$SEQUENCES{DR}->{NAME}.
		"_given_candidates.tab";
	} else {
		$topXfile = $res_folder.$SEQUENCES{ARRAY}->{NAME}."_".$SEQUENCES{DR}->{NAME}.
		"_top$topX.tab";
	}
	
	system("grep Candidate $acc_tab > $topXfile");
	for (my $i=0 ; $i < $topX ; $i++){
		system("awk \'\$1 == \"$topX_candidates_aref->[$i]\"\' $acc_tab >> $topXfile");
	}
	
	## colours
	my @colours = ("red","orange","yellow", "darkgreen", "darkblue", "darkgrey",
		"cyan", "blue","purple", "magenta","lightgrey","lightgreen");
	die "ERROR in $FUNCTION: not enough colours for $topX top hits!\n" 
		if ($topX > @colours);

	my $file_prefix = "";
	if($ms_fas){
		$file_prefix = $res_folder.$SEQUENCES{ARRAY}->{NAME}."_".
		$SEQUENCES{DR}->{NAME}."_given_candidates.";
	} else {
		$file_prefix = $res_folder.$SEQUENCES{ARRAY}->{NAME}."_".
		$SEQUENCES{DR}->{NAME}."_top$topX.";
	}
		
	open(RCMD, ">".$file_prefix."R") || 
		die "couldn't open the file for writing: ".$file_prefix."R \n";
		
	## load library
	print RCMD "library(ggplot2)\n";
	
	## read accuracies file
	print RCMD "topXacc <- read.table(\"$topXfile\", header=T)\n";
	if($ms_fas){
		print RCMD "title_acc <- \"Accuracies of selected DR structures ".
			"in $SEQUENCES{ARRAY}->{NAME}\"\n";
	} else {
		print RCMD "title_acc <- \"Accuracies of top $topX DR structures ".
			"in $SEQUENCES{ARRAY}->{NAME}\"\n";
	}
	
	print RCMD "plot <- ggplot(data=topXacc) + geom_bar(aes(DRnumber, BPaccuracy),".
		"stat=\"identity\", position=\"dodge\")";

	print RCMD " + aes(fill=Sequence_Structure) ";
	print RCMD "+ theme_bw() + labs(title=title_acc) + ";
	print RCMD "scale_fill_manual(values=c(";
	for (my $i=0 ; $i < $topX - 1 ; $i++){
		print RCMD "\"$colours[$i]\",";
	}
	print RCMD "\"".$colours[($topX-1)]."\"))";
	
	print RCMD "+ scale_x_continuous(name=\"DR number in array\") + ";
	
	print RCMD "scale_y_continuous(\"Base-Pair Accuracy\") + ylim(0,1)\n";
	print RCMD "ggsave(plot, filename=\"${file_prefix}pdf\", w=15, h=5)\n";
	print RCMD "ggsave(plot, filename=\"${file_prefix}png\", w=15, h=5)\n";

	
	close(RCMD);
	
	print "running R for plotting...\n";
	system("R --vanilla --slave -f ".$file_prefix."R");
	
}

sub create_2D_structures{
	
}

sub averaged_DR_dotplot_output{
	my ($bp_prob_href, $sortby, $str_format) = @_;
	my $FUNCTION = $SCRIPTNAME."::averaged_DR_dotplot_output";
	
	my @candidates = keys (%CANDIDATES);
	
	die "ERROR in $FUNCTION: There are no candidates!" unless (@candidates);
	

	unless(@SORTED_CANDIDATES){
		@SORTED_CANDIDATES = sort {$CANDIDATES{$b}->{$sortby} <=> 
			$CANDIDATES{$a}->{$sortby}} (@candidates);
	}
	
	## check that the dir for the folded arrays exists, or create it
	my $res_folder = $i_out."/FOLDED_pre-crRNA/";
	(-e $res_folder) or system("mkdir $res_folder");
	my $dp_file = $res_folder.$SEQUENCES{ARRAY}->{NAME}."_".$SEQUENCES{DR}->{NAME}.
		"_average_dp.ps";
	
	$CANDIDATES{MFE}->{STRUCTURE} = $CANDIDATES{$SORTED_CANDIDATES[0]}->{STRUCTURE}
		unless (defined $CANDIDATES{MFE}->{STRUCTURE});
	print "test: ".$CANDIDATES{MFE}->{STRUCTURE}."\n";
	
	
	## create average dotplot
	StructureLibrary::Structure::average_crispr_dotplot(
		$bp_prob_href, $SEQUENCES{ARRAY}->{SEQUENCE}, 
		$SEQUENCES{DR}->{SEQUENCE}, $CANDIDATES{$SORTED_CANDIDATES[0]}->{STRUCTURE}, 
		$CANDIDATES{MFE}->{STRUCTURE},$str_format, $dp_file);
		
	return $dp_file;
}


sub candidate_number{
	my $candidate = shift;
	my $FUNCTION = $SCRIPTNAME."::candidate_number";
	
	die "ERROR in $FUNCTION: No candidate was given!\n" unless ($candidate);
	my $mfe_num = -5;
	if($candidate =~ /MFE(\d?)/){
		if($1){
			return ($mfe_num + $1);
		} else {
			return $mfe_num;
		}
	} elsif($candidate =~ /SUBOPT(\d+)/){
		return $1;
	} elsif ($candidate =~ /(\d+)/){
		return $1;
	} else {
		return 0;
	}
}


##################################################################################
sub clean_RNA_sequence{
	my $rna_seq = shift;

	$rna_seq =~ s/[\s\n]//g; #remove newlines and spaces	
	$rna_seq = uc($rna_seq); ## upper case
	$rna_seq =~ tr/T/U/; ## RNA characters
	$rna_seq =~ s/[^ACGU]/N/g; ## change unkown characters into Ns
	
	return $rna_seq;
}
