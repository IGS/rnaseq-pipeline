#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

cuffmerge.pl - script to merge Cufflinks transcript GTF files with or without reference.

=head1 SYNOPSIS

    cuffmerge.pl --i cufflinks_gtf] [--a annotation_file] [--p outprefix]
    		       [-r seq_dir/seq_fasta] [--o outdir] [--cb cufflinks_bin_dir] [--v] 

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --i <cufflinks_gtf>            = /path/to/cufflinks gtf file for a single sample or 
                                     a listfile of cufflinks gtf files for multiple samples.
                                     The gtf filename should be in the format /path/to/<sample_name>.*.transcripts.gtf

    --a <annotation_file>          = /path/to/annotation file in GTF format. Optional

    --p <prefix>                   = Output files prefix. Optional. [cuffmrg]

    --o <output dir>               = /path/to/output directory. Optional. [present working directory]
	
	--r <reference_seq>            = /path/to/seq_dir or /path/to/seq_fasta. Optional

    --cb <cufflinks_bin_dir>       = /path/to/cufflinks binary directory. Optional. [/usr/local/bin]

    --args <other_params>          = additional Cuffmerge parameters. Optional. Refer Cufflinks manual.

    --v                            = generate runtime messages. Optional

=head1 DESCRIPTION

The script executes the cuffmerge script from the Cufflinks RNA-Seq Analysis package.

=head1 AUTHOR

 Amol Carl Shetty
 Bioinformatics Software Engineer II
 Institute of Genome Sciences
 University of Maryland
 Baltimore, Maryland 21201

=cut

################################################################################

use strict;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Spec;

##############################################################################
### Constants
##############################################################################

use constant FALSE => 0;
use constant TRUE  => 1;

use constant BIN_DIR => '/usr/local/bin';

use constant VERSION => '1.0.0';
use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM." version ".VERSION."\n";

GetOptions( \%hCmdLineOption,
            'gtffile|i=s', 'annotation|a=s', 'prefix|p=s', 'outdir|o=s', 
            'cufflinks_bin_dir|cb=s', 'ref_sequence|r=s', 
			'args=s', 
            'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

my (@aGtfFile, $sGtfFile, $sGtfFileList);
my ($sOutDir, $sSampleName, $sGroupName, $sFile, $sPrefix);
my ($fpLST);
my ($sCmd);
my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;
my $s_flag = 0;	# Flag to determine if cuffCompare will run from CuffMerge

################################################################################
### Main
################################################################################

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'gtffile'} ...\n" : ();

$sOutDir = File::Spec->curdir();
if (defined $hCmdLineOption{'outdir'}) {
    $sOutDir = $hCmdLineOption{'outdir'};

    if (! -e $sOutDir) {
        mkdir($hCmdLineOption{'outdir'}) ||
            die "ERROR! Cannot create output directory\n";
    }
    elsif (! -d $hCmdLineOption{'outdir'}) {
            die "ERROR! $hCmdLineOption{'outdir'} is not a directory\n";
    }
}
$sOutDir = File::Spec->canonpath($sOutDir);

# Determine if input file is a GTF or list of GTFs and handle appropriately
($_, $_, $sFile) = File::Spec->splitpath($hCmdLineOption{'gtffile'});
# Create a list in this output_directory to use as input
$sGtfFileList = "${sOutDir}/cuffmerge_input.list";
open LIST, ">".$sGtfFileList or die "Error! Cannot open $sGtfFileList for writing: $!\n";
if ($sFile =~ m/\.gtf$/) {
	# If input is GTF file, then symlink to output_dir and write path to list
	$sCmd = "ln -sf $hCmdLineOption{'gtffile'} $sOutDir/$sFile";
	exec_command($sCmd);
	$sPrefix = $hCmdLineOption{'prefix'};
	print LIST "${sOutDir}/${sFile}\n";
}
else {
	# If input is list file, then symlink all GTF files to output_dir and write paths to list
	open ($fpLST, "$hCmdLineOption{'gtffile'}") or die "Error! Cannot open $hCmdLineOption{'gtffile'} for reading !!!\n";
	while (<$fpLST>) {
		$_ =~ s/\s+$//;
		next if ($_ =~ /^#/);
		next if ($_ =~ /^$/);
		$sGtfFile = $_;
		($_, $_, $sFile) = File::Spec->splitpath($sGtfFile);
		$sCmd = "ln -sf $sGtfFile $sOutDir/$sFile";
		exec_command($sCmd);
		print LIST "${sOutDir}/${sFile}\n";
	}
	close($fpLST);
	($_, $_, $sPrefix) = File::Spec->splitpath($hCmdLineOption{'gtffile'});
	$sPrefix =~ s/.list$//; 
}
close LIST;

($bDebug || $bVerbose) ? 
	print STDERR "\nMerging Cufflinks Transcript Analysis ...\n" : ();

if ((defined $hCmdLineOption{'ref_sequence'}) && ($hCmdLineOption{'ref_sequence'} !~ m/^$/)) {
	($bDebug || $bVerbose) ? 
	    print STDERR "Reference FASTA files found.  Will call CuffCompare through CuffMerge.\n" : ();
	$s_flag = 1;
}

($bDebug || $bVerbose) ? 
	print STDERR "\nInitiating Cuffmerge Analysis ...\n" : ();

# Add bin_dir to PATH due to some CuffMerge dependencies requiring this
local $ENV{PATH} = $hCmdLineOption{'cufflinks_bin_dir'}.":$ENV{PATH}";

$sCmd  = $hCmdLineOption{'cufflinks_bin_dir'}."/cuffmerge".
		 " -o ".$sOutDir."/".$sPrefix;
$sCmd .= " -g ".$hCmdLineOption{'annotation'} if ((defined $hCmdLineOption{'annotation'}) && ($hCmdLineOption{'annotation'} !~ m/^$/));
$sCmd .= " -s ".$hCmdLineOption{'ref_sequence'} if $s_flag; 
$sCmd .= " ".$hCmdLineOption{'args'} if (defined $hCmdLineOption{'args'});
$sCmd .= " ". $sGtfFileList;

#print($sCmd);
exec_command($sCmd);

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'gtffile'} ... done\n" : ();


################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input parameters are provided
    if (! (defined $phOptions->{'gtffile'}) ) {
		pod2usage( -msg => $sHelpHeader, -exitval => 1);
	}
	
    ## handle some defaults
    $phOptions->{'cufflinks_bin_dir'} = BIN_DIR if (! (defined $phOptions->{'cufflinks_bin_dir'}) );
    $phOptions->{'prefix'} = "cuffmrg" if ((!(defined $hCmdLineOption{'prefix'})) || ($hCmdLineOption{'prefix'} !~ m/^$/));
}

sub exec_command {
	my $sCmd = shift;
	
	if ((!(defined $sCmd)) || ($sCmd eq "")) {
		die "\nSubroutine::exec_command : ERROR! Incorrect command!\n";
	}
	
	my $nExitCode;
	
	print STDERR "$sCmd\n";
	$nExitCode = system("$sCmd");
	if ($nExitCode != 0) {
		die "\tERROR! Command Failed!\n\t$!\n";
	}
	print STDERR "\n";
}

################################################################################
