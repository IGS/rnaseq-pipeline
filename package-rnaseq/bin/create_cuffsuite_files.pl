#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

create_cuffsuite_files.pl - script to create cuffcompare, cuffdiff, and cuffquant input files.

=head1 SYNOPSIS

    create_cuffsuite_files.pl --i input_samples_file --p cuffcompare|cuffdiff --g comparison_groups
                              [--gtf reference_gtf] [--c gtf_listfile] [--s sam_listfile] [--q cxb_listfile] 
							  [--o outdir] [--v] 

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS

    --i <input_samples_file>       = /path/to/input tab-delimited file with the following sample information.
                                     #Sample_Name    Group_Name    [Alignment_SAM_File/Cufflinks_GTF_File]

    --p <cuffcompare|cuffdiff|cuffquant>     = Cuffsuite program (Cuffcompare, Cuffdiff, or Cuffquant) to create files for.

    --g <comparison_groups>        = string of groups to compare. e.g. "GRP#2vsGRP#1,GRP#3vsGRP#1".
                                     Group Ids SHOULD match Group Names in column 2 of the sample info file

    --gtf <reference_gtf>          = Single reference GTF annotation to be used with all cuffdiff comparisons.

    --c <gtf_listfile>             = /path/to/listfile of cufflinks or cuffcompare GTF files.
                                     The cufflinks filename should be in the format /path/to/<sample_name>*.gtf. Optional if specified in the sample file itself.
                                     The cuffcompare filename should be in the format /path/to/GRP#2vsGRP#1*.gtf.

    --s <sam_listfile>             = /path/to/listfile of alignment SAM files sorted by position.
                                     The SAM alignment filename should be in the format /path/to/<sample_name>*.sam.
                                     Optional if specified in the sample file itself.

	--q <cxb_listfile>             = /path/to/listfile of abundance CXB files obtained from Cuffquant
	                                 Optional if specified in the sample file itself.

    --o <output dir>               = /path/to/output directory. Optional. [present working directory]

    --v                            = generate runtime messages. Optional

=head1 DESCRIPTION

The script creates the input files for the cuffcompare, cuffdiff, and cuffquant ergatis components.

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
use FindBin qw($RealBin);

##############################################################################
### Constants
##############################################################################

use constant FALSE => 0;
use constant TRUE  => 1;

use constant CUFFLINKS_BIN_DIR => '/usr/local/bin';

use constant VERSION => '1.0.0';
use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM." version ".VERSION."\n";

GetOptions( \%hCmdLineOption,
            'infile|i=s', 'gtflist|c=s', 'samlist|s=s', 'cuff_prog|p=s', 'comp_grps|g=s', 'gtf=s', 
            'cxblist|q=s', 'outdir|o=s', 'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

my (%hSamples, %hGroups);
my (@aGTFList, @aComparisons);
my ($sOutDir, $sSampleName, $sGroupName, $sCountsFile, $sInFile, $sFile);
my ($sGrpX, $sGrpY, $sCGrp, $sGTFFile);
my ($fpSMPL, $fpLST, $fpOUT);
my ($sCmd, $nI);
my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

################################################################################
### Main
################################################################################

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'infile'} ...\n" : ();

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

($bDebug || $bVerbose) ? 
	print STDERR "\nGenerating input files for Cuffsuite Analysis ...\n" : ();

open ($fpSMPL, "$hCmdLineOption{'infile'}") or die "Error! Cannot open $hCmdLineOption{'infile'} for reading !!!\n";

while (<$fpSMPL>) {
	$_ =~ s/\s+$//;
	
	next if ($_ =~ /^#/);
	
	next if ($_ =~ /^$/);
	# Assign sample group name and file to a sample name key
	($sSampleName, $sGroupName, $sFile) = (split(/\t/, $_))[0, 1, 2];
	
	$hSamples{$sSampleName} = [$sGroupName, undef];
	
	$hSamples{$sSampleName}[1] = $sFile if (defined $sFile);
	
	# Assign sample names to a group name key
	if (!(defined $hGroups{$sGroupName})) {
		my @aSampleArray = ();
		$hGroups{$sGroupName} = \@aSampleArray;
	}
	push @{$hGroups{$sGroupName}}, $sSampleName;
}

close($fpSMPL);

if ($bDebug) {
	foreach $sSampleName (sort keys %hSamples) {
		print STDERR "\t$sSampleName".
					 "\t$hSamples{$sSampleName}[0]".
					 "\t$hSamples{$sSampleName}[1]\n";
	}
}

if ($hCmdLineOption{'cuff_prog'} =~ m/^cuffcompare$/i) {
	if ((defined $hCmdLineOption{'gtflist'}) && ($hCmdLineOption{'gtflist'} !~ m/^$/)) {
		($bDebug || $bVerbose) ? 
			print STDERR "\nExtracting input GTF files for Cuffcompare Analysis ...\n" : ();
		
		open ($fpLST, "$hCmdLineOption{'gtflist'}") or die "Error! Cannot open $hCmdLineOption{'gtflist'} for reading !!!\n";
		
		while (<$fpLST>) {
			$_ =~ s/\s+$//;
			
			next if ($_ =~ /^#/);
			
			next if ($_ =~ /^$/);
			
			$sInFile = $_;
			
			($_, $_, $sFile) = File::Spec->splitpath($sInFile);
						
			foreach $sSampleName (sort keys %hSamples) {
				($bDebug) ? print STDERR "\t##### $sSampleName\n" : ();
				if ($sFile =~ m/^$sSampleName\./) {
					$hSamples{$sSampleName}[1] = $sInFile;
					($bDebug) ? print STDERR "\t##### $sSampleName\t$sInFile\n" : ();
					last;
				}
			}
		}
		
		close($fpLST);
	}
	
	if ($bDebug) {
		foreach $sSampleName (sort keys %hSamples) {
			print STDERR "\t$sSampleName".
						 "\t$hSamples{$sSampleName}[0]".
						 "\t$hSamples{$sSampleName}[1]\n";
		}
	}
	
	@aComparisons = split(/,/, $hCmdLineOption{'comp_grps'});
	
	open($fpLST, ">$sOutDir/cuffsuite_input_file.list") or die "Error! Cannot open $sOutDir/cuffcompare_sample_info.list for writing: $!";
	
	foreach $sCGrp (@aComparisons) {
		($sGrpX, $sGrpY) = split(/vs/, $sCGrp);
		
		open($fpSMPL, ">$sOutDir/$sCGrp.sample.list") or die "Error! Cannot open $sOutDir/$sCGrp.sample.list for writing: $!";
		
		foreach $sSampleName (@{$hGroups{$sGrpY}}) {
			print $fpSMPL "$hSamples{$sSampleName}[1]\n";
		}
		
		foreach $sSampleName (@{$hGroups{$sGrpX}}) {
			print $fpSMPL "$hSamples{$sSampleName}[1]\n";
		}
		
		close($fpSMPL);
		
		print $fpLST "$sOutDir/$sCGrp.sample.list\n";
	}
	
	close($fpLST);
}
elsif ($hCmdLineOption{'cuff_prog'} =~ m/^cuffdiff$/i) {
	# Determine if SAM or CBX is the input
	die "Error!  Cannot have both a SAM/BAM list and a CBX list as input.  Please choose one!\n" if (defined $hCmdLineOption{'samlist'} && defined $hCmdLineOption{'cxblist'});

	my $input_list;
	if ((defined $hCmdLineOption{'samlist'}) && ($hCmdLineOption{'samlist'} !~ m/^$/)) {
	  $input_list = $hCmdLineOption{'samlist'};
	}
	if ((defined $hCmdLineOption{'cxblist'}) && ($hCmdLineOption{'cxblist'} !~ m/^$/)) {
	  $input_list = $hCmdLineOption{'cxblist'};
	}

	if (defined $input_list) {
		# Open list and associate file with sample in hash
		open ($fpLST, $input_list) or die "Error! Cannot open $input_list for reading !!!\n";
		
		while (<$fpLST>) {
			$_ =~ s/\s+$//;
			
			next if ($_ =~ /^#/);
			
			$sInFile = $_;
			
			($_, $_, $sFile) = File::Spec->splitpath($sInFile);
			
			foreach $sSampleName (sort keys %hSamples) {
				if ($sFile =~ m/^$sSampleName\./) {
					$hSamples{$sSampleName}[1] = $sInFile;
					last;
				}
			}
		}
		
		close($fpLST);
	}
	# If GTF list is provided, instead of single GTF file, open and push GTF files into array
	if (! defined $hCmdLineOption{'gtf'}) {
		open ($fpLST, "$hCmdLineOption{'gtflist'}") or die "Error! Cannot open $hCmdLineOption{'gtflist'} for reading !!!\n";
		
		@aGTFList = ();
		while (<$fpLST>) {
			next if ($_ =~ /^#/);
			push @aGTFList, $_;
		}
		
		close($fpLST);
	}
	
	@aComparisons = split(/,/, $hCmdLineOption{'comp_grps'});
	
	open($fpLST, ">$sOutDir/cuffsuite_input_file.list") or die "Error! Cannot open $sOutDir/cuffsuite_input_file.list for writing: $!";
	# For each comparison group...
	foreach $sCGrp (@aComparisons) {
		($bDebug || $bVerbose) ? 
			print STDERR "\tGenerating input files for $sCGrp ...\n" : ();
		
		# If GTF list was provided, use GTF file associated with current comparison group
		if (defined $hCmdLineOption{'gtf'}) {
			$sGTFFile = $hCmdLineOption{'gtf'};
		} elsif (scalar @aGTFList == 1) {
			# Handle the samtools merged.gtf file case, in particulara
            $sGTFFile = shift @aGTFList;
		} else {
			foreach $sInFile (@aGTFList) {
				
				($_, $_, $sFile) = File::Spec->splitpath($sInFile);
				
				if ($sFile =~ m/^$sCGrp/) {
					$sGTFFile = $sInFile;
					last;
				}
			}
		}
		
		# Split the group into its X and Y parts
		($sGrpX, $sGrpY) = split(/vs/, $sCGrp);
		
		open($fpSMPL, ">$sOutDir/$sCGrp.sample.list") or die "Error! Cannot open $sOutDir/$sCGrp.sample.list for writing: $!";
		
		# Join files of the same group part by comma
		foreach $sSampleName (sort @{$hGroups{$sGrpY}}) {
			print $fpSMPL "$hSamples{$sSampleName}[1],";
		}
		
		print $fpSMPL "\n";
		
		foreach $sSampleName (sort @{$hGroups{$sGrpX}}) {
			print $fpSMPL "$hSamples{$sSampleName}[1],";
		}
		
		print $fpSMPL "\n";
		
		close($fpSMPL);
		# Add the GTF file for the comparison group to the line
		print $fpLST "$sOutDir/$sCGrp.sample.list\t$sGTFFile\n";
	}
	
	close($fpLST);
} 
elsif ( $hCmdLineOption{'cuff_prog'} =~ m/^cuffquant$/i) {
	# If list of SAM files is provided, open list and associate file with sample in hash
	if ((defined $hCmdLineOption{'samlist'}) && ($hCmdLineOption{'samlist'} !~ m/^$/)) {
		open ($fpLST, "$hCmdLineOption{'samlist'}") or die "Error! Cannot open $hCmdLineOption{'samlist'} for reading !!!\n";
		
		while (<$fpLST>) {
			$_ =~ s/\s+$//;
			next if ($_ =~ /^#/);
			$sInFile = $_;
			($_, $_, $sFile) = File::Spec->splitpath($sInFile);
			foreach $sSampleName (sort keys %hSamples) {
				if ($sFile =~ m/^$sSampleName\./) {
					$hSamples{$sSampleName}[1] = $sInFile;
					last;
				}
			}
		}
		close($fpLST);
	}
	# If GTF list is provided, instead of single GTF file, open and push GTF files into array
	if (! defined $hCmdLineOption{'gtf'}) {
		open ($fpLST, "$hCmdLineOption{'gtflist'}") or die "Error! Cannot open $hCmdLineOption{'gtflist'} for reading !!!\n";
		
		@aGTFList = ();
		while (<$fpLST>) {
			next if ($_ =~ /^#/);
			push @aGTFList, $_;
		}
		close($fpLST);
	}
	
	@aComparisons = split(/,/, $hCmdLineOption{'comp_grps'});
	
	open($fpLST, ">$sOutDir/cuffsuite_input_file.list") or die "Error! Cannot open $sOutDir/cuffsuite_input_file.list for writing: $!";
	
    my %seen_groups = ();
	# For each comparison group...
	foreach $sCGrp (@aComparisons) {
		($bDebug || $bVerbose) ? 
			print STDERR "\tGenerating input files for $sCGrp ...\n" : ();
		
		# If using a GTF list, find the GTF associated with the current comparison group
		if (defined $hCmdLineOption{'gtf'}) {
			$sGTFFile = $hCmdLineOption{'gtf'};
		} elsif (scalar @aGTFList == 1) {
			# Handle the samtools merged.gtf file case, in particulara
            $sGTFFile = shift @aGTFList;
		} else {
			foreach $sInFile (@aGTFList) {
				($_, $_, $sFile) = File::Spec->splitpath($sInFile);
				if ($sFile =~ m/^$sCGrp/) {
					$sGTFFile = $sInFile;
					last;
				}
			}
		}
		
		# Split the group into both parts, get the files for both parts, and write to output list.
		($sGrpX, $sGrpY) = split(/vs/, $sCGrp);
		
        # Only print group information if group was seen
        unless (defined $seen_groups{$sGrpY}){
            foreach $sSampleName (sort @{$hGroups{$sGrpY}}) {
                print $fpLST "$hSamples{$sSampleName}[1]\t$sGTFFile\n";
            }
            $seen_groups{$sGrpY} = 1;
        }
		
        unless (defined $seen_groups{$sGrpX}){
            foreach $sSampleName (sort @{$hGroups{$sGrpX}}) {
                print $fpLST "$hSamples{$sSampleName}[1]\t$sGTFFile\n";
            }
            $seen_groups{$sGrpX} = 1;
        }
	}
	
	close($fpLST);
}	
else {
	print STDERR "\nIncorrect Cufflinks program stated .....\n";
	pod2usage( -msg => $sHelpHeader, -exitval => 1);
}

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'infile'} ... done\n" : ();


################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input parameters are provided
    if ((! (defined $phOptions->{'infile'}) ) || 
    	(! (defined $phOptions->{'cuff_prog'}) ) || 
    	(! (defined $phOptions->{'comp_grps'}) )) {
		pod2usage( -msg => $sHelpHeader, -exitval => 1);
	}
	
	if (($hCmdLineOption{'cuff_prog'} !~ m/^cuff(compare|diff|quant)$/i)) {
		print STDERR "\nIncorrect Cufflinks program stated .....\n";
		pod2usage( -msg => $sHelpHeader, -exitval => 1);
	}
}

################################################################################
