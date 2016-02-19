#!/usr/bin/perl
use strict;
use warnings;
$SIG{CHLD} = 'IGNORE';

#MAKE SURE CHILD DIES IF KEYBOARD INTERUP SIGNAL
$SIG{INT} = \&signal_handler;

#MAKE SURE CHILD DIES IF KEYBOARD TERMINATE SIGNAL
$SIG{TERM} = \&signal_handler;

#First argument is the Fasta database file, second argument is the results from SRST2, should start in the same directory as the scores file.
use File::Glob;
use POSIX;
use Getopt::Long;
use Cwd 'abs_path';
use Cwd;

#Get path to script
my $script_path = abs_path($0);
my $script_dir = substr($script_path, 0, rindex($script_path, "/"));
my $working_dir = getcwd();

#Get user input
my $fastq_directory;
my $scoresName;
my $fasta_input;
my $definitions_file;
my $gene_fasta;
my $MLST_input;
my $MLST_definitions_file;
my $recN_input;
my $virulence_input;
my $forward;
my $reverse;
my $SingleOrPaired;

my %opt=qw();
GetOptions(\%opt, "help|h!", "fastq_directory:s", "scoreName:s", "serotype_db:s", "serotype_definitions:s", "cps2K:s", "MLST_db:s", "MLST_definitions:s", "recN_db:s", "Virulence_db:s", "forward:s", "reverse:s", "ends:s");

if (exists($opt{"help"})) {
	print "**********************************************\nImplementation of the emm pipeline:\n**********************************************\nperl Ssuis_serotypingPipeline.pl --fastq_directory /path/to/fastq/directory --scoreName Scores_output_name --serotype_db serotype.fasta -- serotype_definitions serotype_definitions.txt --cps2K cps2K.fasta\n\n";
	print "--fastq_directory\tPath to directory containing paired-end fastq files.\n\t\t\tMust be full path to directory, please do not use '.'\n\t\t\tor '..' to declare path\n\t\t\tIf no path is given, the current working\n\t\t\tdirectory is used\n";
	print "--scoreName\t\tName of SRST2 results file\n\t\t\t[optional: default name 'Results']\n";
	print "--serotype_db\t\tMultifasta file containing the serotype database\n\t\t\t(Ssuis_Serotyping.fasta)\n\t\t\t[If none is provided, Ssuis_Serotyping.fasta is looked\n\t\t\tfor in the directory containing the script]\n";
	print "--serotype_definitions\tText file containing the definitions for the\n\t\t\tserotype database file\n\t\t\t(Ssuis_Serotyping_Definitions.txt)\n\t\t\t[If none is provided, Ssuis_Serotyping_Definitions.txt is looked\n\t\t\tfor in the directory containing the script]\n";
	print "--cps2K\t\tMultifasta file containing the cpsH confirmation database\n\t\t\t(Ssuis_cps2K.fasta)\n\t\t\t[If none is provided, Ssuis_cps2K.fasta is looked\n\t\t\tfor in the directory containing the script]\n";
	print "--MLST_db\t\tMultifasta file containing the MLST database\n\t\t\t(Streptococcus_suis.fasta)\n\t\t\t[If none is provided, Streptococcus_suis.fasta is looked\n\t\t\tfor in the directory containing the script]\n";
	print "--MLST_definitions\tText file containing the definitions for the\n\t\t\tMLST database file\n\t\t\t(ssuis.txt)\n\t\t\t[If none is provided, ssuis.txt is looked\n\t\t\tfor in the directory containing the script]\n";
	print "--recN_db\t\tFasta file containing the recN species specfic gene\n\t\t\t(recN_Full.fasta)\n\t\t\t[If none is provided, recN_full.fasta is looked\n\t\t\tfor in the directory containing the script]\n";
	print "--Virulence_db\t\tMultifasta file containing the Virulence genes\n\t\t\t(Virulence.fasta)\n\t\t\t[If none is provided, Virulence.fasta is looked\n\t\t\tfor in the directory containing the script]\n";
	print "--forward\t\tIndicator delimiting the forward reads file for\n\t\t\tpaired read fastq files\nThis option is ignored if single-end reads is selected\n\t\t\t[optional: default '_R1']\n";
	print "--reverse\t\tIndicator delimiting the reverse reads file for\n\t\t\tpaired read fastq files\nThis option is ignored if single-end reads is selected\n\t\t\t[optional: default '_R2']\n\n";
	print "--ends\t\tIndicates whether the reads are paired-end (pe) or single-end (se)\n\t\t\t[optional: default 'pe']\n\n";
	print "Note: We recommend using paired end reads of at least 100nt in length\nand at least 30X coverage.\nWe have not tested the efficiency of the pipeline with reads shorter than 80nt.";
	exit;
}

if(exists($opt{"fastq_directory"})){
	$fastq_directory=$opt{"fastq_directory"};
}
else{
	$fastq_directory=$working_dir;
	my $any_fastqs=glob("*.fastq");

	if(!defined($any_fastqs)){	
		print "Please provide the full path to the directory containing fastq files to use. See help file [--help]";
		exit;
	}
}

if(exists($opt{"scoreName"})){
	$scoresName=$opt{"scoreName"};
}
else{
	$scoresName = "Results";
}

if(exists($opt{"serotype_db"})){
	$fasta_input=$opt{"serotype_db"};
}
else{
	$fasta_input=$script_dir."/Ssuis_Serotyping.fasta";
	unless(-e $fasta_input){
		print "Please provide the serotype database fasta file. See help file [--help]";
		exit;
	}
}

if(exists($opt{"serotype_definitions"})){
	$definitions_file=$opt{"serotype_definitions"};
}
else{
	$definitions_file=$script_dir."/Ssuis_Serotyping_Definitions.txt";
	unless(-e $definitions_file){
		print "Please provide the serotype definitions file. See help file [--help]";
		exit;
	}
}

if(exists($opt{"cps2K"})){
	$gene_fasta=$opt{"cps2K"};
}
else{
	$gene_fasta=$script_dir."/Ssuis_cps2K.fasta";
	unless(-e $gene_fasta){
		print "Please provide the cps2K database fasta file. See help file [--help]";
		exit;
	}
}

if(exists($opt{"MLST_db"})){
	$MLST_input=$opt{"MLST_db"};
}
else{
	$MLST_input=$script_dir."/Streptococcus_suis.fasta";
	unless(-e $MLST_input){
		print "Please provide the MLST database fasta file. See help file [--help]";
		exit;
	}
}

if(exists($opt{"MLST_definitions"})){
	$MLST_definitions_file=$opt{"MLST_definitions"};
}
else{
	$MLST_definitions_file=$script_dir."/ssuis.txt";
	unless(-e $MLST_definitions_file){
		print "Please provide the MLST definitions file. See help file [--help]";
		exit;
	}
}

if(exists($opt{"recN_db"})){
	$recN_input=$opt{"recN_db"};
}
else{
	$recN_input=$script_dir."/recN_full.fasta";
	unless(-e $recN_input){
		print "Please provide the recN fasta file. See help file [--help]";
		exit;
	}
}

if(exists($opt{"Virulence_db"})){
	$virulence_input=$opt{"Virulence_db"};
}
else{
	$virulence_input=$script_dir."/Virulence.fasta";
	unless(-e $virulence_input){
		print "Please provide the Virulence database fasta file.  See help file [--help]";
		exit;
	}
}

if(exists($opt{"forward"})){
	$forward = $opt{"forward"};
}
else{
	$forward = "_R1";
}
if(exists($opt{"reverse"})){
	$reverse = $opt{"reverse"};
}
else{
	$reverse = "_R2";
}

if(exists($opt{"ends"})){
	$SingleOrPaired = $opt{"ends"};
}
else{
	$SingleOrPaired = "pe";
}

if(($SingleOrPaired ne "pe") and ($SingleOrPaired ne "se")){
	die "Ends must be either pe or se";
}

#CHECK IF STRAIN ACTUALLY BELONGS TO THE SPECIES STREPTOCOCCUS SUIS
mkdir "recN";
chdir "recN";

if($SingleOrPaired eq "pe"){
	system("srst2.py --input_pe $fastq_directory/*.fastq --forward $forward --reverse $reverse --output $scoresName\_recN --log --gene_db $recN_input --forward $forward --reverse $reverse --save_scores");
}
elsif($SingleOrPaired eq "se"){
	system("srst2.py --input_se $fastq_directory/*.fastq --output $scoresName\_recN --log --gene_db $recN_input --save_scores");
}

my $recNResults = glob("$scoresName\_recN__genes*.txt");

open my $recNs, "<$recNResults" or die "Can't open recN results file!";

my @ssuis;
my @nonssuis;
my $recCount = 0;
foreach my $recN (<$recNs>){
	$recN =~ s/\r?\n//;

	if($recCount > 0){
		my @info = split(/\s+/, $recN);

		if(scalar(@info) > 1){
			if(substr($info[1], 0, 4) eq "recN"){
				if($SingleOrPaired eq "pe"){
					push(@ssuis, "$fastq_directory/$info[0]$forward.fastq");
					push(@ssuis, "$fastq_directory/$info[0]$reverse.fastq");
				}
				elsif($SingleOrPaired eq "se"){
					push(@ssuis, "$fastq_directory/$info[0]*.fastq");
				}
			}
			else{push(@nonssuis, $info[0])};
		}
	}

	$recCount++;
}

my $ss = join(' ', @ssuis);

#ORGANIZE SPECIES CONFIRMATION OUTPUT
system("mkdir pileups");
system("mkdir sorted_bam");
system("mkdir scores");
system("mv *.pileup ./pileups");
system("mv *.sorted.bam ./sorted_bam");
system("mv *.scores ./scores");

close($recNs);
system("mv $recNResults $scoresName\_speciesConfirmation.txt");

#Check if there are any S. suis in the dataset
if(scalar(@ssuis) < 1){
	print "No S. suis in the dataset";
	exit;
}

#ORGANIZE SYSTEM COMMANDS
my @commands;
if($SingleOrPaired eq "pe"){
	@commands = ("srst2.py --input_pe $ss --forward $forward --reverse $reverse --output $scoresName\_MLST --log --mlst_db $MLST_input --mlst_definitions $MLST_definitions_file --save_scores", "srst2.py --input_pe $ss --forward $forward --reverse $reverse --output $scoresName\_VirulenceFactors --log --gene_db $virulence_input --save_scores", "srst2.py --input_pe $ss --forward $forward --reverse $reverse --output $scoresName --log --mlst_db $fasta_input --mlst_definitions $definitions_file --save_scores");
}
elsif($SingleOrPaired eq "se"){
	@commands = ("srst2.py --input_se $ss --output $scoresName\_MLST --log --mlst_db $MLST_input --mlst_definitions $MLST_definitions_file --save_scores", "srst2.py --input_se $ss --output $scoresName\_VirulenceFactors --log --gene_db $virulence_input --save_scores", "srst2.py --input_se $ss --output $scoresName --log --mlst_db $fasta_input --mlst_definitions $definitions_file --save_scores");
}

chdir "..";
mkdir "MLST";
chdir "MLST";

my $pid1 = fork();
if ($pid1 == 0){
	setpgrp;
	system($commands[0]);
	exit 0;
}

#WAIT 1 SECOND, THEN CHECK IF MLST IS STILL RUNNING, IF NOT -> KILL PROGRAM
sleep(1);
my $exists = kill 0, $pid1;
if($exists == 0){
	print "Unable to run MLST!";
	exit;
}

chdir "..";
mkdir "Virulence";
chdir "Virulence";

#RUN VIRULENCE PIPELINE
my $pid2 = fork();
if ($pid2 == 0){
	setpgrp;
	system($commands[1]);
	exit 0;
}

#WAIT 1 SECOND, THEN CHECK IF MLST IS STILL RUNNING, IF NOT -> KILL PROGRAM
sleep(1);
my $exists2 = kill 0, $pid2;
if($exists2 == 0){
	print "Unable to run virulence check!";
	exit;
}

chdir "..";
mkdir "Serotyping";
chdir "Serotyping";

#RUN SEROTYPING PIPELINE
system($commands[2]);

my $ResultsName = glob("$scoresName\__mlst__*.txt");

#Checks if output file has any results
if ((-s $ResultsName) < 50){
	print "Can't run Serotyping!";
	&signal_handler;
}

#Organize SRST2 output
system("mkdir pileups");
system("mkdir sorted_bam");
system("mkdir scores");
system("mv *.pileup ./pileups");
system("mv *.sorted.bam ./sorted_bam");
system("mv *.scores ./scores");

#Read in results of SRST2
open my $srst2_results, "<$ResultsName" or die "Can't open Serotyping results file!";

mkdir "Pipeline";
open my $EndResults, ">", "$scoresName\_FinalSerotypingResults.txt" or die $!;

print $EndResults "Strain\tSerotype\n";

#READING IN SEROTYPE CALLING RESULTS#
my $resultCount = 0;
my @furtherAnalysis_1;
my @furtherAnalysis_2;
my @notTesting;
foreach my $result (<$srst2_results>){
	$result =~ s/\r?\n//;

	if($resultCount > 0){
		my @reading = split(/\t/, $result);
		#IF SEROTYPE IS CALLED AS A 1, PUT IN ARRAY TO CHECK IF IT IS A 1 OR A 14#
		if(($reading[1] eq "1") || ($reading[1] eq "1*") || ($reading[1] eq "1?") || ($reading[1] eq "1*?")){
			push(@furtherAnalysis_1, $reading[0]);
		}
		#IF SEROTYPE IS CALLED AS A 2, PUT IN ARRAY TP CHECK IF IT IS A 2 OR A 1/2#
		elsif(($reading[1] eq "2") || ($reading[1] eq "2*") || ($reading[1] eq "2?") || ($reading[1] eq "2*?")){
			push(@furtherAnalysis_2, $reading[0]);
		}
		#IF SEROTYPE IS NEITHER A 1 OR A 2, PRINT RESULTS TO FILE#
		else{
			push(@notTesting, ($reading[0] . "\t" . $reading[1]));
			print $EndResults "$reading[0]\t$reading[1]\n";
		}
	}
	$resultCount++;
}

close($srst2_results);
system("mv $ResultsName $scoresName\_InitialCapsuleResults.txt");

chdir "Pipeline";

#CREATE BOWTIE AND SAMTOOLS INDEX FILES FOR CPSK FASTA FILE#
system("bowtie2-build $gene_fasta $gene_fasta");
system("samtools faidx $gene_fasta");

#FOR ALL OF THE 1s ALIGN TO CPSK FASTA TO SEE IF IT IS A 1 OR A 14#
foreach my $one (@furtherAnalysis_1){
	#Run Bowtie2 against cpsK
	if($SingleOrPaired eq "pe"){
		system("bowtie2 -1 $fastq_directory/$one$forward.fastq $fastq_directory/$one$reverse.fastq -S $one\_vs_cpsK.sam --very-sensitive-local --no-unal -a -x $gene_fasta");
	}
	elsif($SingleOrPaired eq "se"){
		system("bowtie2 -U $fastq_directory/$one*.fastq -S $one\_vs_cpsK.sam --very-sensitive-local --no-unal -a -x $gene_fasta");
	}
	#Run Samtools to call SNPs
	system("samtools view -b -o $one\_vs_cpsK.bam -q 1 -S $one\_vs_cpsK.sam");
	system("samtools sort $one\_vs_cpsK.bam $one\_vs_cpsK.sorted");
	system("samtools mpileup -u -L 1000 -f $gene_fasta -Q 20 -q 1 $one\_vs_cpsK.sorted.bam > $one\_vs_cpsK.pileup");
	system("bcftools view -vcg $one\_vs_cpsK.pileup > $one\_vs_cpsK.raw.vcf");
	system("vcfutils.pl varFilter -Q 20 $one\_vs_cpsK.raw.vcf > $one\_vs_cpsK.vcf");

	open my $vcf, "<$one\_vs_cpsK.vcf" or die "Can't find input file!";
	open my $vcf_out, ">", "$one\_vs_cpsK.snps" or die $!;

	#Re-write SNP file to a more readable format
	foreach my $vcf_line (<$vcf>){
		if(substr($vcf_line, 0, 2) ne "##"){
			print $vcf_out $vcf_line;
		}
	}

	open my $SNP_out, ">", "$one\_SNPeffect.txt" or die $!;

	#Check AminoAcid changes caused by SNPs
	system("perl $script_dir/SNP_AminoAcidChange.pl $one\_vs_cpsK.snps $gene_fasta > $one\_SNPeffect.txt");

	#READ IN SNP EFFFECT FILE TO CHECK IF THERE IS A TRP OR CYS AT AMINO ACID POSITION 161#
	open my $SNP_read, "<$one\_SNPeffect.txt" or die "Can't open input file!";

	my @SNP_array;
	my $SNP_line = 0;
	my $Target = 0;
	my $AminoAcid;
	foreach my $snp(<$SNP_read>){
		if($SNP_line > 0){
		        @SNP_array = split(/\t/, $snp);
			
			if($SNP_array[5] == 161){
				$Target = 1;
				$AminoAcid = $SNP_array[7];
			}
		}
		$SNP_line++;
	}

	#IF THERE WAS A SNP AT POSITION 161, CHECK IF IT IS A C, PRINT TO FILE#	
	if($Target == 1){
		if($AminoAcid eq "C"){
			print $EndResults "$one\t1\n";
			push(@notTesting, ($one . "\t1"));
		}
		else{
			print $EndResults "$one\t14\n";
			push(@notTesting, ($one . "\t14"));
		}
	}
	else{
		print $EndResults "$one\t14\n";
		push(@notTesting, ($one . "\t14"));
	}
	
}

#FOR ALL OF THE 2s ALIGN TO CPSK FASTA TO SEE IF IT IS A 2 OR A 1/2#
foreach my $two (@furtherAnalysis_2){
	#Run Bowtie against cpsK
	if($SingleOrPaired eq "pe"){
		system("bowtie2 -1 $fastq_directory/$two$forward.fastq $fastq_directory/$two$reverse.fastq -S $two\_vs_cpsK.sam --very-sensitive-local --no-unal -a -x $gene_fasta");
	}
	elsif($SingleOrPaired eq "se"){
		system("bowtie2 -U $fastq_directory/$two*.fastq -S $two\_vs_cpsK.sam --very-sensitive-local --no-unal -a -x $gene_fasta");		
	}

	#Run Samtools to call SNPs
	system("samtools view -b -o $two\_vs_cpsK.bam -q 1 -S $two\_vs_cpsK.sam");
	system("samtools sort $two\_vs_cpsK.bam $two\_vs_cpsK.sorted");
	system("samtools mpileup -u -L 1000 -f $gene_fasta -Q 20 -q 1 $two\_vs_cpsK.sorted.bam > $two\_vs_cpsK.pileup");
	system("bcftools view -vcg $two\_vs_cpsK.pileup > $two\_vs_cpsK.raw.vcf");
	system("vcfutils.pl varFilter -Q 20 $two\_vs_cpsK.raw.vcf > $two\_vs_cpsK.vcf");

	open my $vcf, "<$two\_vs_cpsK.vcf" or die "Can't find input file!";
	open my $vcf_out, ">", "$two\_vs_cpsK.snps" or die $!;

	#Re-write SNP file to a more readable format
	foreach my $vcf_line (<$vcf>){
		if(substr($vcf_line, 0, 2) ne "##"){
			print $vcf_out $vcf_line;
		}
	}

	open my $SNP_out, ">", "$two\_SNPeffect.txt" or die $!;

	#Check AminoAcid changes caused by SNPs
	system("perl $script_dir/SNP_AminoAcidChange.pl $two\_vs_cpsK.snps $gene_fasta > $two\_SNPeffect.txt");

	#READ IN SNP EFFFECT FILE TO CHECK IF THERE IS A TRP OR CYS AT AMINO ACID POSITION 161#
	open my $SNP_read, "<$two\_SNPeffect.txt" or die "Can't open input file!";

	my @SNP_array;
	my $SNP_line = 0;
	my $Target = 0;
	my $AminoAcid;
	foreach my $snp(<$SNP_read>){
		if($SNP_line > 0){
		        @SNP_array = split(/\t/, $snp);
			if($SNP_array[5] == 161){
				$Target = 1;
				$AminoAcid = $SNP_array[7];
			}
		}
		$SNP_line++;
	}

	#IF THERE WAS A SNP AT POSITION 161, CHECK IF IT IS A C, PRINT TO FILE#		
	if($Target == 1){
		if($AminoAcid eq "C"){
			print $EndResults "$two\t1/2\n";
			push(@notTesting, ($two . "\t1/2"));
		}
		else{
			print $EndResults "$two\t2\n";
			push(@notTesting, ($two . "\t2"));
		}
	}
	else{
		print $EndResults "$two\t2\n";
		push(@notTesting, ($two . "\t2"));
	}
	
}

#SORT OUTPUT FILES#
system("mkdir sam");
system("mkdir unsorted_bam");
system("mkdir sorted_bam");
system("mkdir pileup");
system("mkdir raw_vcf");
system("mkdir filtered_vcf");
system("mkdir snps");
system("mkdir snp_effect");
system("mv *.sam ./sam");
system("mv *.sorted.bam ./sorted_bam");
system("mv *.bam ./unsorted_bam");
system("mv *.pileup ./pileup");
system("mv *.raw.vcf ./raw_vcf");
system("mv *.vcf ./filtered_vcf");
system("mv *.snps ./snps");
system("mv *SNPeffect.txt ./snp_effect");

#WAIT FOR MLST TO FINISH BEFORE MOVING ON TO ORGANIZATION STEP
wait();

chdir "../../MLST";

#Organize SRST2 output
system("mkdir pileups");
system("mkdir sorted_bam");
system("mkdir scores");
system("mv *.pileup ./pileups");
system("mv *.sorted.bam ./sorted_bam");
system("mv *.scores ./scores");

my $mlst_Name = glob("$scoresName\_MLST__mlst__*.txt");

#Read in results of SRST2
open my $mlst_results, "<$mlst_Name" or die "Can't open MLST results file!";

#READING IN SEROTYPE CALLING RESULTS#
my $mlstresultCount = 0;
my @mlst;
foreach my $stresult (<$mlst_results>){
	$stresult =~ s/\r?\n//;

	if($mlstresultCount > 0){
		my @reading = split(/\t/, $stresult);
		push(@mlst, \@reading);
	}
	$mlstresultCount++;
}

#Change name of MLST output
close $mlst_results;
system("mv $mlst_Name $scoresName\_MLSTResults.txt");

chdir "../Virulence";

#Organize Virulence output
system("mkdir pileups");
system("mkdir sorted_bam");
system("mkdir scores");
system("mv *.pileup ./pileups");
system("mv *.sorted.bam ./sorted_bam");
system("mv *.scores ./scores");

my $vir_name = glob("$scoresName\_VirulenceFactors__genes*.txt");

#Read in Virulence factor gene results
open my $vir_results, "<$vir_name" or die "Can't open Virulence results files!";

#READING IN VIRULENCE RESULTS
my $virulenceResultsCount = 0;
my @virulence;
my @virHeader;
my $virNames;
my @virRead;
my @virStrain;
my @virWrite;
foreach my $virResult (<$vir_results>){
	$virResult =~ s/\r?\n//;

	if($virulenceResultsCount == 0){
		@virHeader = split(/\t/, $virResult);
		shift(@virHeader);
		$virNames = join("\t", @virHeader);
	}
	else{
		@virRead = split(/\t/, $virResult);
		push(@virStrain, $virRead[0]);
		shift(@virRead);
		push(@virWrite, join("\t", @virRead));
	}

	$virulenceResultsCount++;
}

my $factorNum = scalar(split(/\t/, $virNames));

close($vir_results);
system("mv $vir_name $scoresName\_VirulenceFactorResults.txt");

chdir "..";

#OPEN OUTPUT FILE
open my $combined_results, ">", "$scoresName\_FinalResults.txt" or die $!;
print $combined_results "Strain\tSerotype\tST\taroA\tcpn60\tdpr\tgki\tmutS\trecA\tthrA\trecN\t$virNames\n";

#NEW READ IN STRAINS
my $strainNum = 0;
foreach my $readings (@virStrain) {
	print $combined_results $readings;

	#FIND SEROTYPE RESULT MATCHING VIRULENCE STRAIN NAME
	my $seroMatch = 0;
	my $SeroCount = 0;
	my @sero_split;

	while(($seroMatch == 0) && ($SeroCount < scalar(@notTesting))){
		@sero_split = split(/\t/, $notTesting[$SeroCount]);
		if($sero_split[0] eq $readings){
			$seroMatch = 1;
			print $combined_results "\t$sero_split[1]";
		}
		else{$SeroCount++};
	}
	if($seroMatch == 0){print $combined_results "\t-"};

	#FIND MLST RESULT MATCHING VIRULENCE STRAIN NAME
	my $mlst_match = 0;
	my $mlst_count = 0;
	my $mlst_scalar = scalar(@mlst);

	while(($mlst_match == 0) && ($mlst_count < $mlst_scalar)){
		if($mlst[$mlst_count][0] eq $readings){
			$mlst_match = 1;
			print $combined_results "\t$mlst[$mlst_count][1]\t$mlst[$mlst_count][2]\t$mlst[$mlst_count][3]\t$mlst[$mlst_count][4]\t$mlst[$mlst_count][5]\t$mlst[$mlst_count][6]\t$mlst[$mlst_count][7]\t$mlst[$mlst_count][8]";
		}
		else{$mlst_count++};
	}
	if($mlst_match == 0){print $combined_results "\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF"};

	#PRINT recN+ and Virulence results
	if(substr($virWrite[$strainNum], 0, 1) eq "?"){
		my $vir_rep = "\t-" x $factorNum;
		print $combined_results "\t+$vir_rep\n";
	}
	else{
		print $combined_results "\t+\t$virWrite[$strainNum]\n";
	}

	$strainNum++;
}

if(scalar(@nonssuis) > 0){
	my $vir_rep = "\tN/A" x $factorNum;
	foreach my $non (@nonssuis){
		print $combined_results "$non\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\t-$vir_rep\n";
	}
}

#MAKES SURE CHILD DIES IF PARENT IS INTERUPTED
sub signal_handler {
	my $exists_kill = kill 0, $pid1;
	if ($exists_kill == 1){
		kill -9, $pid1;
	}
	die;
};
