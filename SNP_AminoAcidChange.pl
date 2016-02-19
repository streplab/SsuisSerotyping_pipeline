use Bio::Perl;
use Bio::SeqIO;
use Bio::Seq;
use Bio::DB::RefSeq;


##READ SNP FILE NAME FROM COMMANDLINE ARGUMENTS
$snp_file=$ARGV[0];

##READ REFERENCE SEQUENCE FROM COMMAND LINE ARGUMENTS
$filename=$ARGV[1];


$gb= new Bio::DB::RefSeq;

$input_sequence = Bio::SeqIO->new( -file => "<$filename");

##READ SEQ OBJECT FROM STREAM##
my $seq_object = $input_sequence->next_seq;

###OPEN INPUT SNP FILE
open (SNPFILE, $snp_file) or die "Can't find SNP file!";

###########LOAD SNPFILE CONTENTS INTO ARRAY################
foreach $line (<SNPFILE>) {
	if(substr($line, 0 , 1) ne "#"){
		@row=(split(/\t/, $line));
	        push @SNParray, [@row];
	}
};

print "Position\tRef_base\tSample_base\tRef_Codon\tSample_Codon\tAA_Position\tRef_AA\tSample_AA\tSubstitution\n";

foreach $line(@SNParray) {
	$codon="";
	$effect="";
	$amino_acid="";
	$ref_amino_acid="";
	@row=@$line;

	$sample_base = $row[4];
	$ref_base = $row[3];
	$gene_pos = $row[1];
	$codon_pos = ($gene_pos % 3);
	$aminoAcid_pos;

	$len=(length $sample_base);
	$ref_len=(length $ref_base);
	if(($len==1) && ($ref_len==1)){	
		if($codon_pos==0) {
			$codon_seq=$seq_object->subseq(($gene_pos-2),($gene_pos));
			$mut_seq=($seq_object->subseq(($gene_pos-2),($gene_pos-1))).($sample_base);
			$codon_obj=Bio::Seq->new(-seq =>$mut_seq, -alphabet => 'dna' );
			$codon=$codon_obj->seq;
			$aminoAcid_pos = ($gene_pos/3);
		};
		if($codon_pos==1) {
			$codon_seq=$seq_object->subseq(($gene_pos),($gene_pos+2));
			$mut_seq=($sample_base).($seq_object->subseq(($gene_pos+1),($gene_pos+1))).($seq_object->subseq(($gene_pos+2),($gene_pos+2)));
			$codon_obj=Bio::Seq->new(-seq =>$mut_seq, -alphabet => 'dna' );
			$codon=$codon_obj->seq;
			$aminoAcid_pos = (($gene_pos+2)/3);
			};
		if($codon_pos==2) {
			$codon_seq=$seq_object->subseq(($gene_pos-1),($gene_pos+1));
			$mut_seq=($seq_object->subseq(($gene_pos-1),($gene_pos-1))).($sample_base).($seq_object->subseq(($gene_pos+1),($gene_pos+1)));
			$codon_obj=Bio::Seq->new(-seq =>$mut_seq, -alphabet => 'dna' );
			$codon=$codon_obj->seq;
			$aminoAcid_pos = (($gene_pos+1)/3);
		};
		###TRANSLATE CODON TO AMINO ACID CODES
		$amino_acid=translate_as_string($codon);
		$ref_amino_acid=translate_as_string($codon_seq);
		if ($amino_acid eq $ref_amino_acid) {
			$effect="synonymous";
		
		}
		else {
			$effect="nonsynonymous";
		};
	}
	print "$gene_pos\t$ref_base\t$sample_base\t$codon_seq\t$codon\t$aminoAcid_pos\t$ref_amino_acid\t$amino_acid\t$effect\n"
}
