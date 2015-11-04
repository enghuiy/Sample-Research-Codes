#!/usr/bin/perl

# use module
use XML::Simple;
use Data::Dumper;

# script takes in a joblist containing (pdbfile, and chain1,chain2 to search)
# and a XML output file from PISA website and outputs a summary file 
# and a list per joblist line,per chain listing the interface residues

(@ARGV == 2 ) || die "@ARGV; Usage: [interface.pisa file] [ joblist (pdb,chain1,chain2) ]\n";

# open joblist that contains the pdb, and chain1,chain2 to search
my $xmlfile = $ARGV[0];
my $jobfile = $ARGV[1];
my $summaryfile = "summary.out";

# create object and read XML file
$xml = new XML::Simple;
$data = $xml->XMLin($xmlfile);


open (SUMFILE,">$summaryfile");

foreach my $myPdbEntry (@{$data->{pdb_entry}}) {

  my $pdbCode = $myPdbEntry->{pdb_code}; 
  my $nInterfaces = $myPdbEntry->{n_interfaces};
  print "parsing $pdbCode in XML file ...\n";

# check if this PDB is in the joblist
  open(JOBFILE,$jobfile);
  if (! grep{/$pdbCode/} <JOBFILE>) {close JOBFILE; next; } else {close JOBFILE; }

  for (my $n=1; $n<=$nInterfaces; $n++) {
    $myInterface = $myPdbEntry->{interface}->{$n}; 


    my $molecule1 = $myInterface->{molecule}->{1};
    my $molecule2 = $myInterface->{molecule}->{2};

# skip if the symmetry operation on molecule 2 is not XYZ
    if($molecule2->{symop} !~ /[X|x],[Y|y],[Z|z]/) {next;}

# skip if this PDB/chain1-chain2 combi is not in the joblist

    my $chainID1 = $molecule1->{chain_id};
    my $chainID2 = $molecule2->{chain_id};

    #print "interface $n: chains $chainID1 $chainID2\n";

    open(JOBFILE,$jobfile);
    @jobdata=<JOBFILE>; close JOBFILE;
    my ($found1,$found2) = (0,0);
    $found1 =  grep(/$pdbCode $chainID1 $chainID2/, @jobdata);
    if (! $found1) { $found2 =  grep(/$pdbCode $chainID2 $chainID1/,@jobdata); }
    if (! $found1 && ! $found2 ) {next;}   

# write the interface stats to a summary file
# write out the residue index and type to a separate file

    my $nIntRes1 = $molecule1->{int_nres};
    my $nIntRes2 = $molecule2->{int_nres};

    my ($interfaceFile1,$interfaceFile2);
    if($found1) {
      print SUMFILE "$pdbCode $n $chainID1 $chainID2 $nIntRes1 $nIntRes2\n";
      $interfaceFile1 = "$pdbCode.$chainID1-$chainID2.$chainID1.dat";
      $interfaceFile2 = "$pdbCode.$chainID1-$chainID2.$chainID2.dat";
    }
    else {
      print SUMFILE"$pdbCode $n $chainID2 $chainID1 $nIntRes2 $nIntRes1\n";
      $interfaceFile1 = "$pdbCode.$chainID2-$chainID1.$chainID1.dat";
      $interfaceFile2 = "$pdbCode.$chainID2-$chainID1.$chainID2.dat";
    }

# interface residues = those with non-zero BSA

    open(INTFILE, ">>$interfaceFile1");
    foreach my $myRes (@{$molecule1->{residues}->{residue} }) {
      if($myRes->{bsa} > 0) {
	my $seqNum = $myRes->{seq_num};
	my $ASA = $myRes->{asa};
	my $resName =  $myRes->{resname};
	print INTFILE "$seqNum $resName $ASA\n";	
      }
    }
    close INTFILE;

    open(INTFILE, ">>$interfaceFile2");
    foreach my $myRes (@{$molecule2->{residues}->{residue} }) {
      if($myRes->{bsa} > 0) {
	my $seqNum = $myRes->{seq_num};
	my $ASA = $myRes->{asa};
	my $resName =  $myRes->{resname};
	print INTFILE "$seqNum $resName $ASA\n";	
      }
    }
    close INTFILE;

    
  } # end-interface-n

}# end-pdb_entry

close SUMFILE;
# print output
#print Dumper($data);

# uniq interface files using command line
open(JOBFILE,$jobfile);
while (<JOBFILE>) {
  my ($pdbCode, $chainID1, $chainID2) = split(" ");
  $file1 = "$pdbCode.$chainID1-$chainID2.$chainID1.dat";
  $file2 = "$pdbCode.$chainID1-$chainID2.$chainID2.dat";
  system ("sort -n $file1 | uniq > temp.out; mv temp.out $file1; ");
  system ("sort -n $file2 | uniq > temp.out; mv temp.out $file2; ");
}
