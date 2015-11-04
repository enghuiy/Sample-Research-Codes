#!/usr/bin/env perl

# script by enghui yap - Jan 2013
# score candidate pdb file against a interface template file
#
# i) combinatorial searches through all possible combination of
# ii) matching candidate atoms, check pairwise distances and
# iii) perform lsqfit and clash check against receptor

use strict;
use Benchmark;
use IO::Handle;

use extractByMatchCode;
my $t0 = Benchmark->new;

# main program
(@ARGV==3) || die("Usage: [pdbfile] [templatefile][recPQR]\n") ;

my $pdbfile      = $ARGV[0];
my $templatefile = $ARGV[1];
my $receptorPQRfile = $ARGV[2];

my statdir="./";
my $cutofftol    = 4.0;
my $index        = "x";
my $candidatePQRfile = $pdbfile;
my $clash_tolerance = 3.0;
my $tstarts= "0,0,0,0,0";
my $bestistartstring = "-1,-1,-1,-1,-1";

my ($t0start,$t1start,$t2start,$t3start,$t4start,$t5start,$t6start,$t7start,$t8start,$t9start)=split(",",$tstarts);
my @besti_start = split(",",$bestistartstring);

my $bVerbose = 0; 
my $bTime=0;

# load interface template 
open(TEMFILE, $templatefile) or die ("Can't open $templatefile\n");
my (@xi,@yi,@zi,@fatomi); 
my $dum;
my $n=0;
while (<TEMFILE>) {
    ($xi[$n],$yi[$n],$zi[$n],$fatomi[$n],$dum) = split(" ");
    $n++;
}
my $nFuncDef=$n;
my $denominator = $nFuncDef * ($nFuncDef-1) / 2; 

# check which lines have same masks
my @indicesSameMask =  (0) x $nFuncDef; 
for(my $n=0;$n<$nFuncDef;$n++) {
  my @indices_n = ();
  for(my $m=0;$m<$nFuncDef;$m++) {
    if ($n == $m) {next;}
    if($fatomi[$n] =~ /^$fatomi[$m]$/) { 
      push(@indices_n, $m);
    }
  }
  @indicesSameMask[$n] = [ @indices_n ];
}

# open status file for restart logging
system("mkdir -p $statdir");
open (my $statfh, ">","$statdir/status.r$nFuncDef.c$cutofftol.$index.$pdbfile.out");
autoflush $statfh 1;

print $statfh  "status dir = $statdir\n";
print $statfh  "status mindiststart = $mindiststart \n";
print $statfh  "status tstarts = $t0start,$t1start,$t2start,$t3start,$t4start,$t5start,$t6start,$t7start,$t8start,$t9start\n";
print $statfh  "status besti = ",join(" ",@besti_start ),"\n";

# populate n-def x n-def pairwise distances in reference distmat
my (@distTemplate, @distTemplatebuf );
my (@distArray,@rankI,@rankJ);
my $nRank=0;
for(my $n=0;$n<$nFuncDef;$n++) {
    push(@distTemplate,   [ (0) x $nFuncDef ] );
    push(@distTemplatebuf,[ (0) x $nFuncDef ] );
    for(my $m=$n+1;$m<$nFuncDef;$m++) {
	 $distTemplate[$n][$m]    =  &rmsd( $xi[$n],$yi[$n],$zi[$n],$xi[$m],$yi[$m],$zi[$m]);
	 $distTemplatebuf[$n][$m] =  $distTemplate[$n][$m] + $cutofftol;
	 $distArray[$nRank] = $distTemplate[$n][$m];
	 $rankI[$nRank] = $n;
	 $rankJ[$nRank] = $m;
	 $nRank++;
    }
}

# sort distArray to see which combi is the most restrictive (for culling)
my (@orderedRanks); 
my @list_order = sort { $distArray[$a] <=> $distArray[$b] } 0 .. $#distArray;
$orderedRanks[0] = $rankI[$list_order[0]];
$orderedRanks[1] = $rankJ[$list_order[0]];

for (my $i=1;$i<$#list_order; $i++) {
  my $rank = $rankI[$list_order[$i]];  if ( ! grep( /^$rank$/, @orderedRanks ) ) { push( @orderedRanks,  $rank); }
  else {
    my $rank = $rankJ[$list_order[$i]];
    if ( ! grep( /^$rank$/, @orderedRanks ))   { push( @orderedRanks,  $rank);}
      }
}
#for my $i (0..$#orderedRanks) { print $orderedRanks[$i],"\n";}
#for my $i (0..$#rankI) { print "$rankI[$i] $rankJ[$i] $distArray[$i]\n";}

#============================================================================
# extract instances of matching coordinates from pdbfile
my @nInstances = (0) x $nFuncDef; 
my (@xns,@yns,@zns,@resnosMatched,@residuesMatched,@atomsMatched);

for(my $n=0;$n<$nFuncDef;$n++) {
    my (@x,@y,@z,@resnos,@residues,@atoms);
    if    ( $fatomi[$n] =~ /RC/ )    { &extractFA_RC ( $pdbfile, \@x,\@y,\@z,\@resnos,\@residues,\@atoms, $fatomi[$n]); }
    else  { &extractFA ( $pdbfile, \@x,\@y,\@z,\@resnos,\@residues,\@atoms, $fatomi[$n]); }

    @xns[$n] = [ @x ] ;
    @yns[$n] = [ @y ] ;
    @zns[$n] = [ @z ] ;
    @resnosMatched[$n]   = [ @resnos ];
    @residuesMatched[$n] = [ @residues ];
    @atomsMatched[$n] = [ @atoms ];
    $nInstances[$n] = scalar(@x);
    if($bVerbose) {print "Number of instances for template line $n ( $fatomi[$n] ): $nInstances[$n] \n";}
    if ($nInstances[$n] == 0) { print "NA\n"; exit; }

} # n-def loop

# ===============================================================================================
# loop through all valid combinations of instances to get best dRMS

# stores the indices of the best-matching instances for each func def
my @besti;
for (my $n=0;$n<$nFuncDef;$n++) {
    $besti[$n] = $besti_start[$n];
} 

my $mindist_inf  = 1000000;
my $mindist = $mindist_inf;
my $nFiltered = 0;

for (my $t=0; $t<$nFuncDef; $t++) {
    print  $statfh "instances: $t $nInstances[ $orderedRanks[$t]] $fatomi[$orderedRanks[$t] ]\n";
}

my $drms;
if($nFuncDef == 5)     { $mindist = getBestDRMS_5N($statfh); }
else {print "Error: cannot handle $nFuncDef no. of template lines!\n"; exit; }
close $statfh;
#=======================================================================
if ($mindist < $mindist_inf) {

    $drms = sqrt ( $mindist / $denominator);
    printf "%10.5f |", $drms;
    &printBestMatchPDB("bestmatch.r$nFuncDef.c$cutofftol.$index.$pdbfile.out.pdb", \@besti);
}
else {
    printf "%10.5s |", "NA";
}
for (my $i=0;$i<$nFuncDef;$i++) { printf "%3d ", $nInstances[$i]; }
printf "| %12d | %8.3f \n",$nFiltered,$cutofftol;

my $tf = Benchmark->new;
my $td = timediff($tf, $t0); 
if ($bTime) {
    my $pstring=" computation time=".timestr($td)."\n";
    printf "%30s",$pstring;
}

#=========================================================================================================================
# SUBROUTINES
#=========================================================================================================================

sub rmsd {
    my $x1 = $_[0];
    my $y1 = $_[1];
    my $z1 = $_[2];
    my $x2 = $_[3];
    my $y2 = $_[4];
    my $z2 = $_[5];
    return sqrt( ($x1-$x2)*($x1-$x2) +  ($y1-$y2)*($y1-$y2)  + ($z1-$z2)*($z1-$z2)  );
}

sub checkSameMask {
  my $p = $_[0];
  my $r = $_[1];
  my $refSameMask = $_[2];

  my $bSame = 0;
  my @indicesSameMask_p = @{ $refSameMask->[$p] };
  for (my $i=0; $i<scalar( @indicesSameMask_p ); $i++) {
   
    if ($indicesSameMask_p[$i] == $r) { $bSame = 1; last; }

  }
  return $bSame;
}

sub checkSkip {

  my $pi    = $_[0];
  my $refx = $_[1];
  my $refy = $_[2];
  my $refz = $_[3];
  my $reft = $_[4];
  my $refSameMask = $_[5];
  my $refdisttempbuf = $_[6];
  my $refdisttemp    = $_[7];
  my $reforderedRanks = $_[8];
  my $refDiffsqP = $_[9];

  $$refDiffsqP = 0;
  my $bSkip = 0;
  my $p = $reforderedRanks->[$pi];

  for(my $ri=0;$ri<$pi;$ri++) {
    my $r = $reforderedRanks->[$ri];

# check same mask
    if ( &checkSameMask($p,$r,$refSameMask) == 1 && ($reft->[$p] == $reft->[$r]) ) {
      $bSkip = 1; 
      last; 
    }
# check distance
    my $dist = &rmsd( $refx->[$p],$refy->[$p],$refz->[$p],   $refx->[$r],$refy->[$r],$refz->[$r]  );
    my ($refdistbuf,$refdist); 
    if ($r < $p) {$refdist = $refdisttemp->[$r][$p];}
    else {$refdist = $refdisttemp->[$p][$r];}
    my $diff = $dist - $refdist;
    if (abs($diff) > $cutofftol ) {
     # print "$p $r $refx->[$p] $refx->[$r] || $dist $refdistbuf \n";
      $bSkip = 1; last;}
    else {
      $$refDiffsqP += $diff * $diff; }
    }
return $bSkip;
}

sub printBestMatchPDB {
  my $outfile = shift;
  my $ref_bestis = shift;

  # print out best match
  open (OUTPDBFILE,">$outfile");
  for (my $n=0; $n<$nFuncDef; $n++) {
    my $besti = $ref_bestis->[$n];
    printf OUTPDBFILE "%-6s %4d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f %5.2f %5.2f\n",
      "ATOM",$n,$atomsMatched[$n][$besti], $residuesMatched[$n][$besti]," ",$resnosMatched[$n][$besti],$xns[$n][$besti],$yns[$n][$besti],$zns[$n][$besti],1.0,0.0;
  }
  close OUTPDBFILE;
}

# to generate a corresponding template to bestmatch for alignment
sub printTemplatePDB_usingBestMatchAtoms {
  my $outfile = shift;
  my $ref_bestis = shift;

  # print out best match
  open (OUTPDBFILE,">$outfile");
  for (my $n=0; $n<$nFuncDef; $n++) {
    my $besti = $ref_bestis->[$n];
    printf OUTPDBFILE "%-6s %4d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f %5.2f %5.2f\n",
      "ATOM",$n,$atomsMatched[$n][$besti], $residuesMatched[$n][$besti]," ",$resnosMatched[$n][$besti],$xi[$n],$yi[$n],$zi[$n],1.0,0.0;
  }
  close OUTPDBFILE;
}

sub getBestDRMS_5N() {

    my $statfh = shift;

  $nFiltered = 0;
  my $p0=$orderedRanks[0]; my $p1=$orderedRanks[1]; 
  my $p2=$orderedRanks[2]; my $p3=$orderedRanks[3];
  my $p4=$orderedRanks[4]; 
  
  my (@xpos,@ypos,@zpos,@tind);
  for (my $t0=$t0start;$t0<$nInstances[$p0]; $t0++) {
      $xpos[$p0] = $xns[$p0][$t0]; $ypos[$p0] = $yns[$p0][$t0]; $zpos[$p0] = $zns[$p0][$t0]; $tind[$p0] = $t0; 
      
      my $bstart1 = ($t0==$t0start ? 1 : 0);
      my $t1start_tmp = ($bstart1 ? $t1start : 0);
      for (my $t1=$t1start_tmp;$t1<$nInstances[$p1]; $t1++) {
	  $xpos[$p1] = $xns[$p1][$t1]; $ypos[$p1] = $yns[$p1][$t1]; $zpos[$p1] = $zns[$p1][$t1]; $tind[$p1] = $t1; 
	  my $diffsqP1;
	  if( &checkSkip(1, \@xpos,\@ypos,\@zpos, \@tind, \@indicesSameMask,\@distTemplatebuf,\@distTemplate,\@orderedRanks,\$diffsqP1) ) {next;};
	  
	  my $bstart2 = ($t1==$t1start ? $bstart1 : 0);
	  my $t2start_tmp = ($bstart2 ? $t2start : 0);
	  for (my $t2=$t2start_tmp;$t2<$nInstances[$p2]; $t2++) {
	      $xpos[$p2] = $xns[$p2][$t2]; $ypos[$p2] = $yns[$p2][$t2]; $zpos[$p2] = $zns[$p2][$t2]; $tind[$p2] = $t2; 
	      my $diffsqP2;
	      if( &checkSkip(2, \@xpos,\@ypos,\@zpos, \@tind, \@indicesSameMask,\@distTemplatebuf,\@distTemplate,\@orderedRanks,\$diffsqP2) ) {next;};
	      
	      my $bstart3 = ($t2==$t2start ? $bstart2 : 0);
	      my $t3start_tmp = ($bstart3 ? $t3start : 0);
	      for (my $t3=$t3start_tmp;$t3<$nInstances[$p3]; $t3++) {
		  $xpos[$p3] = $xns[$p3][$t3]; $ypos[$p3] = $yns[$p3][$t3]; $zpos[$p3] = $zns[$p3][$t3]; $tind[$p3] = $t3;
		  my $diffsqP3;
		  if( &checkSkip(3, \@xpos,\@ypos,\@zpos, \@tind, \@indicesSameMask,\@distTemplatebuf,\@distTemplate,\@orderedRanks,\$diffsqP3) ) {next;};
				  
		  my $bstart4 = ($t3==$t3start ? $bstart3 : 0);
		  my $t4start_tmp = ($bstart4 ? $t4start : 0);
		  for (my $t4=$t4start_tmp;$t4<$nInstances[$p4]; $t4++) {
		      $xpos[$p4] = $xns[$p4][$t4]; $ypos[$p4] = $yns[$p4][$t4]; $zpos[$p4] = $zns[$p4][$t4]; $tind[$p4] = $t4;
		      my $diffsqP4;
		      if( &checkSkip(4, \@xpos,\@ypos,\@zpos, \@tind, \@indicesSameMask,\@distTemplatebuf,\@distTemplate,\@orderedRanks,\$diffsqP4) ) {next;};
		      
		      my $sum = $diffsqP1 + $diffsqP2 + $diffsqP3 + $diffsqP4;
		      if ($sum < $mindist) {
			  
			  # run clash check before accepting this combi
			  my @besti_trial=();
			  $besti_trial[$p0]=$t0;$besti_trial[$p1]=$t1;$besti_trial[$p2]=$t2;$besti_trial[$p3]=$t3; 
			  $besti_trial[$p4]=$t4;
			  &printBestMatchPDB("bestmatch.tmp.pdb", \@besti_trial);
			  &printTemplatePDB_usingBestMatchAtoms("template.tmp.pdb", \@besti_trial);
			  my $bClash = `run_align.clashCheck.sh template.tmp.pdb bestmatch.tmp.pdb $candidatePQRfile $receptorPQRfile $clash_tolerance | tail -1 | awk '{print $1}'`; 
			  
			  if ($bClash==1) { next; }
			  $mindist = $sum; 
			  @besti=@besti_trial;

# print the next indices to debug file
			  my ($d0,$d1,$d2,$d3,$d4);
			  $d4=$t4+1; $d3=$t3; $d2=$t2; $d1=$t1;  $d0=$t0;
			  if( $d4 > $nInstances[$p4] ) { $d4=0; $d3++; }
			  if( $d3 > $nInstances[$p3] ) { $d3=0; $d2++; }   
			  if( $d2 > $nInstances[$p2] ) { $d2=0; $d1++; }
			  if( $d1 > $nInstances[$p1] ) { $d1=0; $d0++; }

			  printf $statfh "%8.3f %4d %4d %4d %4d %4d ",$mindist, $d0,$d1,$d2,$d3,$d4;
			  for (my $n=0;$n<$nFuncDef;$n++) { printf $statfh "%4d ",$besti[$n];} 
			  print $statfh "\n";

		      }
		      
		      $nFiltered++;
		      
		  }
	      }
	  }
      }
      # at the end of current t0 loop, print out the current best for restart purpose
      
      printf $statfh "%8.3f %4d %4d %4d %4d %4d ",$mindist, $t0,$nInstances[$p2]-1,$nInstances[$p3]-1,$nInstances[$p4]-1;
      for (my $n=0;$n<$nFuncDef;$n++) { printf $statfh "%4d ",$besti[$n];}
      print $statfh " CURRENT\n";

  }
  print $statfh "DONE\n";
  return $mindist; 
}
