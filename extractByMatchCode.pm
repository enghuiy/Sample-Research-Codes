#!/usr/bin/env perl

# enghui yap - Jan 2013
# extract atom position according to functional atom definition

package extractByMatchCode;

use strict;
use warnings;
use Exporter;
use loadCoordFiles;

our @ISA= qw( Exporter );
our @EXPORT = qw(extractFA_RC extractFA_RC_atoms extractFA extractFA_allcases extractFA_allcases_fromArrays);

my %AA = (ALA=>'A',TYR=>'Y',MET=>'M',LEU=>'L',CYS=>'C',GLY=>'G',
         ARG=>'R',ASN=>'N',ASP=>'D',GLN=>'Q',GLU=>'E',HIS=>'H',TRP=>'W',
         LYS=>'K',PHE=>'F',PRO=>'P',SER=>'S',THR=>'T',ILE=>'I',VAL=>'V');
 
my %aa = reverse %AA;

# wrapper subroutine to handle comma-separated matchcode formats
sub extractFA_allcases_fromArrays {
    my $refx_in = shift;
    my $refy_in = shift;
    my $refz_in = shift;
    my $refAtomname_in = shift;
    my $refResno_in = shift;
    my $refResidue_in = shift;

    my $refx_out = shift;
    my $refy_out = shift;
    my $refz_out = shift;
    my $refresno_out = shift;
    my $refresidue_out = shift;
    my $refatom_out  = shift;

    my $matchCodeIn = shift;
    my $bRC_atom = shift;

    my @matchCodes = split(",", $matchCodeIn);
    
    @$refx_out=(); @$refy_out=(); @$refz_out=(); @$refresno_out=(); @$refresidue_out=(); @$refatom_out=();

    foreach my $matchCode (@matchCodes) {

	my (@x,@y,@z,@resno,@residue,@atom);
	
	if    ( $matchCode =~ /RC/ ) { 
	    if ( $bRC_atom ) { extractFA_RC_atoms_fromArrays ( $refx_in,$refy_in,$refz_in,$refAtomname_in,$refResno_in,$refResidue_in,
							       \@x,\@y,\@z,\@resno,\@residue,\@atom,$matchCode); }
	    else { extractFA_RC_fromArrays ($refx_in,$refy_in,$refz_in,$refAtomname_in,$refResno_in,$refResidue_in,
					    \@x,\@y,\@z,\@resno,\@residue,\@atom,$matchCode); }
	}
	elsif ( $matchCode =~ /MCP/ ) {
	    extractFA_MCNO_fromArrays ( $refx_in,$refy_in,$refz_in,$refAtomname_in,$refResno_in,$refResidue_in,
					\@x,\@y,\@z,\@resno,\@residue,\@atom,$matchCode); }
	
	else  { extractFA_fromArrays ( $refx_in,$refy_in,$refz_in,$refAtomname_in,$refResno_in,$refResidue_in,
				       \@x,\@y,\@z,\@resno,\@residue,\@atom,$matchCode); }
    
	push (@$refx_out,@x); 
	push (@$refy_out,@y); 	
	push (@$refz_out,@z); 
	push (@$refresno_out,@resno);
	push (@$refresidue_out,@residue);
	push (@$refatom_out,@atom);

    }
}

# subroutines to handle each kind of functional atom definition
sub extractFA_RC_fromArrays {

    my $refx_in = shift;
    my $refy_in = shift;
    my $refz_in = shift;
    my $refAtomname_in = shift;
    my $refResno_in = shift;
    my $refResidue_in = shift;

    my $refx_out = shift;
    my $refy_out = shift;
    my $refz_out = shift;
    my $refresno_out =shift;
    my $refresidue_out = shift;
    my $refatom_out  = shift;

    my $matchCode = shift;

    my $nInstance = 0;

# parse the matchCode
    $matchCode =~ s/^.*_|X//g; 
    ( $matchCode =~ /F|W|Y|P|H/) || die "Error: ring matchCode $matchCode does not contain F,W,Y,P, or H\n";
    my $matchResline="";
    foreach my $char (split //, $matchCode) { 
	if ($char=~/H/) { $matchResline.="HIS HIE HID ";}
	else {$matchResline.="$aa{$char} "; } 

    }
#    print "matching $matchResline\n";

    my $nRing = ($matchCode =~ /F|W|Y/ ? 6 : 5);

# find instances of matches
    my $ringatomcount=0;
    my $currentResno=-1;
    my $ringx=0;my $ringy=0;my $ringz=0;

    for (my $k=0;$k<scalar(@$refx_in); $k++) {

	my $res     = $refResidue_in->[$k]; 

	if ( $matchResline !~ /$res/) {next; } # <=======
	my $atom    = $refAtomname_in->[$k]; 
	my $resno   = $refResno_in->[$k];

	if ($ringatomcount == 0 ) { $currentResno = $resno; }
	if ( $resno == $currentResno &&  (($res =~ /TRP/ && $atom =~ /C[EZ]|CD2|CH2/ ) || ($res !~ /TRP/ && $atom =~ /C[DEGZ]/) || ($res=~/PRO/ && $atom=~/C[ABGD]|N/) || ($res=~/HI/ && $atom=~/CG|[CN]D|[CN]E/)) ){

	    $ringx += $refx_in->[$k];
	    $ringy += $refy_in->[$k];
	    $ringz += $refz_in->[$k];
	    $ringatomcount ++;
	    if($ringatomcount == $nRing ) {
		my $avg; 
		$refx_out->[$nInstance] = sprintf '%.3f', $ringx / $nRing;
		$refy_out->[$nInstance] = sprintf '%.3f', $ringy / $nRing;
		$refz_out->[$nInstance] = sprintf '%.3f', $ringz / $nRing;

		$refresno_out->[$nInstance] = $resno;
		$refresidue_out->[$nInstance] = $res;
		$refatom_out->[$nInstance] = " RC ";
		#print "$nInstance $refx->[$nInstance] $refy->[$nInstance] $refz->[$nInstance]\n";
		$nInstance++;

		# reset
		$ringatomcount = 0; 
		$currentResno = -1;
		$ringx = 0; $ringy =0 ; $ringz =0;
	    } 
	}
	elsif($resno != $currentResno) {
	    # reset
	    $ringatomcount = 0; 
	    $currentResno = -1;
	    $ringx = 0; $ringy =0 ; $ringz =0;
	}
    } # k
    
}

#==================================
sub extractFA_RC_atoms_fromArrays {

    my $refx_in = shift;
    my $refy_in = shift;
    my $refz_in = shift;
    my $refAtomname_in = shift;
    my $refResno_in = shift;
    my $refResidue_in = shift;

    my $refx_out = shift;
    my $refy_out = shift;
    my $refz_out = shift;
    my $refresno_out = shift;
    my $refresidue_out =shift;
    my $refatom_out  = shift;

    my $matchCode = shift;

    my $nInstance = 0;

# parse the matchCode
    $matchCode =~ s/^.*_|X//g; 
    ( $matchCode =~ /F/ || $matchCode =~ /W/ || $matchCode =~ /Y/ ) || die "Error: ring matchCode $matchCode does not contain F,W, or Y\n";
    my $matchResline="";
    foreach my $char (split //, $matchCode) { $matchResline.="$aa{$char} "; }
#    print "matching $matchResline\n";

# find instances of matches

    for (my $k=0;$k<scalar(@$refx_in); $k++) {

	my $res     = $refResidue_in->[$k]; 
	if ( $matchResline !~ /$res/) {next; } # <=======
	my $atom    = $refAtomname_in->[$k]; 
	my $resno   = $refResno_in->[$k];

	if ( ($res =~ /TRP/ && $atom =~ /C[EZ]|CD2|CH2/ ) || ($res !~ /TRP/ && $atom =~ /C[DEGZ]/) ) {

	    $refx_out->[$nInstance] = $refx_in->[$k];
	    $refy_out->[$nInstance] = $refy_in->[$k];
	    $refz_out->[$nInstance] = $refz_in->[$k];
	    $refresno_out->[$nInstance] = $resno;
	    $refresidue_out->[$nInstance] = $res;
	    $refatom_out->[$nInstance] = $atom;
	    #print "$n  $aai[$n] $fatomi[$n] $x[$nInstance] $y[$nInstance] $z[$nInstance]\n";
	    $nInstance++;
	}
    } # k
    
}

# geometric center of mainchain O-N 
# NO1 = gc of N1,O1
# NO  = gc of O1,N2
# NO2 = gc of N2,O2
sub extractFA_MCNO_fromArrays {


    my $refx_in = shift;
    my $refy_in = shift;
    my $refz_in = shift;
    my $refAtomname_in = shift;
    my $refResno_in = shift;
    my $refResidue_in = shift;

    my $refx_out = shift;
    my $refy_out = shift;
    my $refz_out = shift;
    my $refresno_out =shift;
    my $refresidue_out = shift;
    my $refatom_out  = shift;

    my $matchCode = shift;

    ( $matchCode =~ /NO[12]_|^NO_/) || die "unrecognized matchcode $matchCode\n"; 
    (my $fa =  $matchCode) =~ s/_.*$//;

    
    my $nInstance = 0;

    # figure out residue no range
    my $residno_offset=$refResno_in->[0];
    my $residno_max = $refResno_in->[-1];
    my $residno_arraySize = $residno_max - $residno_offset + 1;

   # indexArray to track N,O atom entries
    # initialize to -1 so if the atom doesn't exist, we know
    my @index_N = (-1) x $residno_arraySize; 
    my @index_O = (-1) x $residno_arraySize;
    my @index_res = (-1) x $residno_arraySize;

# extract out N,O, and residue indices
    for (my $k=0;$k<scalar(@$refx_in); $k++) {

	my $res     = $refResidue_in->[$k]; 
	my $atom    = $refAtomname_in->[$k]; $atom =~ s/^\s+|\s+$//g;
	my $resno   = $refResno_in->[$k];

	my $kk = $resno - $residno_offset; # use resid as index

	if ($atom =~ /^O$/) { $index_O[$kk] = $k; $index_res[$kk] = $k; } 
	if ($atom =~ /^N$/) { $index_N[$kk] = $k; $index_res[$kk] = $k; } 
    }

# find instances of matches
    my $cx=0;my $cy=0;my $cz=0;

    my ($rStart,$rEnd); 

    for (my $r=0; $r<scalar(@index_res)-1;$r++) {
	
	my ($kN,$kO); 
	if ( $matchCode =~ /NO1_/) {
	    $kO = $index_O[$r];
	    $kN = $index_N[$r];
	}
	elsif  ( $matchCode =~ /^NO_/) {
	    $kO = $index_O[$r];
	    $kN = $index_N[$r+1];
	}
	elsif  ( $matchCode =~ /^NO2_/) {
            $kO = $index_O[$r+1];
            $kN = $index_N[$r+1];
        }
	else { print "error in MCNO\n"; exit; }
	# skip if N or O atom(s) missing
        if ( $kO == -1 || $kN == -1 ) { next; }

	my $residueO = $refResidue_in->[ $kO ];
	my $residueN = $refResidue_in->[ $kN ];
	my $resno_O = $refResno_in->[$kO];

	# do not include MCON where N is from PROLINE
	if ( $residueN =~ /PRO/) { next; }

	$cx = ($refx_in->[$kO] + $refx_in->[$kN] ) / 2;
	$cy = ($refy_in->[$kO] + $refy_in->[$kN] ) / 2 ; 
	$cz = ($refz_in->[$kO] + $refz_in->[$kN] ) / 2;
	
	$refx_out->[$nInstance] = sprintf '%.3f', $cx;
	$refy_out->[$nInstance] = sprintf '%.3f', $cy;
	$refz_out->[$nInstance] = sprintf '%.3f', $cz;
	
	$refresno_out->[$nInstance] = $resno_O;
	$refresidue_out->[$nInstance] = $residueO;
	$refatom_out->[$nInstance] = sprintf "%4s",$fa;
	
	$nInstance++;
	
    } # r
    
}

#==========================================
sub extractFA_fromArrays {

    my $refx_in = shift;
    my $refy_in = shift;
    my $refz_in = shift;
    my $refAtomname_in = shift;
    my $refResno_in = shift;
    my $refResidue_in = shift;

    my $refx_out = shift;
    my $refy_out = shift;
    my $refz_out = shift;
    my $refresno_out = shift;
    my $refresidue_out = shift;
    my $refatom_out  = shift;

    my $matchCode = shift;

    my $nInstance = 0;

# parse the matchCode
    my (@matchAtom,@matchRes);
    (my $resCode = $matchCode ) =~ s/^.*\_//;

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DEFINE MATCHCODES HERE

    # if only one residue, use the exact atom name given
    if (length($resCode) == 1 ) { 
	( $matchAtom[0] = $matchCode ) =~ s/\_.*$//;
	if ($resCode =~ /H/) {
	    $matchRes[0] = "HI"; # allows for HIS/HIE/HID
	}
	else {	
	    $matchRes[0] = $aa{$resCode};
	}
    }

    elsif($matchCode =~ /NH_RX/)  {
        $matchAtom[0]="NH1";  $matchRes[0]="ARG";
        $matchAtom[1]="NH2";  $matchRes[1]="ARG";
    }
    elsif($matchCode =~ /NS_RX/)  {
        $matchAtom[0]="NE";   $matchRes[0]="ARG";
        $matchAtom[1]="NH1";  $matchRes[1]="ARG";
        $matchAtom[2]="NH2";  $matchRes[2]="ARG";
    }

    elsif($matchCode =~ /NS_HX/)  {
        $matchAtom[0]="ND1";  $matchRes[0]="HI";
        $matchAtom[1]="NE2";  $matchRes[1]="HI";
    }
    elsif($matchCode =~ /CG_VX/)  {
        $matchAtom[0]="CG1";  $matchRes[0]="VAL";
        $matchAtom[1]="CG2";  $matchRes[1]="VAL";
    }

    elsif($matchCode =~ /OS_DX/)  {
	$matchAtom[0]="OD1";  $matchRes[0]="ASP";
        $matchAtom[1]="OD2";  $matchRes[1]="ASP";
    }
    elsif($matchCode =~ /OS_EX/)  {
	$matchAtom[0]="OE1";  $matchRes[0]="GLU";
        $matchAtom[1]="OE2";  $matchRes[1]="GLU";
    }
#===============================================
# LEVEL 1
    elsif($matchCode =~ /N_KRX/)  {
	$matchAtom[0]="NE";   $matchRes[0]="ARG";
	$matchAtom[1]="NH1";  $matchRes[1]="ARG";
	$matchAtom[2]="NH2";  $matchRes[2]="ARG";
	$matchAtom[3]="NZ";   $matchRes[3]="LYS";
    }

    elsif($matchCode =~ /N_NQX/)  {
        $matchAtom[0]="ND2";  $matchRes[0]="ASN";
        $matchAtom[1]="NE2";  $matchRes[1]="GLN";
    }
    elsif($matchCode =~ /O_NQX/)  {
	$matchAtom[0]="OD1";  $matchRes[0]="ASN";
	$matchAtom[1]="OE1";  $matchRes[1]="GLN";
    }
    elsif($matchCode =~ /O_DEX/)  { 
	$matchAtom[0]="OD1";  $matchRes[0]="ASP";
	$matchAtom[1]="OD2";  $matchRes[1]="ASP";
	$matchAtom[2]="OE1";  $matchRes[2]="GLU";
	$matchAtom[3]="OE2";  $matchRes[3]="GLU";
    }
    elsif($matchCode =~ /C_ILVX/)  {
	$matchAtom[0]="CG1";  $matchRes[0]="ILE";
	$matchAtom[1]="CG";   $matchRes[1]="LEU";
        $matchAtom[2]="CG1";  $matchRes[2]="VAL";
        $matchAtom[3]="CG2";  $matchRes[3]="VAL";
    }
    elsif($matchCode =~ /OG_STX/)  {
        $matchAtom[0]="OG";   $matchRes[0]="SER";
        $matchAtom[1]="OG1";  $matchRes[1]="THR";
    }

    elsif($matchCode =~ /O_MCX/)  {
	$matchAtom[0]="O";   $matchRes[0]="...";
    }
    elsif($matchCode =~ /N_MCX/)  {
	$matchAtom[0]="N";   $matchRes[0]="...";
    }

#===============================================
# LEVEL 2 (?)
    elsif($matchCode =~ /ODE_DEQ/)  { 
	$matchAtom[0]="OD1";  $matchRes[0]="ASP";
	$matchAtom[1]="OD2";  $matchRes[1]="ASP";
	$matchAtom[2]="OE1";  $matchRes[2]="GLU";
	$matchAtom[3]="OE2";  $matchRes[3]="GLU";
	$matchAtom[4]="OE1";  $matchRes[4]="GLN";
    }

    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    else { print "Error : matchCode $matchCode not defined!!!!\n"; exit; }

# find instances
    for (my $k=0;$k<scalar(@$refx_in); $k++) {

	my $res     = $refResidue_in->[$k]; 
	my $atom    = $refAtomname_in->[$k]; $atom =~ s/^\s+|\s+$//g; 
	my $resno   = $refResno_in->[$k];
	my $bMatch = 0;
	
	for(my $i=0;$i<scalar(@matchRes); $i++) { 
            if ( $res=~ /$matchRes[$i]/ && $atom =~ /^$matchAtom[$i]$/ ) { $bMatch = 1; last; }
	}
	if ( ! $bMatch ) {next; }

	$refx_out->[$nInstance] = $refx_in->[$k]; $refx_out->[$nInstance]  =~ s/^\s+|\s+$//g;
	$refy_out->[$nInstance] = $refy_in->[$k]; $refy_out->[$nInstance]  =~ s/^\s+|\s+$//g;
	$refz_out->[$nInstance] = $refz_in->[$k]; $refz_out->[$nInstance]  =~ s/^\s+|\s+$//g;
	$refresno_out->[$nInstance] = $resno;
	$refresidue_out->[$nInstance] = $res;
	$refatom_out->[$nInstance] = $atom;
	#print "$n  $aai[$n] $fatomi[$n] $x[$nInstance] $y[$nInstance] $z[$nInstance]\n";
	$nInstance++;

	
    } # k
}
#===================================================================================================
# wrappers to use with pdbfile
#===================================================================================================

# handle comma-separated matchcode formats
sub extractFA_allcases {
    
    my $pdbfile = shift;
    my $refx_out = shift;
    my $refy_out = shift;
    my $refz_out = shift;
    my $refresno_out = shift;
    my $refresidue_out = shift;
    my $refatom_out  = shift;
    my $matchCodeIn = shift;
    my $bRC_atom = shift;
    
    my (@x_in,@y_in,@z_in,@atomname_in,@resno_in,@residue_in);
    
    _loadPDB_xyz_atomRes ($pdbfile,\@x_in,\@y_in,\@z_in,\@atomname_in,\@resno_in,\@residue_in);
    
    extractFA_allcases_fromArrays (\@x_in,\@y_in,\@z_in,\@atomname_in,\@resno_in,\@residue_in,
				   $refx_out,$refy_out,$refz_out,$refresno_out,$refresidue_out,$refatom_out,
				   $matchCodeIn,$bRC_atom);
}


sub extractFA_RC {

    my $pdbfile = shift;
    my $refx_out = shift;
    my $refy_out = shift;
    my $refz_out = shift;
    my $refresno_out = shift;
    my $refresidue_out = shift;
    my $refatom_out  = shift;
    my $matchCodeIn = shift;

    my (@x_in,@y_in,@z_in,@atomname_in,@resno_in,@residue_in);
    
    _loadPDB_xyz_atomRes ($pdbfile,\@x_in,\@y_in,\@z_in,\@atomname_in,\@resno_in,\@residue_in);
    
    extractFA_RC_fromArrays (\@x_in,\@y_in,\@z_in,\@atomname_in,\@resno_in,\@residue_in,
		  $refx_out,$refy_out,$refz_out,$refresno_out,$refresidue_out,$refatom_out,
		  $matchCodeIn);
}

sub extractFA {

    my $pdbfile = shift;
    my $refx_out = shift;
    my $refy_out = shift;
    my $refz_out = shift;
    my $refresno_out = shift;
    my $refresidue_out = shift;
    my $refatom_out  = shift;
    my $matchCodeIn = shift;

    my (@x_in,@y_in,@z_in,@atomname_in,@resno_in,@residue_in);
    
    _loadPDB_xyz_atomRes ($pdbfile,\@x_in,\@y_in,\@z_in,\@atomname_in,\@resno_in,\@residue_in);
    
    extractFA_fromArrays (\@x_in,\@y_in,\@z_in,\@atomname_in,\@resno_in,\@residue_in,
			  $refx_out,$refy_out,$refz_out,$refresno_out,$refresidue_out,$refatom_out,
			  $matchCodeIn);
}

sub extractFA_RC_atoms {

    my $pdbfile = shift;
    my $refx_out = shift;
    my $refy_out = shift;
    my $refz_out = shift;
    my $refresno_out = shift;
    my $refresidue_out = shift;
    my $refatom_out  = shift;
    my $matchCodeIn = shift;

    my (@x_in,@y_in,@z_in,@atomname_in,@resno_in,@residue_in);
    
    _loadPDB_xyz_atomRes ($pdbfile,\@x_in,\@y_in,\@z_in,\@atomname_in,\@resno_in,\@residue_in);
    
    extractFA_RC_atoms_fromArrays (\@x_in,\@y_in,\@z_in,\@atomname_in,\@resno_in,\@residue_in,
				      $refx_out,$refy_out,$refz_out,$refresno_out,$refresidue_out,$refatom_out,
				      $matchCodeIn);
}
1;
