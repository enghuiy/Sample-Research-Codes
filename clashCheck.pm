#!/usr/bin/env perl

# functions to do clash checking
# script by enghui yap - Apr 2014
# modified Jun27 2014 - add pseudo point for N=3 alignment
# 
package clashCheck;
use strict;
use Exporter;

our @ISA= qw( Exporter );
our @EXPORT = qw( _align_and_clashcheck _translateRotate);


sub _align_and_clashcheck {

  my $ref_refatoms_x=shift;
  my $ref_refatoms_y=shift;
  my $ref_refatoms_z=shift;

  my $ref_cmpatoms_x=shift;
  my $ref_cmpatoms_y=shift;
  my $ref_cmpatoms_z=shift;

  my $ref_candatoms_x=shift;
  my $ref_candatoms_y=shift;
  my $ref_candatoms_z=shift;
  my $ref_candatoms_r=shift;

  my $ref_recatoms_x=shift;
  my $ref_recatoms_y=shift;
  my $ref_recatoms_z=shift;
  my $ref_recatoms_r=shift;

  my $ref_outcandatoms_x=shift;
  my $ref_outcandatoms_y=shift;
  my $ref_outcandatoms_z=shift;

  # compute transform from ref/cmp position and apply to inpdbatoms
  my $bStatus = _translateRotate($ref_refatoms_x,$ref_refatoms_y,$ref_refatoms_z,  $ref_cmpatoms_x,$ref_cmpatoms_y,$ref_cmpatoms_z, $ref_candatoms_x,$ref_candatoms_y,$ref_candatoms_z, $ref_outcandatoms_x,$ref_outcandatoms_y, $ref_outcandatoms_z);

  # =================================================================================
  # if it doesn't work, add a dummy atom

  if ($bStatus == -1 ) { 
#    print "Nonunitary Matrix\n";
    my @refAtomsX_tmp = @{$ref_refatoms_x};
    my @refAtomsY_tmp = @{$ref_refatoms_y};
    my @refAtomsZ_tmp = @{$ref_refatoms_z};
    my @cmpAtomsX_tmp = @{$ref_cmpatoms_x};
    my @cmpAtomsY_tmp = @{$ref_cmpatoms_y};
    my @cmpAtomsZ_tmp = @{$ref_cmpatoms_z};
    my $inc = 0.01; my $maxCyc = 5;
    my $ct = 1;
    
    while ( $bStatus == -1 && $ct < $maxCyc) {
	print "adding 4th pt with inc = $inc\n";
      $refAtomsX_tmp[3] = $refAtomsX_tmp[2] + $inc;
      $refAtomsY_tmp[3] = $refAtomsY_tmp[2] + $inc;
      $refAtomsZ_tmp[3] = $refAtomsZ_tmp[2] + $inc;
      $cmpAtomsX_tmp[3] = $cmpAtomsX_tmp[2] + $inc;
      $cmpAtomsY_tmp[3] = $cmpAtomsY_tmp[2] + $inc;
      $cmpAtomsZ_tmp[3] = $cmpAtomsZ_tmp[2] + $inc;
      $bStatus = _translateRotate(\@refAtomsX_tmp,\@refAtomsY_tmp,\@refAtomsZ_tmp,\@cmpAtomsX_tmp,\@cmpAtomsY_tmp,\@cmpAtomsZ_tmp, $ref_candatoms_x,$ref_candatoms_y,$ref_candatoms_z, $ref_outcandatoms_x,$ref_outcandatoms_y, $ref_outcandatoms_z);
      $inc *= 1.5;
    }
    ($bStatus != -1) || die "Die: Cannot recover from non-unitary matrix\n"; 
  }
  # =================================================================================



# clashcheck
  my $bClash;
  $bClash = _checkPQRsInteract ( $ref_outcandatoms_x,$ref_outcandatoms_y, $ref_outcandatoms_z,$ref_candatoms_r,  $ref_recatoms_x,$ref_recatoms_y,$ref_recatoms_z,$ref_recatoms_r);
  
  return $bClash;
}

# check extent of overlap between 2 pqr files 
sub _checkPQRsInteract{
  my $refx1 = $_[0];
  my $refy1 = $_[1];
  my $refz1 = $_[2];
  my $refr1 = $_[3];

  my $refx2 = $_[4];
  my $refy2 = $_[5];
  my $refz2 = $_[6];
  my $refr2 = $_[7];

# if at least one overlap, return true
  my $natom1 = scalar(@{$refx1});
  my $natom2 = scalar(@{$refx2});
    for(my $i=0;$i<$natom1;$i++) {
      for(my $j=0;$j<$natom2;$j++) {
	my $dx = $refx1->[$i]-$refx2->[$j];
	my $dy = $refy1->[$i]-$refy2->[$j];
	my $dz = $refz1->[$i]-$refz2->[$j];
	my $radsum = $refr1->[$i] + $refr2->[$j];
	my $distsq = ($dx*$dx + $dy*$dy + $dz*$dz);
	
	if ( $radsum*$radsum - $distsq > 0 ) {
	  my $diff= $radsum*$radsum - $distsq;
	  print "clash $refx1->[$i] $refy1->[$i] $refz1->[$i] :  $refx2->[$j] $refy2->[$j] $refz2->[$j]: $distsq $radsum $diff\n";
	  return 1;
	}
      }
    }

  return 0;
}

# ==============================================
# takes 2 xyzs for fitting (generate rotation matrix and translation vector) and
# apply rotation+translation to full input pdbfile
# !!!!! assume ref and fit xyzs  have correponding atoms in order !!!
# modified from mmtsb Analyze.pm
sub _translateRotate {
  my $ref_refatoms_x=shift;
  my $ref_refatoms_y=shift;
  my $ref_refatoms_z=shift;

  my $ref_cmpatoms_x=shift;
  my $ref_cmpatoms_y=shift;
  my $ref_cmpatoms_z=shift;

  my $ref_inpdbatoms_x=shift;
  my $ref_inpdbatoms_y=shift;
  my $ref_inpdbatoms_z=shift;

  my $ref_outpdbatoms_x=shift;
  my $ref_outpdbatoms_y=shift;
  my $ref_outpdbatoms_z=shift;

#  for (my $i=0;$i<scalar( @{$ref_refatoms_x});$i++) {     print $ref_refatoms_x->[$i]," ",$ref_refatoms_y->[$i]," ",$ref_refatoms_z->[$i],"\n";    }
#  for (my $i=0;$i<scalar( @{$ref_cmpatoms_x});$i++ ){     print $ref_cmpatoms_x->[$i]," ",$ref_cmpatoms_y->[$i]," ",$ref_cmpatoms_z->[$i],"\n";    }
#  for (my $i=0;$i<1;$i++) {     print $ref_refatoms_x->[$i]," ",$ref_refatoms_y->[$i]," ",$ref_refatoms_z->[$i],"\n";    }
#  for (my $i=0;$i<1;$i++ ){     print $ref_cmpatoms_x->[$i]," ",$ref_cmpatoms_y->[$i]," ",$ref_cmpatoms_z->[$i],"\n";    }
  
  # compute the geometric centers of refatoms and cmpatoms
  my $nref=0;
  my ($cxref,$cyref,$czref);
  my ($cxcmp,$cycmp,$czcmp);
  $cxref=$cyref=$czref=$cxcmp=$cycmp=$czcmp=0.0;

  for (my $ik=0; $ik<scalar( @{$ref_refatoms_x} ); $ik++) {
    
    $cxref+=$ref_refatoms_x->[$ik];
    $cyref+=$ref_refatoms_y->[$ik];
    $czref+=$ref_refatoms_z->[$ik];
    
    $cxcmp+=$ref_cmpatoms_x->[$ik];
    $cycmp+=$ref_cmpatoms_y->[$ik];
    $czcmp+=$ref_cmpatoms_z->[$ik];
    $nref++;
  }

#  print "nref = $nref\n";

  $cxref/=$nref;
  $cyref/=$nref;
  $czref/=$nref;
  
  $cxcmp/=$nref;
  $cycmp/=$nref;
  $czcmp/=$nref;

# rotation matrix
  my $r=();
  for (my $i=1; $i<=3; $i++) {
    $r->[$i]=();
    for (my $j=1; $j<=3; $j++) {
      $r->[$i]->[$j]=0.0;
    }
  }

  # translate coordinates to respective geometric centers
  for (my $ik=0; $ik<$nref; $ik++) {
      my $cx=$ref_cmpatoms_x->[$ik] - $cxcmp;
      my $cy=$ref_cmpatoms_y->[$ik] - $cycmp;
      my $cz=$ref_cmpatoms_z->[$ik] - $czcmp;
      my $rx=$ref_refatoms_x->[$ik] - $cxref;
      my $ry=$ref_refatoms_y->[$ik] - $cyref;
      my $rz=$ref_refatoms_z->[$ik] - $czref;

      $r->[1]->[1]+=$cx*$rx;
      $r->[2]->[1]+=$cx*$ry;
      $r->[3]->[1]+=$cx*$rz;

      $r->[1]->[2]+=$cy*$rx;
      $r->[2]->[2]+=$cy*$ry;
      $r->[3]->[2]+=$cy*$rz;

      $r->[1]->[3]+=$cz*$rx;
      $r->[2]->[3]+=$cz*$ry;
      $r->[3]->[3]+=$cz*$rz;
  }


#  print "$r->[1]->[1] $r->[1]->[2] $r->[1]->[3]\n";
#  print "$r->[2]->[1] $r->[2]->[2] $r->[2]->[3]\n";
#  print "$r->[3]->[3] $r->[3]->[2] $r->[3]->[3]\n";

  my $u=_frotu($r);
  if ( $u == -1 ) { return -1; }

#  print "$u->[1]->[1] $u->[1]->[2] $u->[1]->[3]\n";
#  print "$u->[2]->[1] $u->[2]->[2] $u->[2]->[3]\n";
#  print "$u->[3]->[3] $u->[3]->[2] $u->[3]->[3]\n";

# ===============================================
# rotate and translate the full inpdb coordinates
# ===============================================
  for (my $ik=0; $ik<scalar(@{$ref_inpdbatoms_x}); $ik++) {
    my $x=$ref_inpdbatoms_x->[$ik] - $cxcmp;
    my $y=$ref_inpdbatoms_y->[$ik] - $cycmp;
    my $z=$ref_inpdbatoms_z->[$ik] - $czcmp;
    
    my $tx=$u->[1]->[1]*$x+$u->[1]->[2]*$y+$u->[1]->[3]*$z;
    my $ty=$u->[2]->[1]*$x+$u->[2]->[2]*$y+$u->[2]->[3]*$z;
    my $tz=$u->[3]->[1]*$x+$u->[3]->[2]*$y+$u->[3]->[3]*$z;
    
    $ref_outpdbatoms_x->[$ik]=$tx+$cxref;
    $ref_outpdbatoms_y->[$ik]=$ty+$cyref;
    $ref_outpdbatoms_z->[$ik]=$tz+$czref;
  }  

  return 0;

}    


sub _frotu {
  my $r=shift;

  my ($i,$j,$k);

  my $det=0.0;
  for ($i=1; $i<=3; $i++) {
    my $i1=$i+1;
    $i1-=3 if ($i1>3);
    my $i2=$i+2;
    $i2-=3 if ($i2>3);
    $det+=$r->[$i]->[1]*($r->[$i1]->[2]*$r->[$i2]->[3]-$r->[$i2]->[2]*$r->[$i1]->[3]);
  }

  my $ipt=0;
  my $w;

  for ($i=1; $i<=3; $i++) {
    for ($j=$i; $j<=3; $j++) {
      $ipt++;
      $w->[$ipt]=0.0;
      for ($k=1; $k<=3; $k++) {
	$w->[$ipt]+=$r->[$j]->[$k]*$r->[$i]->[$k];
      }
    }
  }

  my $trace=$w->[1]+$w->[4]+$w->[6];

  my $u=();

  if ($trace<3.0E-6) {
    for ($i=1; $i<=3; $i++) {
      $u->[$i]=();
      for ($j=1; $j<=3; $j++) {
	$u->[$i]->[$j]=0.0;
      }
      $u->[$i]->[$i]=1.0;
    }
    return $u;
  }

  my ($vec,$ev)=_diagq(3,3,$w);

  my $a=();
  $ipt=1;
  for ($i=1; $i<=3; $i++) {
    for ($j=1; $j<=3; $j++) {
      $a->[$j]->[$i]=$vec->[$ipt];
      $ipt++;
    }
  }

  for ($i=1; $i<=3; $i++) {
    $ev->[$i]=sqrt(abs($ev->[$i]));
    $ev->[$i]=1.0E-6 if ($ev->[$i]<1.0E-6);
  }

  $ev->[1]=-$ev->[1] if ($det<0.0);

  my $b=();
  for ($j=1; $j<=3; $j++) {
    $b->[$j]=();
  }

  for ($j=1; $j<=3; $j++) {
    my $evs=$ev->[$j];
    for ($i=1; $i<=3; $i++) {
      $b->[$i]->[$j]=0.0;
      for ($k=1; $k<=3; $k++) {
	$b->[$i]->[$j]+=$r->[$k]->[$i]*$a->[$k]->[$j]/$evs;
      }
    }
  }

  $det=0.0;
  for ($i=1; $i<=3; $i++) {
    my $i1=$i+1;
    $i1-=3 if ($i1>3);
    my $i2=$i+2;
    $i2-=3 if ($i2>3);
    $det+=$a->[$i]->[1]*($a->[$i1]->[2]*$a->[$i2]->[3]-$a->[$i2]->[2]*$a->[$i1]->[3]);
  }

  for ($j=1; $j<=3; $j++) {
    if (abs($ev->[$j]) <= 1.0E-6) {
      my $jp=$j+1;
      my $jq=$j+2;
      $jp-=3 if ($jp>3);
      $jq-=3 if ($jq>3);
      for ($k=1; $k<=3; $k++) {
	my $kp=$k+1;
	my $kq=$k+2;
	$kp-=3 if ($kp>3);
	$kq-=3 if ($kq>3);
	$b->[$k]->[$j]=$b->[$kp]->[$jp]*$b->[$kq]->[$jq]-$b->[$kp]->[$jq]*$b->[$kq]->[$jp];
	$b->[$k]->[$j]=-$b->[$k]->[$j] if ($det<0.0);
      }
    }

    my $c=0.0;
    for ($k=1; $k<=3; $k++) {
      $c+=$b->[$k]->[$j]*$b->[$k]->[$j];
    }

    $c=($c>1.0E-10)?1.0/sqrt($c):0.0;

    for ($k=1; $k<=3; $k++) {
      $b->[$k]->[$j]*=$c;
    }
  }

  for ($j=1; $j<=3; $j++) {
    $u->[$i]=();
    for ($i=1; $i<=3; $i++) {
      $u->[$i]->[$j]=0.0;
      for ($k=1; $k<=3; $k++) {
	$u->[$i]->[$j]+=$a->[$i]->[$k]*$b->[$j]->[$k];
      }
    }
  }

  for ($j=1; $j<=3; $j++) {
    my $c=0.0;
    for ($k=1; $k<=3; $k++) {
      $c+=$u->[$k]->[$j]*$u->[$k]->[$j];
    }
    $c=($c>1.0E-10)?1.0/sqrt($c):0.0;

    for ($k=1; $k<=3; $k++) {
      $u->[$k]->[$j]*=$c;
    }
  }

  $det=0.0;
  
  for ($i=1; $i<=3; $i++) {
    my $i1=$i+1;
    $i1-=3 if ($i1>3);
    my $i2=$i+2;
    $i2-=3 if ($i2>3);
    $det+=$u->[$i]->[1]*($u->[$i1]->[2]*$u->[$i2]->[3]-$u->[$i2]->[2]*$u->[$i1]->[3]);
  }
 
# ORIGINAL 
#  printf STDERR "non-unitary rotation matrix, determinant: %f\n",$det
#    if (abs($det-1.0)>1.0E-4);
# ENGHUI
  if (abs($det-1.0)>1.0E-4) {return -1; }
     
  return $u;
}

sub _diagq {
  my $nx=shift;
  my $nfrqx=shift;
  my $dd=shift;

  my $vec=();
  my $ev=();

  my @a;
  my @b;
  my @p;
  my @w;
  my @ta;
  my @tb;
  my @y;
  
  my $nadd=0;

  my $eta=2.22045E-16;
  my $theta=4.4923E+307;

  my $n=$nx;
  my $nev=$nfrqx;
  my $nevadd=$nev+$nadd; 

  my $del1=$eta/100.0;
  my $delta=$eta*$eta*100.0;
  my $small=$eta*$eta/100.0;
  my $delbig=$theta*$delta/1000.0;
  my $theta1=1000.0/$theta;
  my $toler=100.0*$eta;
  my $rpower=8388608.0;
  my $rpow1=$rpower*0.50;
  my $rand1=$rpower-3.0;
  my $dunity=1.0;

  my $factor=0.0;
  my $ntot=int(($n*($n+1))/2);

  my ($i,$j,$k,$l,$m);

  for ($i=1; $i<=$ntot; $i++) {
    $factor=abs($dd->[$i]) if ($factor<abs($dd->[$i]));
  }

  if ($factor<$theta1) {
    printf STDERR "zero matrix passed in _diagq\n";

    for ($i=1; $i<=$nev; $i++) {
      $ev->[$i]=0.0;
      my $ipt=($i-1)*$n;
      for ($j=1; $j<$n; $j++) {
	$ipt++;
	$vec->[$ipt]=($i+$nadd == $j)?1.0:0.0;
      }
    }
    return ($vec,$ev);
  }

  my $ij=0;
  my $anorm=0.0;
  
  for ($i=1; $i<=$n; $i++) {
    for ($j=$i; $j<=$n; $j++) {
      $ij++;
      my $u=($dd->[$ij]/$factor)*($dd->[$ij]/$factor);
      $u*=0.5 if ($i == $j);
      $anorm+=$u;
    }
  }
  $anorm=sqrt($anorm+$anorm)*$factor;
  my $anormr=$dunity/$anorm;

  for ($i=1; $i<=$ntot; $i++) {
    $dd->[$i]*=$anormr;
  }

  my $nn=$n-1;
  my $mi=0;
  my $mi1=$n-1;

  for ($i=1; $i<=$nn; $i++) {
    my $sum1=0.0;
    $b[$i]=0.0;
    my $ji=$i+1;
    my $ipt=$mi+$i;
    $a[$i]=$dd->[$ipt];
    $ipt++;
    my $bx=$dd->[$ipt];
    my $ji2=$ji+1;
    if ($ji == $n) {
      $b[$i]=$bx;
      $dd->[$mi+$ji]=0.0;
      $mi+=$mi1;
      $mi1--;
    } else {
      for ($j=$ji2; $j<=$n; $j++) {
	$ipt++;
	$sum1+=$dd->[$ipt]*$dd->[$ipt];
      }
      
      if ($sum1<=$small) {
	$b[$i]=$bx;
	$dd->[$mi+$ji]=0.0;
	$mi+=$mi1;
	$mi1--;
      } else {
	my $s=sqrt($sum1+$bx*$bx);
	my $sgn=($bx>=0.0)?abs($dunity):-abs($dunity);
	my $temp=abs($bx);
	$w[$ji]=sqrt(0.5*($dunity+($temp/$s)));
	$ipt=$mi+$ji;
	$dd->[$ipt]=$w[$ji];
	my $ii=$i+2;
	if ($ii<=$n) {
	  $temp=$sgn/(2.0*$w[$ji]*$s);
	  for ($j=$ii; $j<=$n; $j++) {
	    $ipt++;
	    $w[$j]=$temp*$dd->[$ipt];
	    $dd->[$ipt]=$w[$j];
	  }
	}

	$b[$i]=-$sgn*$s;
	
	for ($j=$ji; $j<=$n; $j++) {
	  $p[$j]=0.0;
	}

	my $ml=$mi+$mi1;
	my $ml1=$mi1-1;


	for ($l=$ji; $l<=$n; $l++) {
	  $ipt=$ml+$l;
	  for ($m=$l; $m<=$n; $m++) {
	    $bx=$dd->[$ipt];
	    $p[$l]+=$bx*$w[$m];
	    $p[$m]+=$bx*$w[$l] if ($l!=$m);
	    $ipt++;
	  }
	  $ml+=$ml1;
	  $ml1--;
	}
	my $xkap=0.0;
	
	for ($k=$ji; $k<=$n; $k++) {
	  $xkap+=$w[$k]*$p[$k];
	}

	for ($l=$ji; $l<=$n; $l++) {
	  $p[$l]-=$xkap*$w[$l];
	}
	
	my $mj=$mi+$mi1;
	my $mj1=$mi1-1;

	for ($j=$ji; $j<=$n; $j++) {
	  for ($k=$j; $k<=$n; $k++) {
	    my $expr=($p[$j]*$w[$k])+($p[$k]*$w[$j]);
	    $dd->[$mj+$k]-=$expr+$expr;
	  }
	  $mj+=$mj1;
	  $mj1--;
	}
	
	$mi+=$mi1;
	$mi1--;
      }
    }
  }

  $a[$n]=$dd->[$mi+$n];
  $b[$n]=0.0;

  my $alimit=1.0;
  for ($i=1; $i<=$n; $i++) {
    $w[$i]=$b[$i];
    $b[$i]*=$b[$i];
  }
  
  for ($i=1; $i<=$nevadd; $i++) {
    $ev->[$i]=$alimit;
  }
  my $rootl=-$alimit;

  for ($i=1; $i<=$nevadd; $i++) {  
    my $rootx=$alimit;
    for ($j=$i; $j<=$nevadd; $j++) {
      $rootx=$ev->[$j] if ($ev->[$j]<$rootx);
    }
    $ev->[$i]=$rootx;

    my $trial=($rootl+$ev->[$i])*0.5;

    while(abs($trial-$rootl) > 1.0E-15 && abs($trial-$ev->[$i]) > 1.0E-15) {
      my $nomtch=$n;
      $j=1;
      
      do {
	my $f0=$a[$j]-$trial;
	
	while ($j<=$n && abs($f0)>=$theta1) {
	  $nomtch-- if ($f0>=0.0);
	  $j++;
	  $f0=$a[$j]-$trial-$b[$j-1]/$f0;
	}  
	if ($j<=$n) {
	  $j+=2;
	  $nomtch--;
	}
      }	while ($j<=$n);

      if ($nomtch>=$i) {
	$ev->[$i]=$trial;
	my $nom=($nevadd<=$nomtch)?$nevadd:$nomtch;
	$ev->[$nom]=$trial;
      } else {
	$rootl=$trial;
      }

      $trial=($rootl+$ev->[$i])*0.5;      
    }
  }

  for ($i=1; $i<=$nev; $i++) {
    $ev->[$i]=$ev->[$i+$nadd];
  }

  my $ia=0;
  
  for ($i=1; $i<=$nev; $i++) {
    my $aroot=$ev->[$i];
    for ($j=1; $j<=$n; $j++) {
      $y[$j]=1.0;
    }
    $ia=-1 if ($i==1 || abs($ev->[$i-1]-$aroot)>=$toler);
    $ia++;

    my $elim1=$a[1]-$aroot;
    my $elim2=$w[1];

    for ($j=1; $j<=$nn; $j++) {
      my $temp;
      if (abs($elim1)<=abs($w[$j])) {
	$ta[$j]=$w[$j];
	$tb[$j]=$a[$j+1]-$aroot;
	$p[$j]=$w[$j+1];
	$temp=(abs($w[$j])>$theta1)?$elim1/$w[$j]:1.0;
	$elim1=$elim2-$temp*$tb[$j];
	$elim2=-$temp*$w[$j+1];
      } else {
	$ta[$j]=$elim1;
	$tb[$j]=$elim2;
	$p[$j]=0.0;
	$temp=$w[$j]/$elim1;
	$elim1=$a[$j+1]-$aroot-$temp*$elim2;
	$elim2=$w[$j+1];
      }
      $b[$j]=$temp;
    }

    $ta[$n]=$elim1;
    $tb[$n]=0.0;
    $p[$n]=0.0;
    $p[$nn]=0.0;
    my $iter=1;
    
    if ($ia!=0) { 
      for ($j=1; $j<=$n; $j++) {
	my $rand1=(4099.0*$rand1 % $rpower);
	$y[$j]=$rand1/$rpow1-1.0;
      }
    }

    do {
      $l=$n+1;
      
      for ($j=1; $j<=$n; $j++) {
	$l--;
	do {
	  if (($n-$l-1)<0) {
	    $elim1=$y[$l];
	  } elsif (($n-$l-1)==0) {
	    $elim1=$y[$l]-$y[$l+1]*$tb[$l];
	  } else {
	    $elim1=$y[$l]-$y[$l+1]*$tb[$l]-$y[$l+2]*$p[$l];
	  }
	  if ($elim1>$delbig || $elim1<-$delbig) {
	    for ($k=1; $k<=$n; $k++) {
	      $y[$k]=$y[$k]/$delbig;
	    }
	  } 	  
	} while ($elim1>$delbig || $elim1<-$delbig);
	
	my $temp=$ta[$l];
	$temp=$delta if (abs($temp)<$delta);
	$y[$l]=$elim1/$temp;
      }
      
      if ($iter==1) {
	$elim1=$y[1];
	for ($j=1; $j<=$nn; $j++) {
	  if (abs($ta[$j]-$w[$j])<1E-15) {
	    $y[$j]=$y[$j+1];
	    $elim1=$elim1-$y[$j+1]*$b[$j];
	  } else {
	    $y[$j]=$elim1;
	    $elim1=$y[$j+1]-$elim1*$b[$j];
	  }
	}
	$y[$n]=$elim1;
      } 
      $iter++;
    } while ($iter<=2);
    
    my $ipt;
    if ($ia != 0) {
      for (my $j1=1; $j1<=$ia; $j1++) {
	$k=$i-$j1;
	my $temp=0.0;
	$ipt=($k-1)*$n;
	for ($j=1; $j<=$n; $j++) {
	  $ipt++;
	  $temp+=$y[$j]*$vec->[$ipt];
	}
	$ipt=($k-1)*$n;
	for ($j=1; $j<=$n; $j++) {
	  $ipt++;
	  $y[$j]-=$temp*$vec->[$ipt];
	}
      }
    }
    
    $elim1=0.0;
    for ($j=1; $j<=$n; $j++) {
      $elim1=abs($y[$j]) if (abs($y[$j])>$elim1);
    }
    my $temp=0.0;
    for ($j=1; $j<=$n; $j++) {
      $elim2=$y[$j]/$elim1;
      $temp+=$elim2*$elim2;
    }
    $temp=$dunity/(sqrt($temp)*$elim1);
    for ($j=1; $j<=$n; $j++) {
      $y[$j]*=$temp;
      $y[$j]=0.0 if (abs($y[$j])<$del1);
    }
    $ipt=($i-1)*$n;
    for ($j=1; $j<=$n; $j++) {
      $ipt++;
      $vec->[$ipt]=$y[$j];
    }
  }

  my $ipt;
  my $kk;
  for ($i=1; $i<=$nev; $i++) {
    $ipt=($i-1)*$n;
    for ($j=1; $j<=$n; $j++) {
      $ipt++;
      $y[$j]=$vec->[$ipt];
    }
    $l=$n-2;
    my $mk=($n*($n-1))/2-3;
    my $mk1=3;
    
    my $t;
    for ($j=1; $j<=$l; $j++) {
      $t=0.0;
      $k=$n-$j-1;
      $m=$k+1;
      for ($kk=$m; $kk<=$n; $kk++) {
	$t+=$dd->[$mk+$kk]*$y[$kk];
      }
      for ($kk=$m; $kk<=$n; $kk++) {
	my $epr=$t*$dd->[$mk+$kk];
	$y[$kk]-=$epr+$epr;
      }
      $mk-=$mk1;
      $mk1++;
    }

    $t=0.0;
    for ($j=1; $j<=$n; $j++) {
      $t+=$y[$j]*$y[$j];
    }
    
    my $xnorm=sqrt($t);
    my $xnorm1=$dunity/$xnorm;
    for ($j=1; $j<=$n; $j++) {
      $y[$j]*=$xnorm1;
    }
    
    $ipt=($i-1)*$n;
    for ($j=1; $j<=$n; $j++) {
      $ipt++;
      $vec->[$ipt]=$y[$j];
    }
  }

  for ($i=1; $i<=$n; $i++) {
    $ev->[$i]=$ev->[$i]*$anorm;
  }

  return ($vec,$ev);
}
1;
