#!/usr/bin/env perl

# some basic functions to load files (pdbs, etc)
# 
package loadCoordFiles;
use strict;
use warnings;
use Exporter;

our @ISA= qw( Exporter );
our @EXPORT = qw(_loadPQR_xyzr  _loadPDB_xyz _loadPDB_xyz_atomRes);

sub _loadPQR_xyzr {

  my $pdbfile = $_[0];
  my $refx = $_[1];
  my $refy = $_[2];
  my $refz = $_[3];
  my $refr = $_[4];

  open(FILE, $pdbfile) or die ("Can't open $pdbfile\n");
  my $n=0;
  while(my $line=<FILE>){
    if($line !~/^ATOM/) {next; } 
    $refx->[$n] = substr($line,30,8); 
    $refy->[$n] = substr($line,38,8); 
    $refz->[$n] = substr($line,46,8);     
    $refr->[$n] = substr($line,62,7); 
    $n++;
  }
  close FILE;
}

sub _loadPDB_xyz {
  my $pdbfile = $_[0];
  my $refx = $_[1];
  my $refy = $_[2];
  my $refz = $_[3];

  open(FILE, $pdbfile) or die ("Can't open $pdbfile\n");
  my $n=0;
  while(my $line=<FILE>){
    if($line !~/^ATOM/) {next; } 
    $refx->[$n] = substr($line,30,8); 
    $refy->[$n] = substr($line,38,8); 
    $refz->[$n] = substr($line,46,8);     
    $n++;
  }
  close FILE;
}

sub _loadPDB_xyz_atomRes {
  my $pdbfile = shift;
  my $refx = shift;
  my $refy = shift;
  my $refz = shift;
  my $refAtomname = shift;
  my $refResno = shift; 
  my $refResidue = shift;

 open(FILE, $pdbfile) or die ("Can't open $pdbfile\n");
  my $n=0;
  while(my $line=<FILE>){
    if($line !~/^ATOM/) {next; } 
    $refAtomname->[$n] = substr($line,12,4); 
    $refResidue->[$n]  = substr($line,17,3); 
    $refResno->[$n]    = substr($line,22,4); 
    $refx->[$n] = substr($line,30,8); 
    $refy->[$n] = substr($line,38,8); 
    $refz->[$n] = substr($line,46,8);     
    $n++;
  }
  close FILE;
}

1;
