#!/usr/bin/perl

=cut
$query_string = $ENV{'QUERY_STRING'};
print "PICS query: $query_string\n\n";
foreach $item (split(/&/, $query_string)) {
  if ($item =~ /command1=(\S+)/) {
    $ARGV[0] = $1;
  }
  if ($item =~ /command2=(\S+)/) {
    $ARGV[1] = $1;
  }
  if ($item =~ /command3=(\S+)/) {
    $ARGV[2] = $1;
  }
}
=cut

$num_args = 3;
if ($num_args < 2 ) {
  print "Usage: pics.pl <SNP> <pval> <ld>\n";
  exit;
}

$idx = 0;
$ratiosum = 0;
$myrs = $ARGV[0];
#$myrs =~ s/_/:/g;
$myrs =~ s/RS/rs/;
$pval = $ARGV[1];

if (!(($myrs =~ /rs\d+/) || ($myrs =~ /chr\d+:\d+/) || ($myrs =~ /MERGED.*/))) {
  print "error:  invalid SNP\n";
  die;
}
if (!($pval =~ /\d+/)) {
  print "error:  invalid pvalue\n";
  die;
}
if ($pval<0) {
  print "error:  negative pvalue\n";
  die;
}
if ($pval<1) {
  if ($pval==0) {
    $pval = 323;
  }
  else {
    $pval = -log($pval)/log(10);
  }
}
$indexpval = $pval;

$ldpath=$ARGV[2];

$idx++;
$linedata{$idx} = "$myrs\t$myrs\t1.0000\t1.0000\tN,N\t";
$ratiodata{$idx} = 1.0000;
$ratiosum = $ratiosum + $ratiodata{$idx};

print "Index_SNP\tLinked_SNP\tDprime\tRsquare\tPhase\tPICS_probability\n";
$cachefound = 0;
$datamissing = `grep -w $myrs $ldpath/cachemissing`;
if ($datamissing) {
  #print "$myrs\t$myrs\t1.0000\t1.0000\tN,N\t1.0000\n";
  print STDERR "$myrs not in LD database (cached)\n";
  exit;
}
$datacache = `grep -w $myrs $ldpath/cacheindex`;
if ($datacache) {
  $datald = `grep -w $myrs $ldpath/cache | tr ' ' '\t' | sort -nrk4`;
  $cachefound = 1;
}

$dbfound = 0;

if ($cachefound>0) {
  foreach $item (split(/\n/,$datald)) {
     if (!($alreadyseen{$item})) {
      if ($item =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S)\,(\S)/) {
	$rs1 = $1;
        $rs2 = $2;
	$dprime = $3;
        $r2 = $4;
        $phase1 = $5;
        $phase2 = $6;
        $pval1 = $indexpval+(1-r2)*.3010;
        $pval2 = $r2*$indexpval+(1-r2)*.3010;
        $stddev = sqrt(1-($r2**3.2))*($pval1**0.5)/2;
        if ($pval1==$pval2) {
          $prob1 = 0.5;
          $prob2 = 0.5;
        }
        else {
#         $prob1 = 1-Math::CDF::pnorm(($pval1-$pval2)/$stddev);
          $d = ($pval1-$pval2)/$stddev;
          $A1 = 0.31938153;
          $A2 = -0.356563782;
          $A3 = 1.781477937;
          $A4 = -1.821255978;
          $A5 = 1.330274429;
          $RSQRT2PI = 0.39894228040143267793994605993438;
          $K = 1.0 / (1.0 + 0.2316419 * abs($d));
          $cnd = $RSQRT2PI * exp(- 0.5 * $d * $d) * ($K * ($A1 + $K * ($A2 + $K * ($A3 + $K * ($A4 + $K * $A5)))));
          if ($d > 0) {
            $cnd = 1.0 - $cnd;
          }
          $prob1 = 1-$cnd;
          $prob2 = 1-$prob1;
        }
        $ratio = $prob1/$prob2;
	if ($rs1 eq $myrs) {
	  $phase = $phase1 . "," . $phase2;
          $idx++;
          $linedata{$idx} = "$rs1\t$rs2\t$dprime\t$r2\t$phase\t";
          $ratiodata{$idx} = $ratio;
          $ratiosum = $ratiosum + $ratio;
        }
        elsif ($rs2 eq $myrs) {
	  $phase = $phase2 . "," . $phase1;
          $idx++;
          $linedata{$idx} = "$rs2\t$rs1\t$dprime\t$r2\t$phase\t";
          $ratiodata{$idx} = $ratio;
          $ratiosum = $ratiosum + $ratio;
        }
      }
      $alreadyseen{$item} = 1;
    }
  }
}
else {
  print STDERR "$myrs not in LD database - $dbfound\n";
}

for ($i=1; $i<=$idx; $i++) {
  $ratio = int($ratiodata{$i}/$ratiosum*10000)/10000;
  print $linedata{$i} . $ratio . "\t" . $pval . "\n";
}
