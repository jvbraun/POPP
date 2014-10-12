package Finance::Options;

# $Id: Options.pm,v 1.16 2000/11/25 12:47:07 jvbraun Exp $
#
# Copyright (c) 1999 Jerome V. Braun. All rights reserved.
#
# This program is free software; you can redistribute it and/or
# modify it under the same terms as Perl itself.


=head1 NAME

Finance::Options - Options pricing and related calculations

=head1 SYNOPSIS

In general, for the examples, and in the code, I have made use of
the same variable names used in Haug (1998) to ease comparison.

  # usage
  use Finance::Options;

  # calculate a fair price using the Black-Scholes model
  $price = BlackScholes($call_put_flag, $S, $X, $T, $r, $v);

  # determine implied volatility using the bisection algorithm
  $volatility = ImpliedVolatilityB('BlackScholes', $price, $call_put_flag, $S, $X, $T, $r);

=head1 DESCRIPTION

Finance::Options is a set of routines for performing options pricing and
related calculations.

All of the routines here are drawn from Haug, E.G. (1998)
The Complete Guide to Option Pricing Formulas.  McGraw Hill.
This is a very readable and practical guide to implementing a lot
of option models, as well as related material.

=cut

#---------------------------------------------------------------------
# Start-up

use strict;
use vars qw($VERSION @ISA @EXPORT_OK);

require Exporter;

@ISA          = qw(Exporter);
@EXPORT_OK    = qw(&BlackScholes &GBlackScholes &BSAmericanApprox
                   &ImpliedVolatilityB &ImpliedVolatilityNR
                   &CRRBinomial &RollGeskeWhaley &BAWAmericanApprox
                   &sgn &CBND &CND &TrinomialTree &JumpDiffusion);

$VERSION = '0.03';   # also change version down below in the POD

=head1 Vanilla Options

=head2 European Options

=over 4

=item BlackScholes

Routine to implement the Black and Scholes (1973) option pricing formula.

  $price = BlackScholes($call_put_flag, $S, $X, $T, $r, $v);

C<$call_put_flag> is either 'c' or 'p' for a call or put respectively,
C<$S> is the underlying price, C<$X> is the strike price,
C<$T> is calendar days to maturity/calendar days per year,
C<$r> is the risk-free interest rate,
and C<$v> is the volatility in percent per annum.

=cut


sub BlackScholes {
  my ($call_put_flag, $S, $X, $T, $r, $v) = @_;

  # calculate some auxiliary values
  my $d1 = ( log($S/$X) + ($r+$v**2/2)*$T ) / ( $v * $T**0.5 );
  my $d2 = $d1 - $v * $T**0.5;

  if ($call_put_flag eq 'c') {
    return $S * &CND($d1) - $X * exp( -$r * $T ) * &CND($d2);
  }
  else {      # ($call_put_flag eq 'p')
    return  $X * exp( -$r * $T ) * &CND(-$d2) - $S * &CND(-$d1);
  }

}


=item GBlackScholes

Routine to implement the generalized Black-Scholes option pricing formula.

  $price = GBlackScholes($call_put_flag, $S, $X, $T, $r, $b, $v);

C<$call_put_flag> is either 'c' or 'p' for a call or put respectively,
C<$S> is the underlying price, C<$X> is the strike price,
C<$T> is calendar days to maturity/calendar days per year,
C<$r> is the risk-free interest rate, C<$b> is the cost-of-carry term,
and C<$v> is the volatility in percent per annum.

=cut

sub GBlackScholes {
  my ($call_put_flag, $S, $X, $T, $r, $b, $v) = @_;

  # Debug output
  #print "GBlackScholes got: $call_put_flag, $S, $X, $T, $r, $b, $v\n";

  # calculate some auxiliary values
  my $d1 = ( log($S/$X) + ($b+$v**2/2)*$T ) / ( $v * $T**(.5) );
  my $d2 = $d1 - $v * $T**(.5);

  if ($call_put_flag eq 'c') {
    return $S * exp( ($b - $r) * $T ) * &CND($d1) - $X * exp(-$r * $T) * &CND($d2);
  }
  else {          # ($call_put_flag eq 'p')
    return $X * exp(-$r * $T) * &CND(-$d2) - $S * exp( ($b - $r) * $T ) * &CND(-$d1);
  }
}


=item GDelta

Calculate Delta for the generalized Black-Scholes formula.
Inputs as for C<GBlackScholes>.

=cut

sub GDelta {
  my ($call_put_flag, $S, $X, $T, $r, $b, $v) = @_;

  my $d1 = (log($S/$X) + ($b + $v**2 / 2) * $T) / ($v * $T**0.5);

  if ($call_put_flag eq 'c') {
    return exp(($b-$r) * $T) *  &CND($d1);
  }
  else {
    return exp(($b-$r) * $T) * (&CND($d1) - 1);
  }
}


=item GGamma

Gamma for the generalized Black-Scholes formula.
Inputs as for C<GBlackScholes>.

=cut

sub GGamma {
  my ($S, $X, $T, $r, $b, $v) = @_;
  my $d1 = (log($S/$X) + ($b + $v**2 / 2) * $T) / ($v * $T**0.5);
  return exp(($b-$r) * $T) * &ND($d1) / ($S * $v * $T**0.5);
}


=item GTheta

Theta for the generalized Black-Scholes formula.
Inputs as for C<GBlackScholes>.

=cut

sub GTheta {
  my ($call_put_flag, $S, $X, $T, $r, $b, $v) = @_;

  my $d1 = (log($S/$X) + ($b + $v**2 / 2) * $T) / ($v * $T**0.5);
  my $d2 = $d1 - $v * $T**0.5;

  if ($call_put_flag eq 'c') {
    return -$S * exp(($b - $r) * $T) * &ND($d1) * $v / (2 * $T**0.5)
           - ($b - $r) * $S * exp(($b - $r) * $T) * &CND($d1)
           - $r * $X * exp(-$r * $T) * &CND($d2);
  }
  else {
    return -$S * exp(($b - $r) * $T) * &ND($d1) * $v / (2 * $T**0.5)
           + ($b - $r) * $S * exp(($b - $r) * $T) * &CND(-$d1)
           + $r * $X * exp(-$r * $T) * &CND(-$d2);
  }

}


=item Vega

Vega for the generalized Black-Scholes formula.
Inputs as for C<GBlackScholes>.

=cut

sub GVega {
  my ($S, $X, $T, $r, $b, $v) = @_;

  my $d1 = (log($S/$X) + ($b + $v**2 / 2) * $T) / ($v * $T**0.5);
  return $S * exp(($b - $r) * $T) * &ND($d1) * $T**0.5;
}

=item GRho

Rho for the generalized Black-Scholes formula.
Inputs as for C<GBlackScholes>.

=cut

sub GRho {
  my ($call_put_flag, $S, $X, $T, $r, $b, $v) = @_;

  my $d1 = (log($S/$X) + ($b + $v**2 / 2) * $T) / ($v * $T**0.5);
  my $d2 = $d1 - $v * $T**0.5;

  # call
  if ($call_put_flag eq 'c') {
    if (!$b==0) {
      return $T * $X * exp(-$r * $T) * &CND($d2);
    }
    else {
      return -$T * &GBlackScholes($call_put_flag, $S, $X, $T, $r, $b, $v);
    }
  }

  # put
  else {
    if (!$b==0) {
      return -$T * $X * exp(-$r * $T) * &CND(-$d2);
    }
    else {
      return -$T * &GBlackScholes($call_put_flag, $S, $X, $T, $r, $b, $v);
    }
  }

}

=item GCarry

Carry for the generalized Black-Scholes formula.
Inputs as for C<GBlackScholes>.

=cut

sub GCarry {
  my ($call_put_flag, $S, $X, $T, $r, $b, $v) = @_;

  my $d1 = (log($S/$X) + ($b + $v**2 / 2) * $T) / ($v * $T**0.5);

  if ($call_put_flag eq 'c') {
    return $T * $S * exp(($b - $r) * $T) * &CND($d1);
  }
  else {
    return -$T * $S * exp(($b - $r) * $T) * &CND(-$d1);
  }

}


=item French

French (1984) adjusted Black-Scholes model for trading day volatility.
C<$call_put_flag> is either 'c' or 'p' for a call or put respectively, 
C<$S> is the underlying price, C<$X> is the strike price,
C<$T> is calendar days to maturity/calendar days per year,
C<$t1> is trading days to maturity/trading days per year,
C<$r> is the risk-free interest rate, C<$b> is the cost-of-carry term,
and C<$v> is the volatility in percent per annum.

=cut

sub French {
  my ($call_put_flag, $S, $X, $T, $t1, $r, $b, $v) = @_;

  my $d1 = (log($S/$X) + $b * $T + $v**2 / 2 * $t1) / ($v * $t1**0.5);
  my $d2 = $d1 - $v * $t1**0.5;

  if ($call_put_flag eq 'c') {
    return $S * exp(($b - $r) * $T) * &CND($d1) - $X * exp(-$r * $T) * &CND($d2);
  }
  else {
    return $X * exp(-$r * $T) * &CND(-$d2) - $S * exp(($b - $r) * $T) * &CND(-$d1);
  }

}


=item JumpDiffusion

Merton (1976) jump diffusion model.
C<$call_put_flag> is either 'c' or 'p' for a call or put respectively,
C<$S> is the underlying price, C<$X> is the strike price,
C<$T> is calendar days to maturity/calendar days per year,
C<$r> is the risk-free interest rate, C<$b> is the cost-of-carry term,
C<$lambda> is the expected number of jumps per year,
C<$gamma> is percentage of total volatility explained by the jumps,
and C<$v> is the volatility in percent per annum.

=cut

sub JumpDiffusion {
  my ($call_put_flag, $S, $X, $T, $r, $lambda, $gamma, $v) = @_;

  my $delta = ($gamma * $v**2 / $lambda)**0.5;
  my $Z = ($v**2 - $lambda * $delta**2)**0.5;

  my $sum = 0;
  for (my $i=0; $i<=10; $i++) {
    my $vi = ($Z**2 + $delta**2 * ($i / $T))**0.5;
    $sum+= exp(-$lambda * $T) * ($lambda * $T)** $i / &factorial($i)
           * &GBlackScholes($call_put_flag, $S, $X, $T, $r, $r, $vi);
  }

  return $sum;

}


=item Merton73

Merton (1973) Options on stock indices model.
C<$call_put_flag> is either 'c' or 'p' for a call or put respectively,
C<$S> is the underlying price, C<$X> is the strike price,
C<$T> is calendar days to maturity/calendar days per year,
C<$r> is the risk-free interest rate, C<$q> is dividend yield,
and C<$v> is the volatility in percent per annum.

=cut

sub Merton73 {
  my ($call_put_flag, $S, $X, $T, $r, $q, $v) = @_;

  my $d1 = (log($S / $X) + ($r - $q + $v**2 / 2) * $T) / ($v * $T**0.5);
  my $d2 = $d1 - $v * $T**0.5;

  if ($call_put_flag eq 'c') {
    return $S * exp(-$q * $T) * &CND($d1) - $X * exp(-$r * $T) * &CND($d2);
  }
  else {
    return $X * exp(-$r * $T) * &CND(-$d2) - $S * exp(-$q * $T) * &CND(-$d1);
  }
}

=item Black76

Black (1977) Options on futures or forwards.

  $price = Black76($call_put_flag, $F, $X, $T, $r, $v);

C<$call_put_flag> is either 'c' or 'p' for a call or put respectively,
C<$F> is the price on the futures contract, C<$X> is the strike price,
C<$T> is calendar days to maturity/calendar days per year,
C<$r> is the risk-free interest rate,
and C<$v> is the volatility in percent per annum.

=cut

sub Black76 {
  my ($call_put_flag, $F, $X, $T, $r, $v) = @_;

  my $d1 = (log($F / $X) + ($v**2 / 2) * $T) / ($v * $T**0.5);
  my $d2 = $d1 - $v * $T**0.5;

  if ($call_put_flag eq 'c') {
    return exp(-$r * $T) * ($F * &CND($d1) - $X * &CND($d2));
  }
  else {
    return exp(-$r * $T) * ($X * &CND(-$d2) - $F * &CND(-$d1));
  }
}


=item MiltersonSchwarz

Miltersen and Schwartz (1997) model for options on commodity futures.

  $price = MiltersonSchwarz($call_put_flag, $Pt, $FT, $X, $t1, $T2,
                            $vS, $vE, $vf,
                            $rhoSe, $rhoSf, $rhoef, $Kappae, $Kappaf);

C<$call_put_flag> is either 'c' or 'p' for a call or put respectively,
C<$Pt> is a zero-coupon bond that expires on maturity,
C<$FT> is a futures price with time to expiration <$T2>,
C<$X> is the strike price,
C<$t1> is time to maturity of the option,
C<$vS>, C<$vE>, and C<$vf> are volatilities of spot commodity price,
future convenience yield, and forward interest rate respectively,
C<$rhoSe>, C<$rhoSf>, and C<$rhoef> are correlations between
spot commodity price and forward interest rate, between spot commodity
price and convenience yield, and between forward interest rate and
convenience yield respectively,
and C<$Kappae> and C<$Kappaf> are the speeds of mean reversion of
the forward interest rate and covenience yield respectively.

=back

=cut

sub MiltersonSchwarz {
  my ($call_put_flag, $Pt, $FT, $X, $t1, $T2,
      $vS, $vE, $vf, $rhoSe, $rhoSf, $rhoef, $Kappae, $Kappaf) = @_;

  my $vz =
    $vS**2 * $t1
    + 2*$vS * ( $vf * $rhoSf * 1/$Kappaf
                    * ($t1 - 1/$Kappaf * exp(-$Kappaf * $T2) * (exp($Kappaf * $t1)-1))
               -$vE * $rhoSe * 1/$Kappae
                    * ($t1 - 1/$Kappae * exp(-$Kappae * $T2) * (exp($Kappae * $t1)-1))
               )
    + $vE**2 * 1/$Kappae**2
          * ($t1
             + 1/(2*$Kappae) * exp(-2*$Kappae * $T2) * (exp(2*$Kappae * $t1)-1)
             - 2*1/$Kappae * exp(-$Kappae * $T2) * (exp($Kappae * $t1)-1)
            )
    + $vf**2 * 1/$Kappaf**2
          * ($t1
             + 1/(2*$Kappaf) * exp(-2*$Kappaf * $T2) * (exp(2*$Kappaf * $t1) - 1)
             - 2*1/$Kappaf * exp(-$Kappaf * $T2) * (exp($Kappaf * $t1) - 1))
    - 2 * $vE * $vf * $rhoef * 1/$Kappae * 1/$Kappaf
        * ($t1
           - 1/$Kappae * exp(-$Kappae * $T2) * (exp($Kappae * $t1) - 1)
           - 1/$Kappaf * exp(-$Kappaf * $T2) * (exp($Kappaf * $t1) - 1)
           + 1/($Kappae + $Kappaf)
               * exp(-($Kappae + $Kappaf) * $T2) * (exp(($Kappae + $Kappaf) * $t1) - 1)
           );

  my $vxz =
    $vf * 1/$Kappaf
      * ($vS * $rhoSf * ($t1 - 1/$Kappaf * (1 - exp(-$Kappaf * $t1)))
         + $vf * 1/$Kappaf
              * ($t1 - 1/$Kappaf * exp(-$Kappaf * $T2) * (exp($Kappaf * $t1) - 1)
                     - 1/$Kappaf * (1 - exp(-$Kappaf * $t1))
                     + 1/(2 * $Kappaf) * exp(-$Kappaf * $T2)
                                       * (exp($Kappaf * $t1) - exp(-$Kappaf * $t1))
                )
         - $vE * $rhoef * 1/$Kappae
              * ($t1 - 1/$Kappae * exp(-$Kappae * $T2) * (exp($Kappae * $t1) - 1)
                     - 1/$Kappaf * (1 - exp(-$Kappaf * $t1))
                     + 1/($Kappae + $Kappaf) * exp(-$Kappae * $T2)
                                             * (exp($Kappae * $t1) - exp(-$Kappaf * $t1))
                )
        );

  $vz = $vz**0.5;

  my $d1 = (log($FT / $X) - $vxz + $vz**2 / 2) / $vz;
  my $d2 = (log($FT / $X) - $vxz - $vz**2 / 2) / $vz;

   if ($call_put_flag eq 'c') {
     return $Pt * ($FT * exp(-$vxz) * &CND($d1) - $X * &CND($d2));
   }
  else {
     return $Pt * ($X * &CND(-$d2) - $FT * exp(-$vxz) * &CND(-$d1));
  }

}


=head2 American Options

=over 4

=item RollGeskeWhaley

American calls on stocks with known dividends.

  $price = RollGeskeWhaley($S, $X, $t1, $T2, $r, $d, $v);

C<$S> is the underlying price, C<$X> is the strike price,
C<$t1> is the time to dividend payout, C<$T2> is the time to option expiration,
C<$r> is the risk-free interest rate, C<$d> is the cash dividend,
and C<$v> is the volatility in percent per annum.


=cut

#  To do:  1.  Make the magic numbers more standard or parameterized.

sub RollGeskeWhaley {
  my ($S, $X, $t1, $T2, $r, $d, $v) = @_;

  # magic numbers
  my ($infinity, $epsilon) = (100000000, 0.00001);

  my $Sx = $S - $d * exp(-$r * $t1);

  # not optimal to exercise
  if ( $d <= $X * (1 - exp(-$r * ($T2-$t1))) ) {
    return &BlackScholes('c', $Sx, $X, $T2, $r, $v);
  }

  my $ci = &BlackScholes('c', $S, $X, $T2-$t1, $r, $v);
  my $HighS = $S;
    while ( ($ci - $HighS - $d + $X > 0) and ($HighS < $infinity) ) {
      $HighS = $HighS * 2;
      $ci = &BlackScholes('c', $HighS, $X, $T2 - $t1, $r, $v);
    }

  if ($HighS > $infinity) {
      return &BlackScholes('c', $Sx, $X, $T2, $r, $v);
  }

  my $LowS = 0;
  my $I = $HighS * 0.5;
  $ci = BlackScholes('c', $I, $X, $T2 - $t1, $r, $v);

  # Search algorithm to find the critical stock price I
  while (abs($ci - $I - $d + $X) > $epsilon and ($HighS - $LowS > $epsilon)) {
    if ($ci - $I - $d + $X < 0) {
      $HighS = $I;
    }
    else {
      $LowS = $I;
    }
    $I = ($HighS + $LowS) / 2;
    $ci = &BlackScholes('c', $I, $X, $T2 - $t1, $r, $v);
  }

  my $a1 = (log($Sx / $X) + ($r + $v**2 / 2) * $T2) / ($v * $T2**(0.5));
  my $a2 = $a1 - $v * $T2**(0.5);
  my $b1 = (log($Sx / $I) + ($r + $v**2 / 2) * $t1) / ($v * $t1**(0.5));
  my $b2 = $b1 - $v * $t1**(0.5);

  return $Sx * &CND($b1) + $Sx * &CBND($a1, -$b1, -($t1 / $T2)**(0.5))
         - $X * exp(-$r * $T2) * &CBND($a2, -$b2, -($t1 / $T2)**(0.5))
         - ($X - $d) * exp(-$r * $t1) * &CND($b2);

}


=item BAWAmericanApprox

Barone-Adesi and Whaley (1987) American approximation.

  $price = BAWAmericanApprox($call_put_flag, $S, $X, $T, $r, $b, $v);

C<$call_put_flag> is either 'c' or 'p' for a call or put respectively,
C<$S> is the underlying price, C<$X> is the strike price,
C<$T> is calendar days to maturity/calendar days per year,
C<$r> is the risk-free interest rate, C<$b> is the cost-of-carry term,
and C<$v> is the volatility in percent per annum.

=cut

# To do:  Add sanity checking code to the Newton-Raphson iterations.

sub BAWAmericanApprox {
  my ($call_put_flag, $S, $X, $T, $r, $b, $v) = @_;

  if ($call_put_flag eq 'c') {
    return BAWAmericanCallApprox($S, $X, $T, $r, $b, $v);
  }
  else {
    return  BAWAmericanPutApprox($S, $X, $T, $r, $b, $v);
  }
}

sub BAWAmericanCallApprox {
  my ($S, $X, $T, $r, $b, $v) = @_;

  if ($b>=$r) {
    return &GBlackScholes('c', $S, $X, $T, $r, $b, $v);
  }
  else {
    my $Sk = &Kc($X, $T, $r, $b, $v);
    my $n  = 2 * $b / $v**2;
    my $K  = 2 * $r / ($v**2 * (1 - exp(-$r * $T)));
    my $d1 = (log($Sk / $X) + ($b + $v**2 / 2) * $T) / ($v * $T**(0.5));
    my $Q2 = ( -($n - 1) + (($n - 1)**2 + 4 * $K)**(0.5) ) / 2;
    my $a2 = ($Sk / $Q2) * (1 - exp(($b - $r) * $T) * &CND($d1));

    if ($S<$Sk) {
      return GBlackScholes('c', $S, $X, $T, $r, $b, $v) + $a2 * ($S / $Sk)**$Q2;
    }
    else {
      return $S-$X;
    }
  }

}

# Newton Raphson algorithm to solve for the critical commodity price for a Call
sub Kc {
  my ($X, $T, $r, $b, $v) = @_;

  # Calculation of seed value, Si
  my $n   = 2 * $b / $v**2;
  my $m   = 2 * $r / $v**2;
  my $q2u = ( -($n - 1) + (($n - 1)**2 + 4 * $m)**(0.5) ) / 2;
  my $Su  = $X / (1 - 1 / $q2u);
  my $h2  = -($b * $T + 2 * $v * $T**(0.5)) * $X / ($Su - $X);
  my $Si  = $X + ($Su - $X) * (1 - exp($h2));

  my $K   = 2 * $r / ( $v**2 * (1 - exp(-$r * $T)) );
  my $d1  = ( log($Si/$X) + ($b + $v**2 / 2) * $T ) / ($v * $T**(0.5));
  my $Q2  = ( -($n - 1) + (($n - 1)**2 + 4 * $K)**(0.5) ) / 2;
  my $LHS = $Si - $X;
  my $RHS = &GBlackScholes('c', $Si, $X, $T, $r, $b, $v)
            + (1 - exp(($b - $r) * $T) * &CND($d1) ) * $Si / $Q2;
  my $bi  =        exp(($b - $r) * $T) * &CND($d1) * (1 - 1 / $Q2)
            + (1 - exp(($b - $r) * $T) *  &ND($d1) / ($v * $T**(0.5))) / $Q2;

########## DOUBLE-CHECK AS TO WHETHER OR NOT THE &ND SHOULD BE
########## &CND AS ORIGINALLY ABOVE.

  my $epsilon = 0.000001;

  # Newton Raphson algorithm for finding critical price Si
  while ( abs($LHS - $RHS)/$X > $epsilon ) {
    $Si = ($X + $RHS - $bi * $Si) / (1 - $bi);
    $d1 = ( log($Si/$X) + ($b + $v**2 / 2) * $T ) / ( $v * $T**(0.5) );
    $LHS = $Si - $X;
    $RHS = &GBlackScholes('c', $Si, $X, $T, $r, $b, $v)
           + (1 - exp(($b - $r) * $T) * &CND($d1)) * $Si / $Q2;
    $bi =         exp(($b - $r) * $T) * &CND($d1) * (1 - 1 / $Q2)
           + (1 - exp(($b - $r) * $T) *  &ND($d1) / ($v * $T**(0.5))) / $Q2;
  }

  return $Si;

}

sub BAWAmericanPutApprox {
  my ($S, $X, $T, $r, $b, $v) = @_;

  my $Sk = &Kp($X, $T, $r, $b, $v);
  my $n = 2 * $b / $v**2;
  my $K = 2 * $r / ($v**2 * (1 - exp(-$r * $T)));
  my $d1 = (log($Sk / $X) + ($b + $v**2 / 2) * $T) / ($v * $T**(0.5));
  my $Q1 = ( -($n - 1) - (($n - 1)**2 + 4 * $K)**(0.5) ) / 2;
  my $a1 = -($Sk / $Q1) * (1 - exp(($b - $r) * $T) * &CND(-$d1));

  if ($S>$Sk) {
    return &GBlackScholes('p', $S, $X, $T, $r, $b, $v) + $a1 * ($S / $Sk)**$Q1;
  }
  else {
    return $X-$S;
  }
}

# Newton Raphson algorithm to solve for the critical commodity price for a Put
sub Kp {
  my ($X, $T, $r, $b, $v) = @_;

  # Calculation of seed value, Si
  my $n = 2 * $b / $v**2;
  my $m = 2 * $r / $v**2;
  my $q1u = ( -($n - 1) - (($n - 1)**2 + 4 * $m)**(0.5) ) / 2;
  my $Su = $X / (1 - 1 / $q1u);
  my $h1 = ($b * $T - 2 * $v * $T**(0.5)) * $X / ($X - $Su);
  my $Si = $Su + ($X - $Su) * exp($h1);

  my $K = 2 * $r / ($v**2 * (1 - exp(-$r * $T)));
  my $d1 = (log($Si / $X) + ($b + $v**2 / 2) * $T) / ($v * $T**(0.5));
  my $Q1 = ( -($n - 1) - (($n - 1)**2 + 4 * $K)**(0.5) ) / 2;
  my $LHS = $X - $Si;
  my $RHS = &GBlackScholes('p', $Si, $X, $T, $r, $b, $v)
           - (1 - exp(($b - $r) * $T) * &CND(-$d1)) * $Si / $Q1;
  my $bi =      - exp(($b - $r) * $T) * &CND(-$d1) * (1 - 1 / $Q1)
           - (1 + exp(($b - $r) * $T) *  &ND(-$d1) / ($v * $T**(0.5))) / $Q1;
  my $epsilon = 0.000001;

  # Newton-Raphson algorithm for finding critical price Si
  while (abs($LHS - $RHS) / $X > $epsilon) {
    $Si = ($X - $RHS + $bi * $Si) / (1 + $bi);
    $d1 = (log($Si / $X) + ($b + $v**2 / 2) * $T) / ($v * $T**(0.5));
    $LHS = $X - $Si;
    $RHS = &GBlackScholes('p', $Si, $X, $T, $r, $b, $v)
          - (1 - exp(($b - $r) * $T) * &CND(-$d1)) * $Si / $Q1;
    $bi =       -exp(($b - $r) * $T) * &CND(-$d1) * (1 - 1 / $Q1)
          - (1 + exp(($b - $r) * $T) *  &ND(-$d1) / ($v * $T**(0.5))) / $Q1;
  }
########## DOUBLE-CHECK AS TO WHETHER OR NOT THE &ND SHOULD BE
########## &CND AS ORIGINALLY ABOVE.

  return $Si;

}


=item BSAmericanApprox

Bjerksund and Stensland (1993) American option approximation.  From Haug (1998).

  $price = BSAmericanApprox($call_put_flag, $S, $X, $T, $r, $b, $v);

C<$call_put_flag> is either 'c' or 'p' for a call or put respectively,
C<$S> is the underlying price, C<$X> is the strike price,
C<$T> is calendar days to expiration/calendar days per year,
C<$r> is the risk-free interest rate, C<$b> is the cost-of-carry term,
and C<$v> is the volatility in percent per annum.

=back

=cut

sub BSAmericanApprox {
  my ($call_put_flag, $S, $X, $T, $r, $b, $v) = @_;

  if ($call_put_flag eq 'c') {
    return BSAmericanCallApprox($S, $X, $T, $r, $b, $v);
  }
  else {
    # use the Bjerksund and Stensland put-call transformation
    return BSAmericanCallApprox($X, $S, $T, $r-$b, -$b, $v);
  }

}


sub BSAmericanCallApprox {
  my ($S, $X, $T, $r, $b, $v) = @_;

  # never optimal to exercise before maturity
  if ($b >= $r) {
    return GBlackScholes('c', $S, $X, $T, $r, $b, $v);
  }
  else {
    my $beta      = ( 0.5 - $b / $v**2 ) + ( ($b/$v**2 - 0.5)**2 + 2*$r / $v**2 )**(0.5);
    my $binfinity = ( $beta/($beta-1) ) * $X;
    my $b0        = &max( $X, ( $r/($r-$b) ) * $X );
    my $ht        = -( $b*$T + 2*$v * $T**(0.5) ) * $b0 / ($binfinity-$b0);
    my $I         = $b0 + ($binfinity-$b0) * (1-exp($ht));
    my $alpha     = ($I-$X) * $I**(-$beta);

    if ($S >= $I) {
      return $S-$X;
    }
    else {
      return $alpha * $S**$beta
             - $alpha * &phi($S, $T, $beta, $I, $I, $r, $b, $v)
             +          &phi($S, $T, 1, $I, $I, $r, $b, $v)
             -          &phi($S, $T, 1, $X, $I, $r, $b, $v)
             -     $X * &phi($S, $T, 0, $I, $I, $r, $b, $v)
             +     $X * &phi($S, $T, 0, $X, $I, $r, $b, $v);
    }
  }

}

sub phi {
  my ($S, $T, $gamma, $H, $I, $r, $b, $v) = @_;

  my $lambda = (-$r + $gamma*$b + 0.5*$gamma*($gamma-1) * $v**2) * $T;
  my $d = -(log($S/$H) + ($b + ($gamma-0.5)*$v**2) *$T) / ($v*$T**(0.5));
  my $kappa = 2 * $b/($v**2) + (2*$gamma -1);

  return exp($lambda) * $S**$gamma * (&CND($d) 
         - ($I/$S)**$kappa * &CND($d - 2*log($I/$S)/($v*$T**(0.5))));

}



=head1 Exotic Options

=over

=item Executive

Jennergren and Naslund (1993) model for executive stock options.

  $price = Executive($call_put_flag, $S, $X, $T, $r, $b, $lambda, $v);

C<$call_put_flag> is either 'c' or 'p' for a call or put respectively,
C<$S> is the underlying price, C<$X> is the strike price,
C<$T> is calendar days to maturity/calendar days per year,
C<$r> is the risk-free interest rate, C<$b> is the cost-of-carry term,
C<$lambda> is the jump rate per year (negative exponential model),
and C<$v> is the volatility in percent per annum.

=cut

sub Executive {
  my ($call_put_flag, $S, $X, $T, $r, $b, $lambda, $v) = @_;
  return exp(-$lambda * $T) * &GBlackScholes($call_put_flag, $S, $X, $T, $r, $b, $v);
}


=item ForwardStartOption

Rubinstein (1990) formula for forward start option.

  $price = ForwardStartOption($call_put_flag, $S, $X, $t1, $T, $r, $b, $v);

C<$call_put_flag> is either 'c' or 'p' for a call or put respectively,
C<$S> is the underlying price, C<$X> is the strike price,
C<$t1> is forward start time,
C<$T> is calendar days to maturity/calendar days per year,
C<$r> is the risk-free interest rate, C<$b> is the cost-of-carry term,
and C<$v> is the volatility in percent per annum.

=cut

sub ForwardStartOption {
  my ($call_put_flag, $S, $alpha, $t1, $T, $r, $b, $v) = @_;
  return $S * exp(($b - $r) * $t1) 
            * &GBlackScholes($call_put_flag, 1, $alpha, $T - $t1, $r, $b, $v);
}


=item TimeSwitchOption

Pechtl (1995) discrete time-switch option.

  $price = TimeSwitchOption($call_put_flag, $S, $X, $a, $T, $m, $dt, $r, $b, $v);

C<$call_put_flag> is either 'c' or 'p' for a call or put respectively,
C<$S> is the underlying price, C<$X> is the strike price,
C<$a> is payoff per C<$dt> time interval,
C<$T> is calendar days to maturity/calendar days per year,
C<$m> is number of time units where option has fulfilled its condition,
C<$r> is the risk-free interest rate, C<$b> is the cost-of-carry term,
and C<$v> is the volatility in percent per annum.


=back

=cut

sub TimeSwitchOption {
  my ($call_put_flag, $S, $X, $a, $T, $m, $dt, $r, $b, $v) = @_;

  my $n = int($T/$dt);

  my $Z = 1;   # for a call
  if ($call_put_flag eq 'p') { $Z = -1; }

  my $sum = 0;
  for (my $i=1; $i<=$n; $i++) {
    my $d = (log($S /$X) + ($b - $v**2 / 2) * $i * $dt) / ($v * ($i * $dt)**0.5);
    $sum += &CND($Z * $d) * $dt;
  }

  return $a * exp(-$r * $T) * $sum + $dt * $a * exp(-$r * $T) * $m;
}

=head2 Chooser Options

=over 4

=item SimpleChooser

Rubinstein (1991) formula for a chooser option.

  $price = SimpleChooser($S, $X, $t1, $T2, $m, $dt, $r, $b, $v);

C<$S> is the underlying price, C<$X> is the strike price,
C<$t1> is the time to choice,
C<$T2> is calendar days to maturity/calendar days per year,
C<$r> is the risk-free interest rate, C<$b> is the cost of carry term,
and C<$v> is the volatility in percent per annum.

=back

=cut

sub SimpleChooser {
  my ($S, $X, $t1, $T2, $m, $dt, $r, $b, $v) = @_;

  my $d = (log($S / $X) + ($b + $v ** 2 / 2) * $T2) / ($v * $T2**0.5);
  my $y = (log($S / $X) + $b * $T2 + $v ** 2 * $t1 / 2) / ($v * $t1**0.5);

  return   $S * exp(($b - $r) * $T2)
              * &CND($d) - $X * exp(-$r * $T2) * &CND($d - $v * $T2**0.5)
         - $S * exp(($b - $r) * $T2)
              * &CND(-$y) + $X * exp(-$r * $T2) * &CND(-$y + $v * $t1**0.5);
}



=head2 Options on Two Different Assets

=over 4

=item TwoAssetCorrelation

Zhang (1995) formula for payoff for a call of max(S2-X2) if S1>X1, 0 otherwise
(put max(X2-S2) if S1<X1, 0 otherwise).

  $price = TwoAssetCorrelation($call_put_flag, $S1, $S2, $X1, $X2, $T,
                                               $b1, $b2, $r, $v1, $v2, $rho);

C<$call_put_flag> is either 'c' or 'p' for a call or put respectively,
C<$S1> and C<$S2> are the underlying prices,
C<$X1> and C<$X2> are the strike prices,
C<$T> is calendar days to maturity/calendar days per year,
C<$b1> and C<$b2> are the cost of carry terms,
C<$r> is the risk-free interest rate,
C<$v1> and C<$v2> are the volatilities in percent per annum,
and C<$rho> is the correlation between the two stocks.

=cut

sub TwoAssetCorrelation {
  my ($call_put_flag, $S1, $S2, $X1, $X2, $T, $b1, $b2, $r, $v1, $v2, $rho) = @_;

  my $y1 = (log($S1 / $X1) + ($b1 - $v1**2 / 2) * $T) / ($v1 * $T**0.5);
  my $y2 = (log($S2 / $X2) + ($b2 - $v2**2 / 2) * $T) / ($v2 * $T**0.5);

  if ($call_put_flag eq 'c') {
    return   $S2 * exp(($b2 - $r) * $T)
                 * &CBND($y2 + $v2 * $T**0.5, $y1 + $rho * $v2 * $T**0.5, $rho)
           - $X2 * exp(-$r * $T) * &CBND($y2, $y1, $rho);
  }
  else {
    return   $X2 * exp(-$r * $T) * &CBND(-$y2, -$y1, $rho)
           - $S2 * exp(($b2 - $r) * $T)
                 * &CBND(-$y2 - $v2 * $T**0.5, -$y1 - $rho * $v2 * $T**0.5, $rho);
  }

}

=item EuropeanExchangeOption

European option to exchange an asset for another asset.  Formula
of Margrabe (1978).

  $price = EuropeanExchangeOption($S1, $S2, $Q1, $Q2, $T,
                                  $r, $b1, $b2, $v1, $v2, $rho);

C<$S1> and C<$S2> are the underlying prices,
C<$Q1> and C<$Q2> are quantities of assets C<$S1> and C<$S2> respectively,
C<$T> is calendar days to maturity/calendar days per year,
C<$r> is the risk-free interest rate,
C<$b1> and C<$b2> are the cost of carry terms,
C<$v1> and C<$v2> are the volatilities in percent per annum,
and C<$rho> is the correlation between the assets.

=cut

sub EuropeanExchangeOption {
  my ($S1, $S2, $Q1, $Q2, $T, $r, $b1, $b2, $v1, $v2, $rho) = @_;

  my $v = ($v1**2 + $v2**2 - 2 * $rho * $v1 * $v2)**0.5;
  my $d1 = (log($Q1 * $S1 / ($Q2 * $S2)) + ($b1 - $b2 + $v**2 / 2) * $T) / ($v * $T**0.5);
  my $d2 = $d1 - $v * $T**0.5;

  return $Q1 * $S1 * exp(($b1 - $r) * $T) * &CND($d1)
       - $Q2 * $S2 * exp(($b2 - $r) * $T) * &CND($d2);
}


=item AmericanExchangeOption

Bjerksund and Stensland (1993) model for American version.  Arguments
as in C<EuropeanExchangeOption>.

=cut

sub AmericanExchangeOption {
  my ($S1, $S2, $Q1, $Q2, $T, $r, $b1, $b2, $v1, $v2, $rho) = @_;

  my $v = ($v1**2 + $v2**2 - 2 * $rho * $v1 * $v2)**0.5;

  return &BSAmericanApprox('c', $Q1 * $S1, $Q2 * $S2, $T, $r - $b2, $b1 - $b2, $v);
}


=item ExchangeExchangeOption

Carr (1988) model for exchange option on exchange option.

  $price = ExchangeExchangeOption($type_flag, $S1, $S2, $q, $t1, $T2,
                                 $r, $b1, $b2, $v1, $v2, $rho) = @_;

See Haug (1998) for a description of the variables.

=cut

sub ExchangeExchangeOption {
  my ($type_flag, $S1, $S2, $q, $t1, $T2, $r, $b1, $b2, $v1, $v2, $rho) = @_;

  my $v = ($v1**2 + $v2**2 - 2 * $rho * $v1 * $v2)**0.5;
  my $I1 = $S1 * exp(($b1 - $r) * ($T2 - $t1)) / ($S2 * exp(($b2 - $r) * ($T2 - $t1)));

  my $id = 2;
  if ($type_flag==1 or $type_flag==2) {
    $id = 1;
  }
  # else $id==2

  my $I = &CriticalPrice($id, $I1, $t1, $T2, $v, $q);
  my $d1 = (log($S1 / ($I * $S2)) + ($b1 - $b2 + $v**2 / 2) * $t1) / ($v * $t1**0.5);
  my $d2 = $d1 - $v * $t1**0.5;
  my $d3 = (log(($I * $S2) / $S1) + ($b2 - $b1 + $v**2 / 2) * $t1) / ($v * $t1**0.5);
  my $d4 = $d3 - $v * $t1**0.5;

  my $y1 = (log($S1 / $S2) + ($b1 - $b2 + $v**2 / 2) * $T2) / ($v * $T2**0.5);
  my $y2 = $y1 - $v * $T2**0.5;
  my $y3 = (log($S2 / $S1) + ($b2 - $b1 + $v**2 / 2) * $T2) / ($v * $T2**0.5);
  my $y4 = $y3 - $v * $T2**0.5;

  if ($type_flag==1) {
    return - $S2 * exp(($b2 - $r) * $T2) * &CBND($d2, $y2, ($t1 / $T2)**0.5 )
           + $S1 * exp(($b1 - $r) * $T2) * &CBND($d1, $y1, ($t1 / $T2)**0.5 )
           - $q * $S2 * exp(($b2 - $r) * $t1) * &CND($d2);
  }
  elsif ($type_flag==2) {
    return   $S2 * exp(($b2 - $r) * $T2) * &CBND($d3, $y2,-($t1 / $T2)**0.5 )
           - $S1 * exp(($b1 - $r) * $T2) * &CBND($d4, $y1,-($t1 / $T2)**0.5 )
           + $q * $S2 * exp(($b2 - $r) * $t1) * &CND($d3);
  }
  elsif ($type_flag==3) {
    return   $S2 * exp(($b2 - $r) * $T2) * &CBND($d3, $y3, ($t1 / $T2)**0.5 )
           - $S1 * exp(($b1 - $r) * $T2) * &CBND($d4, $y4, ($t1 / $T2)**0.5)
           - $q * $S2 * exp(($b2 - $r) * $t1) * &CND($d3);
  }
  elsif ($type_flag==4) {
    return - $S2 * exp(($b2 - $r) * $T2) * &CBND($d2, $y3,-($t1 / $T2)**0.5 )
           + $S1 * exp(($b1 - $r) * $T2) * &CBND($d1, $y4,-($t1 / $T2)**0.5 )
           + $q * $S2 * exp(($b2 - $r) * $t1) * &CND($d2);
  }
}

# Numerical search algorithm to find critical price I
sub CriticalPrice {
  my ($id, $I1, $t1, $T2, $v, $q) = @_;

  my $Ii = $I1;
  my $yi = &CriticalPart3($id, $Ii, $t1, $T2, $v);
  my $di = &CriticalPart2($id, $Ii, $t1, $T2, $v);

  my $epsilon = 0.00001;
  while (abs($yi-$q) > $epsilon) {
    $Ii = $Ii - ($yi - $q) / $di;
    $yi = &CriticalPart3($id, $Ii, $t1, $T2, $v);
    $di = &CriticalPart2($id, $Ii, $t1, $T2, $v);
  }

  return $Ii;

}

sub CriticalPart2 {
  my ($id, $I, $t1, $T2, $v) = @_;

  if ($id==1) {
    my $z1 = (log($I) + $v**2 / 2 * ($T2 - $t1)) / ($v * ($T2 - $t1)**0.5);
    return &CND($z1);

  }
  elsif ($id==2) {
    my $z2 = (-log($I) - $v**2 / 2 * ($T2 - $t1)) / ($v * ($T2 - $t1)**0.5);
    return -&CND($z2);
  }
}

#  To do:  1.  Check to see whether the duplication calculation is necessary.
sub CriticalPart3 {
  my ($id, $I, $t1, $T2, $v) = @_;

  if ($id==1) {
    my $z1 = (log($I) + $v**2 / 2 * ($T2 - $t1)) / ($v * ($T2 - $t1)**0.5);
    my $z2 = (log($I) - $v**2 / 2 * ($T2 - $t1)) / ($v * ($T2 - $t1)**0.5);
    return $I * &CND($z1) - &CND($z2);
  }
  elsif ($id==2) {
    my $z1 = (-log($I) + $v**2 / 2 * ($T2 - $t1)) / ($v * ($T2 - $t1)**0.5);
    my $z2 = (-log($I) - $v**2 / 2 * ($T2 - $t1)) / ($v * ($T2 - $t1)**0.5);
    return &CND($z1) - $I * &CND($z2);
  }
}


=item OptionsOnTheMaxMin

Options on the maximum or minimum of two risky assets

  $price = OptionsOnTheMaxMin($type_flag, $S1, $S2, $X, $T,
                              $r, $b1, $b2, $v1, $v2, $rho);

See Haug (1998) for a description of the variables.

=cut

sub OptionsOnTheMaxMin {
  my ($type_flag, $S1, $S2, $X, $T, $r, $b1, $b2, $v1, $v2, $rho) = @_;

  my $v    = ($v1**2 + $v2**2 - 2 * $rho * $v1 * $v2)**0.5;
  my $rho1 = ($v1 - $rho * $v2) / $v;
  my $rho2 = ($v2 - $rho * $v1) / $v;
  my $d    = ( log($S1/$S2) + ($b1 - $b2 + $v**2 / 2) * $T ) / ($v * $T**0.5);
  my $y1   = ( log($S1/$X)  + ($b1 + $v1**2 / 2     ) * $T ) / ($v1 * $T**0.5);
  my $y2   = ( log($S2/$X)  + ($b2 + $v2**2 / 2     ) * $T ) / ($v2 * $T**0.5);

  if ($type_flag eq 'cmin') {
    return   $S1 * exp(($b1 - $r) * $T) * &CBND($y1, -$d, -$rho1)
           + $S2 * exp(($b2 - $r) * $T) * &CBND($y2, $d - $v * $T**0.5, -$rho2)
           - $X * exp(-$r * $T) * &CBND($y1 - $v1 * $T**0.5, $y2 - $v2 * $T**0.5, $rho);
  }
  elsif ($type_flag eq 'cmax') {
    return   $S1 * exp(($b1 - $r) * $T) * &CBND($y1, $d, $rho1)
           + $S2 * exp(($b2 - $r) * $T) * &CBND($y2, -$d + $v * $T**0.5, $rho2)
           - $X * exp(-$r * $T) * (1 - &CBND(-$y1 + $v1 * $T**0.5, -$y2 + $v2 * $T**0.5, $rho));
  }
  elsif ($type_flag eq 'pmin') {
    return $X * exp(-$r * $T) - $S1 * exp(($b1 - $r) * $T)
           + &EuropeanExchangeOption($S1, $S2, 1, 1, $T, $r, $b1, $b2, $v1, $v2, $rho)
           + &OptionsOnTheMaxMin('cmin', $S1, $S2, $X, $T, $r, $b1, $b2, $v1, $v2, $rho);
  }
  elsif ($type_flag eq 'pmax') {
    return $X * exp(-$r * $T) - $S2 * exp(($b2 - $r) * $T)
           - &EuropeanExchangeOption($S1, $S2, 1, 1, $T, $r, $b1, $b2, $v1, $v2, $rho)
           + &OptionsOnTheMaxMin('cmax', $S1, $S2, $X, $T, $r, $b1, $b2, $v1, $v2, $rho);
  }

}

=item SpreadApproximation

  $price = SpreadApproximation($call_put_flag, $f1, $f2, $X, $T, $r, $v1, $v2, $rho);

C<$call_put_flag> is either 'c' or 'p' for a call or put respectively,
C<$f1> and C<$f2> are futures contract prices, C<$X> is the strike price,
C<$T> is calendar days to maturity/calendar days per year,
C<$r> is the risk-free interest rate,
C<$v1> and C<$v2> are the volatilities, and
C<$rho> is the correlation between the two futures contracts.

=back

=cut

sub SpreadApproximation {
  my ($call_put_flag, $f1, $f2, $X, $T, $r, $v1, $v2, $rho) = @_;

  my $v = ($v1**2 + ($v2 * $f2 / ($f2 + $X))**2 - 2 * $rho * $v1 * $v2 * $f2 / ($f2 + $X))**0.5;
  my $F = $f1 / ($f2 + $X);

  return GBlackScholes($call_put_flag, $F, 1, $T, $r, 0, $v) * ($f2 + $X);
}


=head1 Numerical Methods

=head2 Binomial Options Pricing

=over 4

=item CRRBinomial

Cox-Ross-Rubinstein (1979) binomial tree model.

  $price = CRRBinomial($ame_eur_flag, $call_put_flag, $S, $X, $T, $r, $b, $n, $v);

C<$ame_eur_flag> is either 'a' or 'e' for American or European options
respectively,
C<$call_put_flag> is either 'c' or 'p' for a call or put respectively,
C<$S> is the underlying price, C<$X> is the strike price,
C<$T> is calendar days to maturity/calendar days per year,
C<$r> is the risk-free interest rate, C<$b> is the cost of carry term,
C<$n> is the number of time steps,
and C<$v> is the volatility in percent per annum.

=back

=cut

sub CRRBinomial {
  my ($ame_eur_flag, $call_put_flag, $S, $X, $T, $r, $b, $n, $v) = @_;

  my $z = -1;      # put
  if ($call_put_flag eq 'c') { $z = 1; }


  my $dt = $T / $n;
  my $u  = exp($v * $dt**(0.5));
  my $d  = 1 / $u;
  my $p  = (exp($b * $dt) - $d) / ($u - $d);
  my $Df = exp(-$r * $dt);

  # starting conditions
  my @option_value = ();
  for (my $i=0; $i<=$n; $i++) {
    $option_value[$i] = &max( 0, $z * ($S * $u**$i * $d**($n-$i) - $X) );
  }

  # dynamic programming
  for (my $j=$n-1; $j>=0; $j--) {
    for (my $i=0; $i<=$j; $i++) {

      # European option
      if ($ame_eur_flag eq 'e') {
        $option_value[$i] =  ($p * $option_value[$i+1] + (1 - $p) * $option_value[$i]) * $Df;
      }

      # American option
      else {
        $option_value[$i] =
          &max( ($z * ($S * $u**$i * $d**(abs($i-$j)) - $X)),
                ($p * $option_value[$i+1] + (1 - $p) * $option_value[$i]) * $Df
              );
      }

    }
  }

  return $option_value[0];
}


=head2 Trinomial Trees

=over 4

=item TrinomialTree

The trinomial tree (Boyle, 1986).

  $price = TrinomialTree($ame_eur_flag, $call_put_flag, $S, $X, $T, $r, $b, $n, $v);

Inputs as in C<CRRBinomial>.

=back

=cut

sub TrinomialTree {
  my ($ame_eur_flag, $call_put_flag, $S, $X, $T, $r, $b, $n, $v) = @_;

  my @option_value = ();

  my $z = -1;        # for put

  if ($call_put_flag eq 'c') { $z = 1; }

  # initialize
  my $dt = $T / $n;
  my $u = exp($v * (2 * $dt)**0.5);
  my $d = exp(-$v * (2 * $dt)**0.5);
  my $pu = (  (exp($b * $dt / 2) - exp(-$v * ($dt / 2)**0.5)) 
            / (exp($v * ($dt / 2)**0.5) - exp(-$v * ($dt / 2)**0.5))
            )**2;
  my $pd = (  (exp($v * ($dt / 2)**0.5) - exp($b * $dt / 2))
            / (exp($v * ($dt / 2)**0.5) - exp(-$v * ($dt / 2)**0.5))
            )**2;
  my $pm = 1 - $pu - $pd;
  my $Df = exp(-$r * $dt);

  # initialize leaf nodes
  for (my $i=0; $i<=(2*$n); $i++) {
    $option_value[$i] = &max(0, $z * ($S * $u**max($i-$n, 0) * $d**max($n*2-$n-$i, 0) - $X));
  }

  for (my $j=$n-1; $j>=0; $j--) {
    for (my $i=0; $i<=($j*2); $i++) {

      if ($ame_eur_flag eq 'e') {
        $option_value[$i] = ($pu * $option_value[$i+2] + $pm * $option_value[$i+1] 
                             + $pd * $option_value[$i]) * $Df;
      }
      else {
        $option_value[$i] =
          &max( ($z * ($S * $u**&max($i-$j, 0) * $d**&max($j*2-$j-$i, 0) - $X)),
                ($pu * $option_value[$i+2] + $pm * $option_value[$i+1]
                                           + $pd * $option_value[$i]          ) * $Df
               );
      }
    }
  }

  return $option_value[0];

}


=head1 Volatility and Correlation

=over 4

=item ImpliedVolatilityB

Uses the bisection algorithm to numerically determing the volatility of an option
given its fair market value.

  $volatility = ImpliedVolatilityB($algorithm, $cm, @args);
  $volatility = ImpliedVolatilityB('BlackScholes', $cm, $call_put_flag, $S, $X, $T, $r);

C<$algorithm> is one of (BlackScholes, GBlackScholes),
C<$cm> is the fair market value of the option,
and C<@args> is the list of arguments to the option pricing algorithm
excepting the final volatility value C<$v> -- order is important!

=cut

# This is one place where it would be very helpful to have named
# arguments for subroutines.  With this, one could write something
# like:
#   ImpliedVolatilityB(ALGORITHM => 'GBlackScholes',
# followed by the rest of the named arguments.  The code would look
# something like:
# sub ImpliedVolatility {
#   my %args = (@_);
#
#   # set up dispatch table
#   my %dispatch = ('GBlackScholes' => \&GBlackScholes,
#                  );
#   my $algorithm = $dispatch{$args{ALGORITHM}};
#
#   # e.g.
#   my $clow  = $algorithm->(%args, VOLATILITY => $vlow);
#   my $chigh = $algorithm->(%args, VOLATILITY => $vhigh);
#   etc. etc. etc.
# }
#
# This would prevent worries about where volatility should
# come in, might be useful elsewhere, and is more self-documenting.
# On the negative side, it's a lot less compact.
#
# For now, though, we will implement something similar but using
# a convention that the volatility is the last parameter in the
# algorithm to be called.  Another fix would be to tabulate the
# position of the volatility parameter in each of the supported
# methods, then splice it in as necessary given the other parameters.


sub ImpliedVolatilityB {
  my $algorithm = shift;   # the option pricing algorithm
  my $cm = shift;          # the fair market value

  # the rest of the arguments will be passed in as is
  my @routine_args = @_;

  # determine the correct algorithm to use
  my %dispatch = ('GBlackScholes' => \&GBlackScholes,
                  'BlackScholes'  => \&BlackScholes,
                  );
  my $routine = $dispatch{$algorithm};

  # note that Gaarder's code uses .01 for vlow, which did choke
  # on a real-life OEX option calculation -- this should
  # probably be made an input for the user or made a default
  # (which again argues for using a hash-style function input
  # all around)
  my ($vlow, $vhigh, $epsilon) = (0.001, 1, 0.000001);

  # get initial estimates which bracket the desired value
  my $clow  = $routine->(@routine_args, $vlow);
  my $chigh = $routine->(@routine_args, $vhigh);
  my $vi = $vlow + ($cm-$clow)*($vhigh-$vlow)/($chigh-$clow);

  # Debug output
  #print "clow $clow chigh $chigh vlow $vlow vhigh $vhigh vi $vi\n";

  # perform the bisection algorithm
  while (abs ($cm - $routine->(@routine_args, $vi)) > $epsilon) {
    # set the end-point
    if ($routine->(@routine_args, $vi) < $cm) {
      $vlow = $vi;
    }
    else {
      $vhigh = $vi;
    }

    # get the new mid-point
    $clow  = $routine->(@routine_args, $vlow);
    $chigh = $routine->(@routine_args, $vhigh);
    $vi = $vlow + ($cm-$clow)*($vhigh-$vlow)/($chigh-$clow);

    # Debug output
    #print "clow $clow chigh $chigh vlow $vlow vhigh $vhigh vi $vi\n";

  }

  return $vi;
}


=item ImpliedVolatilityNR

Uses the Newton-Raphson algorithm to numerically determing the volatility of an
option priced using the generalized Black-Scholes model, given its fair market
value.

  $volatility = ImpliedVolatilityNR($call_put_flag, $S, $X, $T, $r, $cm);

C<$call_put_flag> is either 'c' or 'p' for a call or put respectively,
C<$S> is the underlying price, C<$X> is the strike price,
C<$T> is calendar days to maturity/calendar days per year,
C<$r> is the risk-free interest rate,
and C<$cm> is the fair market value of the option,

=back

=cut

sub ImpliedVolatilityNR {
  my ($call_put_flag, $S, $X, $T, $r, $cm) = @_;

  my $epsilon = 0.00000000001;

  # Manaster and Koehler (1982) seed value
  my $vi = ( abs(log($S/$X) + $r*$T) * 2/$T )**0.5;

  my $ci = GBlackScholes($call_put_flag, $S, $X, $T, $r, $r, $vi);
  my $vegai = GVega($S, $X, $T, $r, $r, $vi);

  while (abs($cm-$ci)>$epsilon) {
    $vi    = $vi - ($ci-$cm)/$vegai;
    $ci    = GBlackScholes($call_put_flag, $S, $X, $T, $r, $r, $vi);
    $vegai = GVega($S, $X, $T, $r, $r, $vi);
  }

  return $vi;

}


=head1 Utility Routines

From E. G. Haug (1998) The Complete Guide to Option Pricing
Formulas.  McGraw-Hill.

=over

=item CND

Approximate the cumulative normal distribution.  That is, the value
of the integral of the standard normal density from minus infinity
to C<$x>.

  $p = &CND($x);

=cut

sub CND {
  my $x = shift;     # the percentile under consideration

  my $Pi = 3.141592653589793238;

  # Taylor series coefficients
  my ($a1, $a2, $a3, $a4, $a5) =
    (0.319381530, -0.356563782, 1.781477937, -1.821255978,  1.330274429);

  # use symmetry to perform the calculation to the right of 0
  my $L = abs($x);

  my $k = 1/( 1 + 0.2316419*$L);

  my $CND = 1 - 1/(2*$Pi)**0.5 * exp(-$L**2/2)
                * ($a1*$k + $a2*$k**2 + $a3*$k**3 + $a4*$k**4 + $a5*$k**5);

  # then return the appropriate value
  return ($x >= 0) ? $CND : 1-$CND;

}


=item ND

The standard normal density evaluated at C<$x>.

  $p = &ND($x);

=cut

sub ND {
  my $x = shift;
  my $Pi = 3.141592653589793238;
  return ( 1 / (2 * $Pi)**(0.5) ) * exp(-$x**2 / 2);
}


=item CBND

Standardized cumulative bivariate normal distribution with covariance
C<$rho> evaluated at C<$a, $b>.  That is, the double integral evaluated
from minus infinity to each argument respectively.

  $p = CBND($a, $b, $rho);

=cut

sub CBND {
  my ($a, $b, $rho) = @_;

  my $Pi = 3.141592653589793238;

  # Taylor series coefficients
  my @X = (0.24840615, 0.39233107, 0.21141819, 0.03324666, 0.00082485334);
  my @y = (0.10024215, 0.48281397, 1.0609498, 1.7797294, 2.6697604);

  my $a1 = $a / (2 * (1-$rho**2))**(0.5);
  my $b1 = $b / (2 * (1-$rho**2))**(0.5);

  if ($a<=0 and $b<=0 and $rho<=0) {

    # do a Taylor series approximation
    my $sum = 0;

    for (my $i=0; $i<=4; $i++) {
      for (my $j=0; $j<=4; $j++) {
        $sum += $X[$i] * $X[$j] * exp( $a1 * (2 * $y[$i] - $a1)
                + $b1 * (2 * $y[$j] - $b1) + 2 * $rho * ($y[$i] - $a1) * ($y[$j] - $b1) );
      }
    }

    return (1 - $rho**2)**(0.5) / $Pi * $sum;

  }

  elsif ($a<=0 and $b>=0 and $rho>=0) {
    return &CND($a) - &CBND($a, -$b, -$rho);
  }

  elsif ($a>=0 and $b<=0 and $rho>=0) {
    return &CND($b) - &CBND(-$a, $b, -$rho);
  }

  elsif ($a>=0 and $b>=0 and $rho<=0) {
    return &CND($a) + &CND($b) - 1 + &CBND(-$a, -$b, $rho);
  }

  else {
    my $rho1 = ($rho * $a - $b) * &sgn($a) / ($a**2 - 2 * $rho * $a * $b + $b**2)**(0.5);
    my $rho2 = ($rho * $b - $a) * &sgn($b) / ($a**2 - 2 * $rho * $a * $b + $b**2)**(0.5);
    my $delta = (1 - &sgn($a) * &sgn($b)) / 4;
    return &CBND($a, 0, $rho1) + &CBND($b, 0, $rho2) - $delta;
  }

}

=item max

Maximum of an array.

  $max = max(@array);

=cut

sub max {
   my $max = shift;     # set $max to something...
   foreach (@_) { $max = $_ if $_ > $max };
   return $max;
}


=item sgn

Sign of a number.  Zero returns zero.

  $sign = sgn($num);

=cut

sub sgn {
  my $x = shift;

  if    ($x>0) { return  1; }
  elsif ($x<0) { return -1; }
  else         { return  0; }
}

=item factorial

Factorial function -- naive implementation.

  print factorial(4);    # prints 24

=back

=cut

sub factorial {
  my $n = shift;
  return $n == 0 ? 1 : $n * factorial($n-1);
}

1;

=head1 AUTHOR

Jerome V. Braun <jerome.braun@kmri.com>

=head1 VERSION

Version 0.03

Document generated from
  $Id: Options.pm,v 1.16 2000/11/25 12:47:07 jvbraun Exp $

=head1 COPYRIGHT

Copyright (c) 1999 Jerome V. Braun. All rights reserved.
This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself (see http://www.perl.com/perl/misc/Artistic.html).

=cut

__END__


