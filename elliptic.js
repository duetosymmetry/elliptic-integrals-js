/* The arithmetic-geometric mean of two non-negative numbers */
Math.agm = function(a0,g0)
{
  var an = (a0+g0)/2;
  var gn = Math.sqrt(a0*g0);

  while (Math.abs(an-gn) > 1e-15)
  {
    a0 = 0.5 * (an + gn);
    g0 = Math.sqrt(an*gn);

    an = a0;
    gn = g0;
  };

  return an;
};

/* EllipticK(m) - The complete elliptic integral of the first type.
 * The argument is the *parameter* m = k^2, where k is the *modulus*.
 * The parameter must satisfy m < 1.
 * In terms of the integral definition, we have
 * K(m) = \int_0^{\pi/2} \frac{d\theta}{\sqrt{1 - m (\sin\theta)^2}}
 * See http://dlmf.nist.gov/19.8.E5 for this method.
 */
Math.EllipticK = function( m )
{
  var kprime = Math.sqrt(1 - m);

  return 0.5 * Math.PI / Math.agm(1, kprime);
};

/* EllipticE(m) - The complete elliptic integral of the second type.
 * The argument is the *parameter* m = k^2, where k is the *modulus*.
 * The parameter must satisfy m < 1.
 * In terms of the integral definition, we have
 * E(m) = \int_0^{\pi/2} \sqrt{1 - m (\sin\theta)^2} d\theta
 * This algorithm comes from:
 * http://scitation.aip.org/content/aip/journal/jap/34/9/10.1063/1.1729771
 * Garrett, Milan Wayne. "Calculation of fields, forces, and mutual
 * inductances of current systems by elliptic integrals." Journal of
 * Applied Physics 34.9 (1963): 2567-2573.
 * See Eqs. (18)-(21)
 */
Math.EllipticE = function( m )
{
  var kprime = Math.sqrt(1. - m);

  var a0 = 1., g0 = kprime;

  var an = a0, gn = g0;

  var twoPow = 0.25;

  var partialSum = 1. - 0.5 * m;

  do {
    partialSum -= twoPow * (an - gn)*(an - gn);
    twoPow *= 2.; 

    a0 = 0.5 * (an + gn);
    g0 = Math.sqrt(an*gn);

    an = a0;
    gn = g0;

  } while (Math.abs(an-gn) > 1e-15);

  return 0.5 * Math.PI * partialSum / an; 
};

/* EllipticPi(n,m) - The complete elliptic integral of the third type.
 * The arguments are the characteristic n,
 * and the *parameter* m = k^2, where k is the *modulus*.
 * The arguments must satisfy n < 1, m < 1.
 * In terms of the integral definition, we have
 * Pi(n,m) = \int_0^{\pi/2} \frac{1}{(1-n(\sin\theta)^2)\sqrt{1 - m (\sin\theta)^2}} d\theta
 * This algorithm comes from:
 * http://scitation.aip.org/content/aip/journal/jap/34/9/10.1063/1.1729771
 * Garrett, Milan Wayne. "Calculation of fields, forces, and mutual
 * inductances of current systems by elliptic integrals." Journal of
 * Applied Physics 34.9 (1963): 2567-2573.
 * See Eqs. (18)-(21)
 */
Math.EllipticPi = function( n, m )
{
  var kprime = Math.sqrt(1. - m);

  var a0 = 1., g0 = kprime, zeta0 = 0.;

  var an = a0, gn = g0,
      deltan = (1. - n)/kprime,
      epsn =  n/(1. - n), zetan = zeta0;

  do {

    zeta0 = 0.5 * (epsn + zetan);
    epsn = (deltan*epsn + zetan)/(1. + deltan);
    zetan = zeta0;

    a0 = 0.5 * (an + gn);
    g0 = Math.sqrt(an*gn);

    an = a0;
    gn = g0;

    deltan = 0.25 * gn / an * (2. + deltan + 1./deltan);

  } while ((Math.abs(an-gn) > 1e-15) || (Math.abs(deltan - 1.) > 1e-15));

  return 0.5 * Math.PI / an * (1. + zetan);
};
