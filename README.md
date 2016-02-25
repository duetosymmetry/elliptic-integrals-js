# elliptic-integrals-js
Complete elliptic integrals in javascript

Here I've implemented complete elliptic integrals of the first, second, and third kind in javascript.
The implementation follows an iteration scheme based on the convergence of the arithmetic-geometric mean, which converges at least quadratically with number of iterations---therefore effectively doubling the number of digits each iteration.
These iteration schemes come from [Garrett, Milan Wayne, Journal of Applied Physics 34.9 (1963): 2567-2573](http://dx.doi.org/10.1063/1.1729771), Eqs. (18)-(21).

The functions extend the `Math` object with:
* `Math.agm(a,g)`: arithmetic-geometric mean of two non-negative numbers
* `Math.EllipticK(m)`: Complete elliptic integral of the first type
* `Math.EllipticE(m)`: Complete elliptic integral of the second type
* `Math.EllipticPi(n,m)`: Complete elliptic integral of the third type

The arguments are the *parameter* `m`, which is related to the *modulus* `k` via `m=k^2`; and the *characteristic* `n`.
The algorithms are valid for `m<1`, `n<1`.

To be completely clear, the functions are computing the following integrals:
* $ K(m) = \int_0^{\pi/2} \frac{d\theta}{\sqrt{1 - m (\sin\theta)^2}} $
* $ E(m) = \int_0^{\pi/2} \sqrt{1 - m (\sin\theta)^2} d\theta $
* $ Pi(n,m) = \int_0^{\pi/2} \frac{1}{(1-n(\sin\theta)^2)\sqrt{1 - m (\sin\theta)^2}} d\theta $

This agrees with the conventions of `Mathematica`.
