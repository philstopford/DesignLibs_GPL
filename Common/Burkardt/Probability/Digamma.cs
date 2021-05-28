using System;

namespace Burkardt.Probability
{
    public static class Digamma
    {
        public static double digamma ( double x )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIGAMMA calculates the digamma or Psi function.
        //
        //  Discussion:
        //
        //    DiGamma ( X ) = d ( log ( Gamma ( X ) ) ) / dX
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 March 2016
        //
        //  Author:
        //
        //    Original FORTRAN77 by J Bernardo.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    J Bernardo,
        //    Algorithm AS 103:
        //    Psi ( Digamma ) Function,
        //    Applied Statistics,
        //    Volume 25, Number 3, pages 315-317, 1976.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the digamma function.
        //    0 < X.
        //
        //    Output, double DIGAMMA, the value of the digamma function at X.
        //
        {
            double c = 8.5;
            double euler_mascheroni = - 0.57721566490153286060;
            double value;
            //
            //  Check the input.
            //
            if ( x <= 0.0 )
            {
                value = 0.0;
                return value;
            }
            //
            //  Use approximation for small argument.
            //
            if ( x <= 0.000001 )
            {
                value = - euler_mascheroni - 1.0 / x + 1.6449340668482264365 * x;
                return value;
            }
            //
            //  Reduce to DIGAMA(X + N).
            //
            value = 0.0;
            double x2 = x;
            while ( x2 < c )
            {
                value = value - 1.0 / x2;
                x2 = x2 + 1.0;
            }
            //
            //  Use Stirling's (actually de Moivre's) expansion.
            //
            double r = 1.0 / x2;
            value = value + Math.Log ( x2 ) - 0.5 * r;

            r = r * r;

            value = value 
                    - r * ( 1.0 / 12.0 
                            - r * ( 1.0 / 120.0 
                                    - r * ( 1.0 / 252.0 
                                            - r * ( 1.0 / 240.0
                                                    - r * ( 1.0 / 132.0 ) ) ) ) );

            return value;
        }
    }
}