namespace Burkardt.PolynomialNS
{
    public static class Chebyshev
    {
        public static void cheb ( int deg, double pt, ref double[] tcheb, int index = 0 )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEB computes normalized Chebyshev polynomials.
            //
            //  Discussion:
            //
            //    This subroutine computes the array TCHEB of normalized Chebyshev 
            //    polynomials from degree 0 to DEG:
            //      T_0(x)=1, 
            //      T_j(x) = sqrt(2) * cos ( j * acos(x) ) 
            //    at the point x = PT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //  
            //  Modified:
            //
            //    14 February 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
            //    Marco Vianello.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Marco Caliari, Stefano de Marchi, Marco Vianello,
            //    Algorithm 886:
            //    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
            //    ACM Transactions on Mathematical Software,
            //    Volume 35, Number 3, October 2008, Article 21, 11 pages.
            //
            //  Parameters:
            //
            //    Input, int DEG, the degree.
            //    0 <= DEG.
            //
            //    Input, double PT, the evaluation point.
            //
            //    Output, double TCHEB[DEG+1], the value of the normalized
            //    Chebyshev polynomials of degrees 0 through DEG at the point PT.
            //
        {
            int j;
            const double sqrt2 = 1.4142135623730951;

            if ( deg < 0 )
            {
                return;
            }

            tcheb[index + 0] = 1.0;

            if ( deg < 1 )
            {
                return;
            }

            tcheb[index + 1] = sqrt2 * pt;
 
            if ( deg < 2 )
            {
                return;
            }

            tcheb[index + 2] = 2.0 * pt * tcheb[index + 1] - sqrt2 * tcheb[index + 0];
            //
            //  Chebyshev recurrence.
            //
            for ( j = 3; j <= deg; j++ )
            {
                tcheb[index + j] = 2.0 * pt * tcheb[index + j-1] - tcheb[index + j-2];
            }
        }
    }
}