using System;

namespace Burkardt.Quadrature
{
    public static class EIQFS
    {
        public static double eiqfs ( int nt, double[] t, double[] wts, Func < double, int, double > f )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EIQFS evaluates a quadrature formula defined by CLIQF or CLIQFS.
        //
        //  Discussion:
        //
        //    This routine evaluates an interpolatory quadrature formula with all knots 
        //    simple and all knots included in the quadrature.  This routine will be used
        //    typically after CLIQF or CLIQFS has been called.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 January 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Sylvan Elhay, Jaroslav Kautsky,
        //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
        //    Interpolatory Quadrature,
        //    ACM Transactions on Mathematical Software,
        //    Volume 13, Number 4, December 1987, pages 399-415.
        //
        //  Parameters:
        //
        //    Input, int NT, the number of knots.
        //
        //    Input, double T[NT], the knots.
        //
        //    Input, double WTS[NT], the weights.
        //
        //    Input, double F ( double X, int I ), the name of a routine which
        //    evaluates the function and some of its derivatives.  The routine
        //    must return in F the value of the I-th derivative of the function
        //    at X.  The value of I will always be 0.  The value X will always be a knot.
        //
        //    Output, double EIQFS, the value of the quadrature formula 
        //    applied to F.
        //
        {
            int j;
            double qfsum;

            qfsum = 0.0;
            for ( j = 0; j < nt; j++ )
            {
                qfsum = qfsum + wts[j] * f ( t[j], 0 );
            }
            return qfsum;
        }
    }
}