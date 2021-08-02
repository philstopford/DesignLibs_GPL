using System;

namespace Burkardt
{
    public static class EIQF
    {
        public static double eiqf(int nt, double[] t, int[] mlt, double[] wts, int nwts, int[] ndx,
        int key, Func <double, int, double > f )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EIQF evaluates an interpolatory quadrature formula.
        //
        //  Discussion:
        //
        //   The knots, weights and integrand are supplied.
        //
        //   All knots with nonzero NDX are used.
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
        //    Input, int MLT[NT], the multiplicity of the knots.
        //
        //    Input, double WTS[NWTS], the weights.
        //
        //    Input, int NWTS, the number of weights.
        //
        //    Input, int NDX[NT], used to index the array WTS.  
        //    If KEY = 1, then NDX need not be preset.  For more details see the 
        //    comments in CAWIQ.
        //
        //    Input, int KEY, indicates the structure of the WTS
        //    array.  It will normally be set to 1.  This will cause the weights to be 
        //    packed sequentially in array WTS.  For more details see the comments 
        //    in CAWIQ.
        //
        //    Input, double F ( double X, int I ), the name of a routine which
        //    evaluates the function and some of its derivatives.  The routine
        //    must return in F the value of the I-th derivative of the function
        //    at X.  The highest value of I will be the maximum value in MLT minus
        //    one.  The value X will always be a knot.
        //
        //    Output, double EIQF, the value of the quadrature formula 
        //    applied to F.
        //
        {
            int i;
            int j;
            int l;
            double p;
            double qfsum;

            l = Math.Abs(key);

            if (l < 1 || 4 < l)
            {
                Console.WriteLine("");
                Console.WriteLine("EIQF - Fatal error!");
                Console.WriteLine("  Magnitude of KEY must be between 1 and 4.");
                return (1);
            }

            qfsum = 0.0;
            for (j = 0; j < nt; j++)
            {
                l = Math.Abs(ndx[j]);
                if (l != 0)
                {
                    p = 1.0;
                    for (i = 0; i < mlt[j]; i++)
                    {
                        qfsum = qfsum + wts[l + i - 1] * f(t[j], i) / p;
                        if (key <= 0)
                        {
                            p = p * (i + 1);
                        }
                    }
                }
            }

            return qfsum;
        }
    }
}