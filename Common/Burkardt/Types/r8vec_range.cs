using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8vec_range(int n, double[] x, double xmin, double xmax, double[] y,
                ref double ymin, ref double ymax)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_RANGE finds the range of Y's within a restricted X range.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    The routine is given a set of pairs of points (X,Y), and a range
            //    XMIN to XMAX of valid X values.  Over this range, it seeks
            //    YMIN and YMAX, the minimum and maximum values of Y for
            //    valid X's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the array.
            //
            //    Input, double X[N], the X array.
            //
            //    Input, double XMIN, XMAX, the range of X values to check.
            //
            //    Input, double Y[N], the Y array.
            //
            //    Output, double *YMIN, *YMAX, the range of Y values whose
            //    X value is within the X range.
            //
        {
            int i;
            const double r8_huge = 1.79769313486231571E+308;

            ymin = r8_huge;
            ymax = -r8_huge;

            for (i = 0; i < n; i++)
            {
                if (xmin <= x[i] && x[i] <= xmax)
                {
                    ymin = Math.Min(ymin, y[i]);
                    ymax = Math.Max(ymax, y[i]);
                }
            }
        }

        public static void r8vec_range_2(int n, double[] a, ref double amin, ref double amax)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_RANGE_2 updates a range to include a new R8VEC
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    Given a range AMIN to AMAX, and an array A, the routine will
            //    decrease AMIN if necessary, or increase AMAX if necessary, so that
            //    every entry of A is between AMIN and AMAX.
            //
            //    However, AMIN will not be increased, nor AMAX decreased.
            //
            //    This routine may be used to compute the maximum and minimum of a
            //    collection of arrays one at a time.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the array.
            //
            //    Input, double A[N], the array.
            //
            //    Input/output, double *AMIN, *AMAX.  On input, the
            //    current legal range of values for A.  On output, AMIN and AMAX
            //    are either unchanged, or else "widened" so that all entries
            //    of A are within the range.
            //
        {
            amax = Math.Max(amax, r8vec_max(n, a));
            amin = Math.Min(amin, r8vec_min(n, a));
        }

    }
}