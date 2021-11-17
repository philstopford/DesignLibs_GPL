using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8vec_mm_to_01(int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_MM_TO_01 shifts and rescales data to lie within [0,1].
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    On input, A contains the original data.  On output, A has been shifted
        //    and scaled so that all entries lie between 0 and 1.
        //
        //    The formula is:
        //
        //      A(I) := ( A(I) - AMIN ) / ( AMAX - AMIN )
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
        //    Input, int N, the number of data values.
        //
        //    Input/output, double A[N], the data to be rescaled.
        //
    {
        int i;

        double amax = r8vec_max(n, a);
        double amin = r8vec_min(n, a);

        if (Math.Abs(amin - amax) <= double.Epsilon)
        {
            for (i = 0; i < n; i++)
            {
                a[i] = 0.5;
            }
        }
        else
        {
            for (i = 0; i < n; i++)
            {
                a[i] = (a[i] - amin) / (amax - amin);
            }
        }
    }

    public static double[] r8vec_mm_to_cd(int n, double[] a, double bmin, double bmax)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_MM_TO_CD shifts and rescales data to lie within a given pair of bounds.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The mininum entry of A is mapped to BMIN, the maximum entry
        //    to BMAX, and values in between are mapped linearly.
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
        //    Input, int N, the number of data values.
        //
        //    Input, double A[N], the data to be remapped.
        //
        //    Input, double BMIN, BMAX, the values to which min(A) and max(A)
        //    are to be assigned.
        //
        //    Output, double R8VEC_MM_TO_CD[N], the remapped data.
        //
    {
        int i;

        double[] b = new double[n];

        if (Math.Abs(bmax - bmin) <= double.Epsilon)
        {
            for (i = 0; i < n; i++)
            {
                b[i] = bmin;
            }

            return b;
        }

        double amax = r8vec_max(n, a);
        double amin = r8vec_min(n, a);

        if (Math.Abs(amin - amax) <= double.Epsilon)
        {
            for (i = 0; i < n; i++)
            {
                b[i] = 0.5 * (bmax + bmin);
            }
        }
        else
        {
            for (i = 0; i < n; i++)
            {
                b[i] = ((amax - a[i]) * bmin
                        + (a[i] - amin) * bmax)
                       / (amax - amin);
            }
        }

        return b;
    }


}