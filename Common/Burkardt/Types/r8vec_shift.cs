using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8vec_shift ( int shift, int n, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SHIFT performs a shift on an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8 values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int SHIFT, the amount by which each entry is to
        //    be shifted.
        //
        //    Input, int N, the length of the vector.
        //
        //    Input/output, double X[N], the vector to be shifted.
        //
    {
        int i;

        double[] y = new double[n];

        for ( i = 0; i < n; i++ )
        {
            y[i] = x[i];
        }

        for ( i = 0; i < n; i++ )
        {
            x[i] = 0.0;
        }

        int ilo = Math.Max ( 0, shift );
        int ihi = Math.Min ( n, n + shift );

        for ( i = ilo; i < ihi; i++ )
        {
            x[i] = y[i-shift];
        }
    }

    public static void r8vec_shift_circular ( int shift, int n, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SHIFT_CIRCULAR performs a circular shift on an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8 values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int SHIFT, the amount by which each entry is to
        //    be shifted.
        //
        //    Input, int N, the length of the vector.
        //
        //    Input/output, double X[N], the vector to be shifted.
        //
    {
        int i;

        double[] y = new double[n];

        for ( i = 0; i < n; i++ )
        {
            y[i] = x[i];
        }

        for ( i = 0; i < n; i++ )
        {
            int j = i4_wrap ( i - shift, 0, n - 1 );
            x[i] = y[j];
        }
    }

}