﻿using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8vec_nint ( int n, ref double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_NINT rounds the entries of an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input/output, double A[N], the vector to be rounded.
        //
    {
        int i;

        for ( i = 0; i < n; i++ )
        {
            int s = a[i] switch
            {
                < 0.0 => -1,
                _ => 1
            };

            a[i] = s * ( int ) ( Math.Abs ( a[i] ) + 0.5 );
        }
    }

    public static double[] r8vec_nint_new ( int n, double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_NINT_NEW rounds the entries of an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, double A[N], the vector to be rounded.
        //
        //    Output, double R8VEC_NINT_NEW[N], the rounded values.
        //
    {
        int i;

        double[] b = new double[n];

        for ( i = 0; i < n; i++ )
        {
            int s = a[i] switch
            {
                < 0.0 => -1,
                _ => 1
            };

            b[i] = s * ( int ) ( Math.Abs ( a[i] ) + 0.5 );
        }

        return b;
    }
}