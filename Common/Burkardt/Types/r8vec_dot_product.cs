﻿namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8vec_dot(int n, double[] a1, double[] a2, int a1Index = 0, int a2Index = 0)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_DOT computes the dot product of a pair of R8VEC's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vectors.
        //
        //    Input, double A1[N], A2[N], the two vectors to be considered.
        //
        //    Output, double R8VEC_DOT, the dot product of the vectors.
        //
    {
        double value = 0.0;
        for (int i = 0; i < n; i++)
        {
            value += a1[a1Index + i] * a2[a2Index + i];
        }

        return value;
    }

    public static double r8vec_dot_product(int n, double[] a1, double[] a2, int a1Index = 0, int a2Index = 0)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
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
        //    03 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vectors.
        //
        //    Input, double A1[N], A2[N], the two vectors to be considered.
        //
        //    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
        //
    {
        double value = 0.0;
        for (int i = 0; i < n; i++)
        {
            value += a1[a1Index + i] * a2[a2Index + i];
        }

        return value;
    }
        
    public static double r8vec_dot_product_affine ( int n, double[] v0, double[] v1, double[] v2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_DOT_PRODUCT_AFFINE computes the affine dot product.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vectors.
        //
        //    Input, double V0[N], the base vector.
        //
        //    Input, double V1[N], V2[N], the two vectors to be considered.
        //
        //    Output, double R8VEC_DOT_PRODUCT_AFFINE, the dot product of the vectors.
        //
    {
        int i;

        double value = 0.0;
        for ( i = 0; i < n; i++ )
        {
            value += ( v1[i] - v0[i] ) * ( v2[i] - v0[i] );
        }
        return value;
    }

        
    public static double r8vec_i4vec_dot_product(int n, double[] r8vec, int[] i4vec)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_I4VEC_DOT_PRODUCT computes the dot product of an R8VEC and an I4VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    An I4VEC is a vector of I4's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vectors.
        //
        //    Input, double R8VEC[N], the first vector.
        //
        //    Input, int I4VEC[N], the second vector.
        //
        //    Output, double R8VEC_I4VEC_DOT_PRODUCT, the dot product of the vectors.
        //
    {
        int i;

        double value = 0.0;
        for (i = 0; i < n; i++)
        {
            value += r8vec[i] * i4vec[i];
        }

        return value;
    }

}