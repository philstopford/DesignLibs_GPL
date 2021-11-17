namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8vec_concatenate ( int n1, double[] a, int n2, double[] b, ref double[] c )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_CONCATENATE concatenates two R8VEC's.
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
        //    22 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, the number of entries in the first vector.
        //
        //    Input, double A[N1], the first vector.
        //
        //    Input, int N2, the number of entries in the second vector.
        //
        //    Input, double B[N2], the second vector.
        //
        //    Output, double C[N1+N2], the concatenated vector.
        //
    {
        int i;

        for ( i = 0; i < n1; i++ )
        {
            c[i] = a[i];
        }
        for ( i = 0; i < n2; i++ )
        {
            c[n1+i] = b[i];
        }
    }

    public static double[] r8vec_concatenate_new ( int n1, double[] a, int n2, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_CONCATENATE_NEW concatenates two R8VEC's.
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
        //    22 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, the number of entries in the first vector.
        //
        //    Input, double A[N1], the first vector.
        //
        //    Input, int N2, the number of entries in the second vector.
        //
        //    Input, double B[N2], the second vector.
        //
        //    Output, double R8VEC_CONCATENATE_NEW[N1+N2], the concatenated vector.
        //
    {
        int i;

        double[] c = new double[n1+n2];

        for ( i = 0; i < n1; i++ )
        {
            c[i] = a[i];
        }
        for ( i = 0; i < n2; i++ )
        {
            c[n1+i] = b[i];
        }

        return c;
    }
}