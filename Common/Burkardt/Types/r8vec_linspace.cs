namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8vec_linspace ( int n, double a_first, double a_last, ref double[] a, int index = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_LINSPACE creates a vector of linearly spaced values.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
        //
        //    In other words, the interval is divided into N-1 even subintervals,
        //    and the endpoints of intervals are used as the points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double A_FIRST, A_LAST, the first and last entries.
        //
        //    Output, double A[N], a vector of linearly spaced data.
        //
    {
        switch (n)
        {
            case 1:
                a[index + 0] = ( a_first + a_last ) / 2.0;
                break;
            default:
            {
                int i;
                for ( i = 0; i < n; i++ )
                {
                    a[index + i] = ( (n - 1 - i) * a_first 
                                     + i * a_last ) 
                                   / (n - 1);
                }

                break;
            }
        }
    }

    public static double[] r8vec_linspace_new(int n, double a_first, double a_last)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
        //
        //    In other words, the interval is divided into N-1 even subintervals,
        //    and the endpoints of intervals are used as the points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double A_FIRST, A_LAST, the first and last entries.
        //
        //    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
        //
    {
        double[] a = new double[n];

        switch (n)
        {
            case 1:
                a[0] = (a_first + a_last) / 2.0;
                break;
            default:
            {
                for (int i = 0; i < n; i++)
                {
                    a[i] = ((n - 1 - i) * a_first
                            + i * a_last)
                           / (n - 1);
                }

                break;
            }
        }

        return a;
    }

    public static void r8vec_linspace2 ( int n, double a_first, double a_last, ref double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_LINSPACE creates a vector of linearly spaced values.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    4 points evenly spaced between 0 and 12 will yield 2, 4, 6, 8, 10.
        //
        //    In other words, the interval is divided into N+1 even subintervals,
        //    and the endpoints of internal intervals are used as the points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double A_FIRST, A_LAST, the first and last entries.
        //
        //    Output, double A[N], a vector of linearly spaced data.
        //
    {
        switch (n)
        {
            case 1:
                a[0] = ( a_first + a_last ) / 2.0;
                break;
            default:
            {
                int i;
                for ( i = 0; i < n; i++ )
                {
                    a[i] = (( n - i     ) * a_first
                            + (     i + 1 ) * a_last) 
                           / ( n     + 1 );
                }

                break;
            }
        }
    }

}