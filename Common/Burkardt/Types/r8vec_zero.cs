namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r8vec_zero_new(int n)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_ZERO_NEW creates and zeroes an R8VEC.
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
            //    10 July 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Output, double R8VEC_ZERO_NEW[N], a vector of zeroes.
            //
        {
            double[] a = new double[n];

            for (int i = 0; i < n; i++)
            {
                a[i] = 0.0;
            }

            return a;
        }

        public static void r8vec_zeros ( int n, ref double[] a )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_ZEROS zeroes an R8VEC.
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
            //    Input, int N, the number of entries in the vector.
            //
            //    Output, double A[N], a vector of zeroes.
            //
        {
            int i;

            for ( i = 0; i < n; i++ )
            {
                a[i] = 0.0;
            }
            return;
        }

        public static double[] r8vec_zeros_new ( int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_ZEROS_NEW creates and zeroes an R8VEC.
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
            //    10 July 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Output, double R8VEC_ZEROS_NEW[N], a vector of zeroes.
            //
        {
            double []a;
            int i;

            a = new double[n];

            for ( i = 0; i < n; i++ )
            {
                a[i] = 0.0;
            }
            return a;
        }
    }
}