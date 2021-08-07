namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8vec_cum ( int n, double[] a, ref double[] a_cum )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_CUM computes the cumulutive sums of an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    Input:
        //
        //      A = { 1.0, 2.0, 3.0, 4.0 }
        //
        //    Output:
        //
        //      A_CUM = { 1.0, 3.0, 6.0, 10.0 }
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double A[N], the vector to be summed.
        //
        //    Output, double A_CUM[N], the cumulative sums.
        //
        {
            int i;

            a_cum[0] = a[0];

            for ( i = 1; i < n; i++ )
            {
                a_cum[i] = a_cum[i-1] + a[i];
            }
        }
        
        public static double[] r8vec_cum_new(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_CUM_NEW computes the cumulutive sums of an R8VEC.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    Input:
            //
            //      A = { 1.0, 2.0, 3.0, 4.0 }
            //
            //    Output:
            //
            //      A_CUM = { 1.0, 3.0, 6.0, 10.0 }
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 May 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, double A[N], the vector to be summed.
            //
            //    Output, double R8VEC_CUM_NEW[N], the cumulative sums.
            //
        {
            double[] a_cum;
            int i;

            a_cum = new double[n];

            a_cum[0] = a[0];

            for (i = 1; i < n; i++)
            {
                a_cum[i] = a_cum[i - 1] + a[i];
            }

            return a_cum;
        }

        public static double[] r8vec_cum0_new(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_CUM0_NEW computes the cumulutive sums of an R8VEC.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    Input:
            //
            //      A = { 1.0, 2.0, 3.0, 4.0 }
            //
            //    Output:
            //
            //      A_CUM = { 0.0, 1.0, 3.0, 6.0, 10.0 }
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 May 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, double A[N], the vector to be summed.
            //
            //    Output, double R8VEC_CUM0_NEW[N+1], the cumulative sums.
            //
        {
            double[] a_cum;
            int i;

            a_cum = new double[n + 1];

            a_cum[0] = 0.0;

            for (i = 1; i <= n; i++)
            {
                a_cum[i] = a_cum[i - 1] + a[i - 1];
            }

            return a_cum;
        }
    }
}