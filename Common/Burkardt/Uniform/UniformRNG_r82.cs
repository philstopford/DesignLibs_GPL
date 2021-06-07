namespace Burkardt.Uniform
{
    public static partial class UniformRNG
    {
        public static void r82_uniform(double[] rlo, double[] rhi, ref int seed, ref double[] r )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R82_UNIFORM returns a random R82 value in a given range.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double RLO[2], RHI[2], the minimum and maximum values.
            //
            //    Input/output, int *SEED, a seed for the random number generator.
            //
            //    Output, double R[2], the randomly chosen value.
            //
        {
            for (int i = 0; i < 2; i++)
            {
                r[i] = rlo[i] + (rhi[i] - rlo[i]) * r8_uniform_01(ref seed);
            }
        }
        
        public static void r82vec_uniform ( int n, double[] alo, double[] ahi, ref int seed,
        ref double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82VEC_UNIFORM returns a random R82VEC in a given range.
        //
        //  Discussion:
        //
        //    A is a two dimensional array of order N by 2, stored as a vector
        //    of rows: A(0,0), A(0,1), // A(1,0), A(1,1) // ...
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double ALO[2], AHI[2], the range allowed for the entries.
        //
        //    Input/output, int *SEED, a seed for the random number generator.
        //
        //    Output, double A[N*2], the vector of randomly chosen integers.
        //
        {
            for ( int i = 0; i < n; i++ )
            {
                for ( int j = 0; j < 2; j++ )
                {
                    a[2*i+j] = r8_uniform ( alo[j], ahi[j], ref seed );
                }
            }
        }
        
    }
}