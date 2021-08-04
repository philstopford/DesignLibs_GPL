using System;
using Burkardt.Types;

namespace Burkardt.CorrelationNS
{
    public static partial class SamplePaths
    {
        public static double[] sample_paths_cholesky(int n, int n2, double rhomax, double rho0,
            Func< int, double[], double, double[] > correlation, ref typeMethods.r8vecNormalData data, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SAMPLE_PATHS_CHOLESKY: sample paths for stationary correlation functions.
        //
        //  Discussion:
        //
        //    This method uses the Cholesky factorization of the correlation matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points on each path.
        //
        //    Input, int N2, the number of paths.
        //
        //    Input, double RHOMAX, the maximum value of RHO.
        //
        //    Input, double RHO0, the correlation length.
        //
        //    Input, double *CORRELATION ( int n, double rho_vec[], double rho0), 
        //    the name of the function which evaluates the correlation.
        //
        //    Input/output, int &SEED, a seed for the random number
        //    generator.
        //
        //    Output, double X[N*N2], the sample paths.
        //
        {
            double[] cor;
            double[] cor_vec;
            int flag = 0;
            int i;
            int j;
            int k;
            double[] l;
            double[] r;
            double[] rho_vec;
            double rhomin;
            double[] x;
            //
            //  Choose N equally spaced sample points from 0 to RHOMAX.
            //
            rhomin = 0.0;
            rho_vec = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            //
            //  Evaluate the correlation function.
            //
            cor_vec = correlation(n, rho_vec, rho0);
            //
            //  Construct the correlation matrix;
            //
            //  From the vector 
            //    [ C(0), C(1), C(2), ... C(N-1) ]
            //  construct the vector
            //    [ C(N-1), ..., C(2), C(1), C(0), C(1), C(2), ...  C(N-1) ]
            //  Every row of the correlation matrix can be constructed by a subvector
            //  of this vector.
            //
            cor = new double[n * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    k = typeMethods.i4_wrap(j - i, 0, n - 1);
                    cor[i + j * n] = cor_vec[k];
                }
            }

            //
            //  Get the Cholesky factorization of COR:
            //
            //    COR = L * L'.
            //
            l = typeMethods.r8mat_cholesky_factor(n, cor, ref flag);
            //
            //  The matrix might not be nonnegative definite.
            //
            if (flag == 2)
            {
                Console.WriteLine("");
                Console.WriteLine("SAMPLE_PATHS_CHOLESKY - Fatal error!");
                Console.WriteLine("  The correlation matrix is not");
                Console.WriteLine("  symmetric nonnegative definite.");
                return null;
            }

            //
            //  Compute a matrix of N by N2 normally distributed values.
            //
            r = typeMethods.r8mat_normal_01_new(n, n2, ref data, ref seed);
            //
            //  Compute the sample path.
            //
            x = typeMethods.r8mat_mm_new(n, n, n2, l, r);

            return x;
        }

        public static double[] sample_paths2_cholesky(int n, int n2, double rhomin, double rhomax,
                double rho0, Func<int, int, double[], double[], double, double[]> correlation2, ref typeMethods.r8vecNormalData data, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SAMPLE_PATHS2_CHOLESKY: sample paths for stationary correlation functions.
            //
            //  Discussion:
            //
            //    This method uses the Cholesky factorization of the correlation matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 November 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points on each path.
            //
            //    Input, int N2, the number of paths.
            //
            //    Input, double RHOMIN, RHOMAX, the range of RHO.
            //
            //    Input, double RHO0, the correlation length.
            //
            //    Input, double *CORRELATION2 ( int m, int n, double s[], double t[], 
            //    double rho0 ), the name of the function which evaluates the correlation.
            //
            //    Input/output, int &SEED, a seed for the random number
            //    generator.
            //
            //    Output, double X[N*N2], the sample paths.
            //
        {
            double[] cor;
            int flag = 0;
            double[] l;
            double[] r;
            double[] s;
            double[] x;
            //
            //  Choose N equally spaced sample points from RHOMIN to RHOMAX.
            //
            s = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            //
            //  Evaluate the correlation function.
            //
            cor = correlation2(n, n, s, s, rho0);
            //
            //  Get the Cholesky factorization of COR:
            //
            //    COR = L * L'.
            //
            l = typeMethods.r8mat_cholesky_factor(n, cor, ref flag);
            //
            //  The matrix might not be nonnegative definite.
            //
            if (flag == 2)
            {
                Console.WriteLine("");
                Console.WriteLine("SAMPLE_PATHS2_CHOLESKY - Fatal error!");
                Console.WriteLine("  The correlation matrix is not");
                Console.WriteLine("  symmetric nonnegative definite.");
                return null;
            }

            //
            //  Compute a matrix of N by N2 normally distributed values.
            //
            r = typeMethods.r8mat_normal_01_new(n, n2, ref data, ref seed);
            //
            //  Compute the sample path.
            //
            x = typeMethods.r8mat_mm_new(n, n, n2, l, r);

            return x;
        }


    }
}