using System;
using Burkardt.FullertonFnLib;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace Burkardt.CorrelationNS;

public static partial class SamplePaths
{
    public static Correlation.CorrelationResult sample_paths_eigen(FullertonLib.BesselData globaldata, FullertonLib.r8BESJ0Data jdata, int n, int n2, double rhomax, double rho0,
            Func<FullertonLib.BesselData, FullertonLib.r8BESJ0Data, int, double[], double, Correlation.CorrelationResult> correlation, ref typeMethods.r8vecNormalData data, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SAMPLE_PATHS_EIGEN: sample paths for stationary correlation functions.
        //
        //  Discussion:
        //
        //    This method uses the eigen-decomposition of the correlation matrix.
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
        int i;
        int j;
        int k;
        //
        //  Choose N equally spaced sample points from 0 to RHOMAX.
        //
        double rhomin = 0.0;
        double[] rho_vec = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
        //
        //  Evaluate the correlation function.
        //
        Correlation.CorrelationResult tr = correlation(globaldata, jdata, n, rho_vec, rho0);
        double[] cor_vec = tr.result;
        globaldata = tr.data;
        jdata = tr.j0data;
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
        double[] cor = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                k = typeMethods.i4_wrap(Math.Abs(i - j), 0, n - 1);
                cor[i + j * n] = cor_vec[k];
            }
        }

        //
        //  Get the eigendecomposition of COR:
        //
        //    COR = V * D * V'.
        //
        //  Because COR is symmetric, V is orthogonal.
        //
        double[] d = new double[n];
        double[] w = new double[n];
        double[] v = new double[n * n];

        TRED2.tred2(n, cor, ref d, ref w, ref v);

        TQL2.tql2(n, ref d, ref w, v);
        //
        //  We assume COR is non-negative definite, and hence that there
        //  are no negative eigenvalues.
        //
        double dmin = typeMethods.r8vec_min(n, d);

        if (dmin < -Math.Sqrt(typeMethods.r8_epsilon()))
        {
            Console.WriteLine("");
            Console.WriteLine("SAMPLE_PATHS_EIGEN - Warning!");
            Console.WriteLine("  Negative eigenvalues observed as low as " + dmin + "");
        }

        for (i = 0; i < n; i++)
        {
            d[i] = Math.Max(d[i], 0.0);
        }

        //
        //  Compute the eigenvalues of the factor C.
        //
        for (i = 0; i < n; i++)
        {
            d[i] = Math.Sqrt(d[i]);
        }

        //
        //  Compute C, such that C' * C = COR.
        //
        double[] c = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                c[i + j * n] = 0.0;
                for (k = 0; k < n; k++)
                {
                    c[i + j * n] += d[k] * v[i + k * n] * v[j + k * n];
                }
            }
        }

        //
        //  Compute N by N2 independent random normal values.
        //
        double[] r = typeMethods.r8mat_normal_01_new(n, n2, ref data, ref seed);
        //
        //  Multiply to get the variables X which have correlation COR.
        //
        double[] x = typeMethods.r8mat_mm_new(n, n, n2, c, r);

        Correlation.CorrelationResult result = new() { result = x, data = globaldata };
            
        return result;
    }

    public static double[] sample_paths2_eigen(int n, int n2, double rhomin, double rhomax,
            double rho0, Func<int, int, double[], double[], double, double[]> correlation2, ref typeMethods.r8vecNormalData data, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SAMPLE_PATHS2_EIGEN: sample paths for stationary correlation functions.
        //
        //  Discussion:
        //
        //    This method uses the eigen-decomposition of the correlation matrix.
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
        int i;
        int j;
        //
        //  Choose N equally spaced sample points from RHOMIN to RHOMAX.
        //
        double[] s = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
        //
        //  Evaluate the correlation function.
        //
        double[] cor = correlation2(n, n, s, s, rho0);
        //
        //  Get the eigendecomposition of COR:
        //
        //    COR = V * D * V'.
        //
        //  Because COR is symmetric, V is orthogonal.
        //
        double[] d = new double[n];
        double[] w = new double[n];
        double[] v = new double[n * n];

        TRED2.tred2(n, cor, ref d, ref w, ref v);

        TQL2.tql2(n, ref d, ref w, v);
        //
        //  We assume COR is non-negative definite, and hence that there
        //  are no negative eigenvalues.
        //
        double dmin = typeMethods.r8vec_min(n, d);

        if (dmin < -Math.Sqrt(typeMethods.r8_epsilon()))
        {
            Console.WriteLine("");
            Console.WriteLine("SAMPLE_PATHS2_EIGEN - Warning!");
            Console.WriteLine("  Negative eigenvalues observed as low as " + dmin + "");
        }

        for (i = 0; i < n; i++)
        {
            d[i] = Math.Max(d[i], 0.0);
        }

        //
        //  Compute the eigenvalues of the factor C.
        //
        for (i = 0; i < n; i++)
        {
            d[i] = Math.Sqrt(d[i]);
        }

        //
        //  Compute C, such that C' * C = COR.
        //
        double[] c = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                c[i + j * n] = 0.0;
                int k;
                for (k = 0; k < n; k++)
                {
                    c[i + j * n] += d[k] * v[i + k * n] * v[j + k * n];
                }
            }
        }

        //
        //  Compute N by N2 independent random normal values.
        //
        double[] r = typeMethods.r8mat_normal_01_new(n, n2, ref data, ref seed);
        //
        //  Multiply to get the variables X which have correlation COR.
        //
        double[] x = typeMethods.r8mat_mm_new(n, n, n2, c, r);

        return x;
    }

}