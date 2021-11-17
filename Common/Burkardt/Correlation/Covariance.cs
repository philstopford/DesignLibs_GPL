using System;
using Burkardt.Types;

namespace Burkardt.CorrelationNS;

public static partial class Correlation
{
    public static double[] correlation_to_covariance(int n, double[] c, double[] sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CORRELATION_TO_COVARIANCE: covariance matrix from a correlation matrix.
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
        //    Input, int N, the order of the matrix.
        //
        //    Input, double C[N*N], the correlation matrix.
        //
        //    Input, double SIGMA[N], the standard deviations.
        //
        //    Output, double K[N*N], the covariance matrix.
        //
    {
        double c_max;
        double c_min;
        double e;
        double error_frobenius;
        int i;
        int j;
        double[] k;
        double tol;

        tol = Math.Sqrt(typeMethods.r8_epsilon());
        //
        //  C must be symmetric.
        //
        error_frobenius = typeMethods.r8mat_is_symmetric(n, n, c);

        if (tol < error_frobenius)
        {
            Console.WriteLine("");
            Console.WriteLine("CORRELATION_TO_COVARIANCE - Fatal error!");
            Console.WriteLine("  Input matrix C fails symmetry test with error " + error_frobenius + "");
            return null;
        }

        //
        //  The diagonal must be 1.
        //
        for (i = 0; i < n; i++)
        {
            e = Math.Abs(c[i + i * n] - 1.0);
            if (tol < e)
            {
                Console.WriteLine("");
                Console.WriteLine("CORRELATION_TO_COVARIANCE - Fatal error!");
                Console.WriteLine("  Input matrix C has non-unit diagonal entries.");
                Console.WriteLine("  Error on row " + i + " is " + e + "");
                return null;
            }
        }

        //
        //  Off-diagonals must be between -1 and 1.
        //
        c_min = typeMethods.r8mat_max(n, n, c);

        if (c_min < -1.0 - tol)
        {
            Console.WriteLine("");
            Console.WriteLine("CORRELATION_TO_COVARIANCE - Fatal error!");
            Console.WriteLine("  Input matrix C has entries less than -1.0");
            return null;
        }

        c_max = typeMethods.r8mat_max(n, n, c);

        if (1.0 + tol < c_max)
        {
            Console.WriteLine("");
            Console.WriteLine("CORRELATION_TO_COVARIANCE - Fatal error!");
            Console.WriteLine("  Input matrix C has entries greater than +1.0");
            return null;
        }

        //
        //  Form K.
        //
        k = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                k[i + j * n] = sigma[i] * c[i + j * n] * sigma[j];
            }
        }

        return k;
    }

    public static void covariance_to_correlation(int n, double[] k, ref double[] c, ref double[] sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COVARIANCE_TO_CORRELATION: correlation matrix from a covariance matrix.
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
        //    Input, int N, the order of the matrix.
        //
        //    Input, double K[N*N], the covariance matrix.
        //
        //    Output, double C[N*N], the correlation matrix.
        //
        //    Output, double SIGMA[N], the standard deviations.
        //
    {
        double e;
        double error_frobenius;
        int i;
        int j;
        double sigma_min;
        double tol;

        tol = Math.Sqrt(typeMethods.r8_epsilon());
        //
        //  K must be symmetric.
        //
        error_frobenius = typeMethods.r8mat_is_symmetric(n, n, k);

        if (tol < error_frobenius)
        {
            Console.WriteLine("");
            Console.WriteLine("COVARIANCE_TO_CORRELATION - Fatal error");
            Console.WriteLine("  Input matrix K fails symmetry test with error " + error_frobenius + "");
            return;
        }

        //
        //  It must be the case that K(I,J)^2 <= K(I,I) * K(J,J).
        //
        e = 0.0;
        for (i = 0; i < n; i++)
        {
            for (j = i + 1; j < n; j++)
            {
                e = Math.Max(e, k[i + j * n] * k[i + j * n] - k[i + i * n] * k[j + j * n]);
            }
        }

        if (tol < e)
        {
            Console.WriteLine("");
            Console.WriteLine("COVARIANCE_TO_CORRELATION - Fatal error");
            Console.WriteLine("  Input matrix K fails K(I,J)^2 <= K(I,I)*K(J,J)");
            return;
        }

        //
        //  Get the diagonal.
        //
        for (i = 0; i < n; i++)
        {
            sigma[i] = k[i + i * n];
        }

        //
        //  Ensure the diagonal is positive.
        //
        sigma_min = typeMethods.r8vec_min(n, sigma);

        switch (sigma_min)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("COVARIANCE_TO_CORRELATION - Fatal error!");
                Console.WriteLine("  Input matrix K has nonpositive diagonal entry = " + sigma_min + "");
                return;
        }

        //
        //  Convert from variance to standard deviation.
        //
        for (i = 0; i < n; i++)
        {
            sigma[i] = Math.Sqrt(sigma[i]);
        }

        //
        //  Form C.
        //
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                c[i + j * n] = k[i + j * n] / sigma[i] / sigma[j];
            }
        }
    }
}