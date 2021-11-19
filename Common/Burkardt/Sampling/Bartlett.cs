using System;
using Burkardt.PDFLib;
using Burkardt.Types;

namespace Burkardt.Sampling;

public static class Bartlett
{
    public static double[] bartlett_sample(int m, int df, double[] sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BARTLETT_SAMPLE samples the Bartlett distribution.
        //
        //  Discussion:
        //
        //    If the matrix T is sampled from the Bartlett distribution, then 
        //    the matrix W = T' * T is a sample from the Wishart distribution.
        //
        //    This function requires functions from the PDFLIB and RNGLIB libraries.
        //
        //    The "initialize()" function from RNGLIB must be called before using
        //    this function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Patrick Odell, Alan Feiveson,
        //    A numerical procedure to generate a sample covariance matrix,
        //    Journal of the American Statistical Association,
        //    Volume 61, Number 313, March 1966, pages 199-203.
        //
        //    Stanley Sawyer,
        //    Wishart Distributions and Inverse-Wishart Sampling,
        //    Washington University,
        //    30 April 2007, 12 pages.
        //
        //  Parameters:
        //
        //    Input, int M, the order of the matrix.
        //
        //    Input, int DF, the number of degrees of freedom.
        //    M <= DF.
        //
        //    Input, double SIGMA[M*M], the covariance matrix, which should be 
        //    a symmetric positive definite matrix.
        //
        //    Output, double BARTLETT_SAMPLE[M*M], the sample matrix from 
        //    the Bartlett distribution.
        //
    {
        int flag = 0;

        if (df < m)
        {
            Console.WriteLine("");
            Console.WriteLine("BARTLETT_SAMPLE - Fatal error!");
            Console.WriteLine("  DF = " + df + " < M = " + m + ".");
            return null;
        }

        //
        //  Get the upper triangular Cholesky factor of SIGMA.
        //
        double[] r = typeMethods.r8mat_cholesky_factor_upper(m, sigma, ref flag);

        if (flag != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("BARTLETT_SAMPLE - Fatal error!");
            Console.WriteLine("  Unexpected error return from R8MAT_CHOLESKY_FACTOR_UPPER.");
            Console.WriteLine("  FLAG = " + flag + "");
            return null;
        }

        //
        //  Sample the unit Bartlett distribution.
        //
        double[] tu = bartlett_unit_sample(m, df);
        //
        //  Construct the matrix T = TU * R.
        //
        double[] t = typeMethods.r8mat_mm_new(m, m, m, tu, r);

        return t;
    }

    public static double[] bartlett_unit_sample(int m, int df)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BARTLETT_UNIT_SAMPLE samples the unit Bartlett distribution.
        //
        //  Discussion:
        //
        //    If the matrix T is sampled from the unit Bartlett distribution, then 
        //    the matrix W = T' * T is a sample from the unit Wishart distribution.
        //
        //    This function requires functions from the PDFLIB and RNGLIB libraries.
        //
        //    The "initialize()" function from RNGLIB must be called before using
        //    this function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Patrick Odell, Alan Feiveson,
        //    A numerical procedure to generate a sample covariance matrix,
        //    Journal of the American Statistical Association,
        //    Volume 61, Number 313, March 1966, pages 199-203.
        //
        //    Stanley Sawyer,
        //    Wishart Distributions and Inverse-Wishart Sampling,
        //    Washington University,
        //    30 April 2007, 12 pages.
        //
        //  Parameters:
        //
        //    Input, int M, the order of the matrix.
        //
        //    Input, int DF, the number of degrees of freedom.
        //    M <= DF.
        //
        //    Output, double BARTLETT_UNIT_SAMPLE[M*M], the sample matrix from the 
        //    unit Bartlett distribution.
        //
    {
        int i;

        if (df < m)
        {
            Console.WriteLine("");
            Console.WriteLine("BARTLETT_UNIT_SAMPLE - Fatal error!");
            Console.WriteLine("  DF = " + df + " < M = " + m + "");
            return null;
        }

        double[] t = new double[m * m];

        for (i = 0; i < m; i++)
        {
            int j;
            for (j = 0; j < i; j++)
            {
                t[i + j * m] = 0.0;
            }

            double df_chi = df - i;
            t[i + i * m] = Math.Sqrt(PDF.r8_chi_sample(df_chi));
            for (j = i + 1; j < m; j++)
            {
                t[i + j * m] = PDF.r8_normal_01_sample();
            }
        }

        return t;
    }
}