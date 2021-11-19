using System;
using Burkardt.PDFLib;
using Burkardt.Types;

namespace Burkardt.Sampling;

public static class Wishart
{
    public static double[] wishart_sample(int m, int df, double[] sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WISHART_SAMPLE samples the Wishart distribution.
        //
        //  Discussion:
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
        //    Output, double WISHART_SAMPLE[M*M], the sample matrix from 
        //    the Wishart distribution.
        //
    {
        int flag = 0;

        if (df < m)
        {
            Console.WriteLine("");
            Console.WriteLine("WISHART_SAMPLE - Fatal error!");
            Console.WriteLine("  DF = " + df + " < M = " + m + "");
            return null;
        }

        //
        //  Get R, the upper triangular Cholesky factor of SIGMA.
        //
        double[] r = typeMethods.r8mat_cholesky_factor_upper(m, sigma, ref flag);

        if (flag != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("WISHART_SAMPLE - Fatal error!");
            Console.WriteLine("  Unexpected error return from R8MAT_CHOLESKY_FACTOR_UPPER.");
            Console.WriteLine("  FLAG = " + flag + "");
            return null;
        }

        //
        //  Get AU, a sample from the unit Wishart distribution.
        //
        double[] au = wishart_unit_sample(m, df);
        //
        //  Construct the matrix A = R' * AU * R.
        //
        double[] aur = typeMethods.r8mat_mm_new(m, m, m, au, r);
        double[] a = typeMethods.r8mat_mtm_new(m, m, m, r, aur);

        return a;
    }

    public static double[] wishart_sample_inverse(int m, int df, double[] sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WISHART_SAMPLE_INVERSE returns the inverse of a sample Wishart matrix.
        //
        //  Discussion:
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
        //    Input, double SIGMA[M*M], the covariance matrix, which should be 
        //    a symmetric positive definite matrix.
        //
        //    Output, double WISHART_SAMPLE[M*M], the inverse of a sample matrix from 
        //    the Wishart distribution.
        //
    {
        int flag = 0;

        if (df < m)
        {
            Console.WriteLine("");
            Console.WriteLine("WISHART_SAMPLE - Fatal error!");
            Console.WriteLine("  DF = " + df + " < M = " + m + "");
            return null;
        }

        //
        //  Get R, the upper triangular Cholesky factor of SIGMA.
        //
        double[] r = typeMethods.r8mat_cholesky_factor_upper(m, sigma, ref flag);

        if (flag != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("WISHART_SAMPLE - Fatal error!");
            Console.WriteLine("  Unexpected error return from R8MAT_CHOLESKY_FACTOR_UPPER.");
            Console.WriteLine("  FLAG = " + flag + "");
            return null;
        }

        //
        //  Get S, the inverse of R.
        //
        double[] s = typeMethods.r8ut_inverse(m, r);
        //
        //  Get UA, the inverse of a sample from the unit Wishart distribution.
        //
        double[] ua = wishart_unit_sample_inverse(m, df);
        //
        //  Construct the matrix A = S * UA * S'.
        //
        double[] uas = typeMethods.r8mat_mmt_new(m, m, m, ua, s);
        double[] a = typeMethods.r8mat_mm_new(m, m, m, s, uas);

        return a;
    }

    public static double[] wishart_unit_sample(int m, int df)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WISHART_UNIT_SAMPLE samples the unit Wishart distribution.
        //
        //  Discussion:
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
        //    Output, double WISHART_UNIT_SAMPLE[M*M], the sample matrix from the 
        //    unit Wishart distribution.
        //
    {
        int i;

        if (df < m)
        {
            Console.WriteLine("");
            Console.WriteLine("WISHART_UNIT_SAMPLE - Fatal error!");
            Console.WriteLine("  DF = " + df + " < M = " + m + ".");
            return null;
        }

        double[] c = new double[m * m];

        for (i = 0; i < m; i++)
        {
            int j;
            for (j = 0; j < i; j++)
            {
                c[i + j * m] = 0.0;
            }

            double df_chi = df - i;
            c[i + i * m] = Math.Sqrt(PDF.r8_chi_sample(df_chi));
            for (j = i + 1; j < m; j++)
            {
                c[i + j * m] = PDF.r8_normal_01_sample();
            }
        }

        double[] a = typeMethods.r8mat_mtm_new(m, m, m, c, c);

        return a;
    }

    public static double[] wishart_unit_sample_inverse(int m, int df)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WISHART_UNIT_SAMPLE_INVERSE inverts a unit Wishart sample matrix.
        //
        //  Discussion:
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
        //    Output, double WISHART_UNIT_SAMPLE[M*M], the inverse of a
        //    sample matrix from the unit Wishart distribution.
        //
    {
        int i;

        if (df < m)
        {
            Console.WriteLine("");
            Console.WriteLine("WISHART_UNIT_SAMPLE_INVERSE - Fatal error!");
            Console.WriteLine("  DF = " + df + " < M = " + m + ".");
            return null;
        }

        double[] c = new double[m * m];

        for (i = 0; i < m; i++)
        {
            int j;
            for (j = 0; j < i; j++)
            {
                c[i + j * m] = 0.0;
            }

            double df_chi = df - i;
            c[i + i * m] = Math.Sqrt(PDF.r8_chi_sample(df_chi));
            for (j = i + 1; j < m; j++)
            {
                c[i + j * m] = PDF.r8_normal_01_sample();
            }
        }

        //
        //  Compute B, the inverse of C.
        //
        double[] b = typeMethods.r8ut_inverse(m, c);
        //
        //  The inverse of the Wishart sample matrix C'*C is inv(C) * C'.
        //
        double[] a = typeMethods.r8mat_mmt_new(m, m, m, b, b);

        return a;
    }
}