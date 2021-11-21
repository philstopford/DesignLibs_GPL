using System;

namespace Burkardt.AppliedStatistics;

public static partial class Algorithms
{
    public static double[] wshrt(double[] d, int n, int np, ref int seed )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WSHRT returns a random Wishart variate.
        //
        //  Discussion:
        //
        //    This routine is a Wishart variate generator.  
        //
        //    On output, SA is an upper-triangular matrix of size NP * NP,
        //    written in linear form, column ordered, whose elements have a 
        //    Wishart(N, SIGMA) distribution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 April 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by William Smith, Ronald Hocking.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    William Smith, Ronald Hocking,
        //    Algorithm AS 53, Wishart Variate Generator,
        //    Applied Statistics,
        //    Volume 21, Number 3, pages 341-345, 1972.
        //
        //  Parameters:
        //
        //    Input, double D[NP*(NP+1)/2], the upper triangular array that
        //    represents the Cholesky factor of the correlation matrix SIGMA.
        //    D is stored in column-major form.
        //
        //    Input, int N, the number of degrees of freedom.
        //    1 <= N <= NP.
        //
        //    Input, int NP, the size of variables.
        //
        //    Input/output, int &SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double WSHART[NP*(NP+1)/2], a sample from the 
        //    Wishart distribution.
        //
    {
        double u1 = 0;
        double u2 = 0;

        int k = 0;
        int nnp = np * (np + 1) / 2;
        //
        //  Load SB with independent normal (0, 1) variates.
        //
        double[] sb = new double[nnp];

        while (k < nnp)
        {
            rnorm(ref seed, ref u1, ref u2);

            sb[k] = u1;
            k += 1;

            if (k >= nnp)
            {
                continue;
            }

            sb[k] = u2;
            k += 1;
        }

        //
        //  Load diagonal elements with square root of chi-square variates.
        //
        int ns = 0;

        for (int i = 1; i <= np; i++)
        {
            double df = np - i + 1;
            ns += i;
            u1 = 2.0 / (9.0 * df);
            u2 = 1.0 - u1;
            u1 = Math.Sqrt(u1);
            //
            //  Wilson-Hilferty formula for approximating chi-square variates:
            //  The original code did not take the absolute value!
            //
            sb[ns - 1] = Math.Sqrt(df * Math.Abs(Math.Pow(u2 + sb[ns - 1] * u1, 3)));
        }

        double[] sa = new double[nnp];

        double rn = n;
        int nr = 1;

        for (int i = 1; i <= np; i++)
        {
            nr = nr + i - 1;
            for (int j = i; j <= np; j++)
            {
                int ip = nr;
                int nq = j * (j - 1) / 2 + i - 1;
                double c = 0.0;
                for (k = i; k <= j; k++)
                {
                    ip = ip + k - 1;
                    nq += 1;
                    c += sb[ip - 1] * d[nq - 1];
                }

                sa[ip - 1] = c;
            }
        }

        for (int i = 1; i <= np; i++)
        {
            int ii = np - i + 1;
            int nq = nnp - np;
            for (int j = 1; j <= i; j++)
            {
                int ip = ii * (ii - 1) / 2;
                double c = 0.0;
                for (k = i; k <= np; k++)
                {
                    ip += 1;
                    nq += 1;
                    c += sa[ip - 1] * sa[nq - 1];
                }

                sa[nq - 1] = c / rn;
                nq = nq - 2 * np + i + j - 1;
            }
        }

        return sa;
    }



}