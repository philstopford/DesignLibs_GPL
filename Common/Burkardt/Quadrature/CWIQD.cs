using System;

namespace Burkardt.Quadrature;

public static class CWIQD
{
    public static double[] cwiqd(int m, int nm, int l, double v, double[] xk, int nstar,
            double[] phi, double[] a, double[] r )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CWIQD computes all the weights for a given knot.
        //
        //  Discussion:
        //
        //    The variable names correspond to the 1982 reference, and explanations of
        //    some of the terminology may be found there.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 January 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Sylvan Elhay, Jaroslav Kautsky,
        //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
        //    Interpolatory Quadrature,
        //    ACM Transactions on Mathematical Software,
        //    Volume 13, Number 4, December 1987, pages 399-415.
        //
        //    Jaroslav Kautsky, Sylvan Elhay,
        //    Calculation of the Weights of Interpolatory Quadratures,
        //    Numerische Mathematik,
        //    Volume 40, 1982, pages 407-422.
        //
        //  Parameters:
        //
        //    Input, int M, the multiplicity of the knot in question.
        //
        //    Input, int NM, is equal to max ( N - M, 1 ), where N is 
        //    the number of knots used, counted according to multiplicity.
        //
        //    Input, int L, min ( M, N - M + 1), where N is the number
        //    of knots used, counted according to multiplicity.
        //
        //    Input, double V,  the knot in question.
        //
        //    Input, double XK[NM], all but the last M entries in the 
        //    diagonal of K-hat.
        //
        //    Input, int NSTAR, the dimension of the Jacobi matrix.
        //
        //    Input, double PHI[NSTAR], the eigenvalues of the Jacobi matrix.
        //
        //    Input, double A[NSTAR], the square of the first row of the 
        //    orthogonal matrix that diagonalizes the Jacobi matrix.
        //
        //    Input, double R[L], used to compute the right 
        //    principal vectors.
        //
        //    Output, double CWIQD[M], the weights.
        //
    {
        int i;
        int j;
        int k;
        double sum;

        double[] d = new double[m];
        double[] wf = new double[nstar];
        double[] y = new double[m];
        double[] z = new double[m];
        //
        //  Compute products required for Y-hat.
        //
        for (j = 0; j < nstar; j++)
        {
            wf[j] = a[j];
            for (i = 0; i < nm; i++)
            {
                wf[j] *= (phi[j] - xk[i]);
            }
        }

        //
        //  Compute Y-hat.
        //
        for (i = 0; i < m; i++)
        {
            sum = 0.0;
            for (j = 0; j < nstar; j++)
            {
                sum += wf[j];
                wf[j] *= (phi[j] - v);
            }

            y[i] = sum;
        }

        //
        //  If N = 1 the right principal vector is already in R.
        //  Otherwise compute the R-principal vector of grade M-1.
        //
        for (i = 1; i <= nm; i++)
        {
            double tmp = v - xk[i - 1];

            int last = Math.Min(l, i + 1);
            int jr;
            for (jr = 2; jr <= last; jr++)
            {
                j = last - jr + 2;
                r[j - 1] = tmp * r[j - 1] + r[j - 2];
            }

            r[0] = tmp * r[0];
        }

        //
        //  Compute left principal vector(s) and weight for highest derivative.
        //  The following statement contains the only division in this
        //  routine.  Any test for overflow should be made after it.
        //
        d[m - 1] = y[m - 1] / r[0];

        switch (m)
        {
            case 1:
                return d;
        }

        //
        //  Compute left principal vector.
        //
        z[0] = 1.0 / r[0];
        for (i = 2; i <= m; i++)
        {
            sum = 0.0;
            int minil = Math.Min(i, l);
            for (j = 2; j <= minil; j++)
            {
                k = i - j + 1;
                sum += r[j - 1] * z[k - 1];
            }

            z[i - 1] = -sum * z[0];
        }

        //
        //  Accumulate weights.
        //
        for (i = 2; i <= m; i++)
        {
            sum = 0.0;
            for (j = 1; j <= i; j++)
            {
                k = m - i + j;
                sum += z[j - 1] * y[k - 1];
            }

            k = m - i + 1;
            d[k - 1] = sum;
        }

        return d;
    }
}