namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static double eval_pol ( double[] a, int n, double x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EVAL_POL evaluates a polynomial at X.
        //
        //  Discussion:
        //
        //    EVAL_POL = A(0) + A(1)*X + ... + A(N)*X^N
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2021
        //
        //  Author:
        //
        //    Barry Brown, James Lovato, Kathy Russell.
        //
        //  Input:
        //
        //    double A(0:N), coefficients of the polynomial.
        //
        //    int *N, length of A.
        //
        //    double *X, the point at which the polynomial
        //    is to be evaluated.
        //
        //  Output:
        //
        //    double EVAL_POL, the value of the polynomial at X.
        //
    {
        int i;

        double term = a[n-1];
        for ( i = n-1-1; i >= 0; i-- )
        {
            term = a[i]+term*x;
        }

        double devlpl = term;
        return devlpl;
    }
}