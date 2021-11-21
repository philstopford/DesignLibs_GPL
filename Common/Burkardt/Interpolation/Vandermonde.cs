using System;
using Burkardt.SolveNS;
using Burkardt.Types;

namespace Burkardt.Interpolation;

public static class Vandermonde
{
    public static double[] vandermonde_coef_1d(int nd, double[] xd, double[] yd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VANDERMONDE_COEF_1D computes coefficients of a 1D Vandermonde interpolant.
        //
        //  Discussion:
        //
        //    We assume the interpolant has the form
        //
        //      p(x) = c1 + c2 * x + c3 * x^2 + ... + cn * x^(n-1).
        //
        //    We have n data values (x(i),y(i)) which must be interpolated:
        //
        //      p(x(i)) = c1 + c2 * x(i) + c3 * x(i)^2 + ... + cn * x(i)^(n-1) = y(i)
        //
        //    This can be cast as an NxN linear system for the polynomial
        //    coefficients:
        //
        //      [ 1 x1 x1^2 ... x1^(n-1) ] [  c1 ] = [  y1 ]
        //      [ 1 x2 x2^2 ... x2^(n-1) ] [  c2 ] = [  y2 ]
        //      [ ...................... ] [ ... ] = [ ... ]
        //      [ 1 xn xn^2 ... xn^(n-1) ] [  cn ] = [  yn ]
        //
        //    and if the x values are distinct, the system is theoretically
        //    invertible, so we can retrieve the coefficient vector c and
        //    evaluate the interpolant.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ND, the number of data points.
        //
        //    Input, double XD[ND], YD[ND], the data values.
        //
        //    Output, double VANDERMONDE_COEF_1D[ND], the coefficients of the 
        //    interpolating polynomial.
        //
    {
        double[] ad = vandermonde_matrix_1d(nd, xd);

        double[] cd = QRSolve.qr_solve(nd, nd, ad, yd);

        return cd;
    }

    public static double[] vandermonde_matrix_1d(int nd, double[] xd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VANDERMONDE_MATRIX_1D computes a Vandermonde 1D interpolation matrix.
        //
        //  Discussion:
        //
        //    We assume the interpolant has the form
        //
        //      p(x) = c1 + c2 * x + c3 * x^2 + ... + cn * x^(n-1).
        //
        //    We have n data values (x(i),y(i)) which must be interpolated:
        //
        //      p(x(i)) = c1 + c2 * x(i) + c3 * x(i)^2 + ... + cn * x(i)^(n-1) = y(i)
        //
        //    This can be cast as an NxN linear system for the polynomial
        //    coefficients:
        //
        //      [ 1 x1 x1^2 ... x1^(n-1) ] [  c1 ] = [  y1 ]
        //      [ 1 x2 x2^2 ... x2^(n-1) ] [  c2 ] = [  y2 ]
        //      [ ...................... ] [ ... ] = [ ... ]
        //      [ 1 xn xn^2 ... xn^(n-1) ] [  cn ] = [  yn ]
        //
        //    and if the x values are distinct, the matrix A is theoretically
        //    invertible (though in fact, generally badly conditioned).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ND, the number of data points.
        //
        //    Input, double XD[ND], the data values.
        //
        //    Output, double VANDERMONDE_MATRIX_1D[ND*ND], the Vandermonde matrix for X.
        //
    {
        int i;
        int j;

        double[] ad = new double[nd * nd];

        for (i = 0; i < nd; i++)
        {
            ad[i + 0 * nd] = 1.0;
        }

        for (j = 1; j < nd; j++)
        {
            for (i = 0; i < nd; i++)
            {
                ad[i + j * nd] = ad[i + (j - 1) * nd] * xd[i];
            }
        }

        return ad;
    }

    public static double[] vandermonde_value_1d(int nd, double[] cd, int ni, double[] xi)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VANDERMONDE_VALUE_1D evaluates a Vandermonde interpolant.
        //
        //  Discussion:
        //
        //    The polynomial 
        //
        //      p(x) = cd0 + cd1 * x + cd2 * x^2 + ... + cd(nd-1) * x^(nd-1)
        //
        //    is to be evaluated at the vector of values X.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ND, the number of data values.
        //
        //    Input, double CD[ND], the polynomial coefficients.  
        //    CD[I] is the coefficient of X^I.
        //
        //    Input, int NI, the number of interpolation points.
        //
        //    Input, double XI[NI], the interpolation points.
        //
        //    Output, double VANDERMONDE_VALUE_1D[NI], the interpolation values.
        //
    {
        int i;
        int j;

        double[] yi = new double[ni];

        for (j = 0; j < ni; j++)
        {
            yi[j] = cd[nd - 1];
        }

        for (i = nd - 2; 0 <= i; i--)
        {
            for (j = 0; j < ni; j++)
            {
                yi[j] = yi[j] * xi[j] + cd[i];
            }
        }

        return yi;
    }

    public static double[] vandermonde_interp_2d_matrix(int n, int m, double[] x, double[] y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VANDERMONDE_INTERP_2D_MATRIX computes a Vandermonde 2D interpolation matrix.
        //
        //  Discussion:
        //
        //    We assume the approximating function has the form of a polynomial
        //    in X and Y of total degree M.
        //
        //      p(x,y) = c00 
        //             + c10 * x                + c01 * y
        //             + c20 * x^2   + c11 * xy + c02 * y^2
        //             + ...
        //             + cm0 * x^(m) + ...      + c0m * y^m.
        //
        //    If we let T(K) = the K-th triangular number 
        //            = sum ( 1 <= I <= K ) I
        //    then the number of coefficients in the above polynomial is T(M+1).
        //
        //    We have n data locations (x(i),y(i)) and values z(i) to approximate:
        //
        //      p(x(i),y(i)) = z(i)
        //
        //    and we assume that N = T(M+1).
        //
        //    This can be cast as an NxN linear system for the polynomial
        //    coefficients:
        //
        //      [ 1 x1 y1  x1^2 ... y1^m ] [ c00 ] = [  z1 ]
        //      [ 1 x2 y2  x2^2 ... y2^m ] [ c10 ] = [  z2 ]
        //      [ 1 x3 y3  x3^2 ... y3^m ] [ c01 ] = [  z3 ]
        //      [ ...................... ] [ ... ] = [ ... ]
        //      [ 1 xn yn  xn^2 ... yn^m ] [ c0n ] = [  zn ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of data points.  It is necessary 
        //    that N = T(M+1), where T(K) is the K-th triangular number.
        //
        //    Input, int M, the degree of the polynomial.
        //
        //    Input, double X[N], Y[N], the data locations.
        //
        //    Output, double VANDERMONDE_INTERP_2D_MATRIX[N*N], the Vandermonde matrix for X.
        //
    {
        int s;

        int tmp1 = typeMethods.triangle_num(m + 1);

        if (n != tmp1)
        {
            Console.WriteLine("");
            Console.WriteLine("VANDERMONDE_INTERP_2D_MATRIX - Fatal error!");
            Console.WriteLine("  For interpolation, we need N = T(M+1).");
            Console.WriteLine("  But we have N = " + n + "");
            Console.WriteLine("  M = " + m + "");
            Console.WriteLine("  and T(M+1) = " + tmp1 + "");
            return null;
        }

        double[] a = new double[n * n];
        int j = 0;

        for (s = 0; s <= m; s++)
        {
            int ex;
            for (ex = s; 0 <= ex; ex--)
            {
                int ey = s - ex;
                int i;
                for (i = 0; i < n; i++)
                {
                    a[i + j * n] = Math.Pow(x[i], ex) * Math.Pow(y[i], ey);
                }

                j += 1;
            }
        }

        return a;
    }
}