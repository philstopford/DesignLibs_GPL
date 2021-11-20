using Burkardt.MatrixNS;
using Burkardt.SolveNS;
using Burkardt.Types;

namespace Burkardt.PolynomialNS;

public static class Vandermonde
{
    public static double[] vandermonde_approx_1d_coef(int n, int m, double[] x, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VANDERMONDE_APPROX_1D_COEF computes a 1D polynomial approximant.
        //
        //  Discussion:
        //
        //    We assume the approximating function has the form
        //
        //      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m.
        //
        //    We have n data values (x(i),y(i)) which must be approximated:
        //
        //      p(x(i)) = c0 + c1 * x(i) + c2 * x(i)^2 + ... + cm * x(i)^m = y(i)
        //
        //    This can be cast as an Nx(M+1) linear system for the polynomial
        //    coefficients:
        //
        //      [ 1 x1 x1^2 ... x1^m ] [  c0 ] = [  y1 ]
        //      [ 1 x2 x2^2 ... x2^m ] [  c1 ] = [  y2 ]
        //      [ .................. ] [ ... ] = [ ... ]
        //      [ 1 xn xn^2 ... xn^m ] [  cm ] = [  yn ]
        //
        //    In the typical case, N is greater than M+1 (we have more data and equations
        //    than degrees of freedom) and so a least squares solution is appropriate,
        //    in which case the computed polynomial will be a least squares approximant
        //    to the data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of data points.
        //
        //    Input, int M, the degree of the polynomial.
        //
        //    Input, double X[N], Y[N], the data values.
        //
        //    Output, double VANDERMONDE_APPROX_1D_COEF[M+1], the coefficients of 
        //    the approximating polynomial.  C(0) is the constant term, and C(M) 
        //    multiplies X^M.
        //
    {
        double[] a = VandermondeMatrix.vandermonde_approx_1d_matrix(n, m, x);

        double[] c = QRSolve.qr_solve(n, m + 1, a, y);

        return c;
    }

    public static double[] vandermonde_approx_2d_coef(int n, int m, double[] x, double[] y,
            double[] z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VANDERMONDE_APPROX_2D_COEF computes a 2D polynomial approximant.
        //
        //  Discussion:
        //
        //    We assume the approximating function has the form of a polynomial
        //    in X and Y of total degree M.
        //
        //      p(x,y) = c00 
        //             + c10 * x                + c01 *  y
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
        //    This can be cast as an NxT(M+1) linear system for the polynomial
        //    coefficients:
        //
        //      [ 1 x1 y1  x1^2 ... y1^m ] [ c00 ] = [  z1 ]
        //      [ 1 x2 y2  x2^2 ... y2^m ] [ c10 ] = [  z2 ]
        //      [ 1 x3 y3  x3^2 ... y3^m ] [ c01 ] = [  z3 ]
        //      [ ...................... ] [ ... ] = [ ... ]
        //      [ 1 xn yn  xn^2 ... yn^m ] [ c0m ] = [  zn ]
        //
        //    In the typical case, N is greater than T(M+1) (we have more data and 
        //    equations than degrees of freedom) and so a least squares solution is 
        //    appropriate, in which case the computed polynomial will be a least squares
        //    approximant to the data.
        //
        //    The polynomial defined by the T(M+1) coefficients C could be evaluated 
        //    at the Nx2-vector x by the command
        //
        //      pval = r8poly_value_2d ( m, c, n, x )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of data points.
        //
        //    Input, int M, the maximum degree of the polynomial.
        //
        //    Input, double X[N], Y[N] the data locations.
        //
        //    Input, double Z[N], the data values.
        //
        //    Output, double VANDERMONDE_APPROX_2D_COEF[T(M+1)], the 
        //    coefficients of the approximating polynomial.  
        //
    {
        int tm = typeMethods.triangle_num(m + 1);

        double[] a = VandermondeMatrix.vandermonde_approx_2d_matrix(n, m, tm, x, y);

        double[] c = QRSolve.qr_solve(n, tm, a, z);

        return c;
    }
}