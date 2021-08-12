using Burkardt.SolveNS;
using Burkardt.Types;

namespace Burkardt.Lagrange
{
    public static class Lagrange1D
    {
        public static double[] lagrange_approx_1d(int m, int nd, double[] xd, double[] yd,
        int ni, double[] xi )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_APPROX_1D evaluates the Lagrange approximant of degree M.
        //
        //  Discussion:
        //
        //    The Lagrange approximant L(M,ND,XD,YD)(X) is a polynomial of
        //    degree M which approximates the data (XD(I),YD(I)) for I = 1 to ND.
        //
        //    We can represent any polynomial of degree M+1 as the sum of the Lagrange 
        //    basis functions at the M+1 Chebyshev points.
        //
        //      L(M)(X) = sum ( 1 <= I <= M+1 ) C(I) LB(M,XC)(X)
        //
        //    Given our data, we can seek the M+1 unknown coefficients C which minimize
        //    the norm of || L(M)(XD(1:ND)) - YD(1:ND) ||.
        //
        //    Given the coefficients, we can then evaluate the polynomial at the
        //    points XI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the polynomial degree.
        //
        //    Input, int ND, the number of data points.
        //    ND must be at least 1.
        //
        //    Input, double XD[ND], the data points.
        //
        //    Input, double YD[ND], the data values.
        //
        //    Input, int NI, the number of interpolation points.
        //
        //    Input, double XI[NI], the interpolation points.
        //
        //    Output, double LAGRANGE_APPROX_1D[NI], the interpolated values.
        //
        {
            double a;
            double b;
            double[] ld;
            double[] li;
            int nc;
            double[] xc;
            double[] yc;
            double[] yi;

            nc = m + 1;
            //
            //  Evaluate the Chebyshev points.
            //
            a = -1.0;
            b = +1.0;
            xc = typeMethods.r8vec_cheby_extreme_new(nc, a, b);
            //
            //  Evaluate the Lagrange basis functions for the Chebyshev points 
            //  at the data points.
            //
            ld = lagrange_basis_1d(nc, xc, nd, xd);
            //
            //  The value of the Lagrange approximant at each data point should
            //  approximate the data value: LD * YC = YD, where YC are the unknown
            //  coefficients.
            //
            yc = QRSolve.qr_solve(nd, nc, ld, yd);
            //
            //  Now we want to evaluate the Lagrange approximant at the "interpolant
            //  points": LI * YC = YI
            //
            li = lagrange_basis_1d(nc, xc, ni, xi);

            yi = typeMethods.r8mat_mv_new(ni, nc, li, yc);

            return yi;
        }
        
        public static double[] lagrange_base_1d ( int nd, double[] xd, int ni, double[] xi, int xdIndex = 0, int xiIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_BASE_1D evaluates the Lagrange basis polynomials.
        //
        //  Discussion:
        //
        //    Given ND distinct abscissas, XD(1:ND),
        //    the I-th Lagrange basis polynomial LB(I)(T) is defined as the polynomial of
        //    degree ND - 1 which is 1 at XD(I) and 0 at the ND - 1
        //    other abscissas.
        //
        //    A formal representation is:
        //
        //      LB(I)(T) = Product ( 1 <= J <= ND, I /= J )
        //       ( T - T(J) ) / ( T(I) - T(J) )
        //
        //    This routine accepts a set of NI values at which all the Lagrange
        //    basis polynomials should be evaluated.
        //
        //    Given data values YD at each of the abscissas, the value of the
        //    Lagrange interpolating polynomial at each of the interpolation points
        //    is then simple to compute by matrix multiplication:
        //
        //      YI(1:NI) = LB(1:NI,1:ND) * YD(1:ND)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ND, the number of data points.
        //    ND must be at least 1.
        //
        //    Input, double XD[ND], the data points.
        //
        //    Input, int NI, the number of interpolation points.
        //
        //    Input, double XI[NI], the interpolation points.
        //
        //    Output, double LAGRANGE_BASIS[NI*ND], the values
        //    of the Lagrange basis polynomials at the interpolation points.
        //
        {
            int i;
            int j;
            int k;
            double[] lb;
            //
            //  Evaluate the polynomial.
            //
            lb = new double[ni * nd];

            for (j = 0; j < nd; j++)
            {
                for (i = 0; i < ni; i++)
                {
                    lb[i + j * ni] = 1.0;
                }
            }

            for (i = 0; i < nd; i++)
            {
                for (j = 0; j < nd; j++)
                {
                    if (j != i)
                    {
                        for (k = 0; k < ni; k++)
                        {
                            lb[k + i * ni] = lb[k + i * ni] * (xi[xiIndex + k] - xd[xdIndex + j]) / (xd[xdIndex + i] - xd[xdIndex + j]);
                        }
                    }
                }
            }

            return lb;
        }


        public static double[] lagrange_basis_1d(int nd, double[] xd, int ni, double[] xi )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_BASIS_1D evaluates the Lagrange basis polynomials.
        //
        //  Discussion:
        //
        //    Given ND distinct abscissas, XD(1:ND),
        //    the I-th Lagrange basis polynomial LB(I)(T) is defined as the polynomial of
        //    degree ND - 1 which is 1 at XD(I) and 0 at the ND - 1
        //    other abscissas.
        //
        //    A formal representation is:
        //
        //      LB(I)(T) = Product ( 1 <= J <= ND, I /= J )
        //       ( T - T(J) ) / ( T(I) - T(J) )
        //
        //    This routine accepts a set of NI values at which all the Lagrange
        //    basis polynomials should be evaluated.
        //
        //    Given data values YD at each of the abscissas, the value of the
        //    Lagrange interpolating polynomial at each of the interpolation points
        //    is then simple to compute by matrix multiplication:
        //
        //      YI(1:NI) = LB(1:NI,1:ND) * YD(1:ND)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ND, the number of data points.
        //    ND must be at least 1.
        //
        //    Input, double XD[ND], the data points.
        //
        //    Input, int NI, the number of interpolation points.
        //
        //    Input, double XI[NI], the interpolation points.
        //
        //    Output, double LAGRANGE_BASIS[NI*ND], the values
        //    of the Lagrange basis polynomials at the interpolation points.
        //
        {
            int i;
            int j;
            int k;
            double[] lb;
            //
            //  Evaluate the polynomial.
            //
            lb = new double[ni * nd];

            for (j = 0; j < nd; j++)
            {
                for (i = 0; i < ni; i++)
                {
                    lb[i + j * ni] = 1.0;
                }
            }

            for (i = 0; i < nd; i++)
            {
                for (j = 0; j < nd; j++)
                {
                    if (j != i)
                    {
                        for (k = 0; k < ni; k++)
                        {
                            lb[k + i * ni] = lb[k + i * ni] * (xi[k] - xd[j]) / (xd[i] - xd[j]);
                        }
                    }
                }
            }

            return lb;
        }
        
        public static double lagrange_basis_function_1d ( int mx, double[] xd, int i, double xi ) 

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_BASIS_FUNCTION_1D evaluates a 1D Lagrange basis function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int MX, the degree of the basis function.
        //
        //    Input, double XD[MX+1], the interpolation nodes.
        //
        //    Input, int I, the index of the basis function.
        //    0 <= I <= MX.
        //
        //    Input, double XI, the evaluation point.
        //
        //    Output, double LAGRANGE_BASIS_FUNCTION_1D, the value of the I-th Lagrange 1D 
        //    basis function for the nodes XD, evaluated at XI.
        //
        {
            int j;
            double yi;

            yi = 1.0;

            if ( xi != xd[i] )
            {
                for ( j = 0; j < mx + 1; j++ )
                {
                    if ( j != i )
                    {
                        yi = yi * ( xi - xd[j] ) / ( xd[i] - xd[j] );
                    }
                }
            }

            return yi;
        }
        
        public static double[] lagrange_value_1d ( int nd, double[] xd, double[] yd, int ni, 
        double[] xi )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_VALUE_1D evaluates the Lagrange interpolant.
        //
        //  Discussion:
        //
        //    The Lagrange interpolant L(ND,XD,YD)(X) is the unique polynomial of
        //    degree ND-1 which interpolates the points (XD(I),YD(I)) for I = 1
        //    to ND.
        //
        //    The Lagrange interpolant can be constructed from the Lagrange basis
        //    polynomials.  Given ND distinct abscissas, XD(1:ND), the I-th Lagrange 
        //    basis polynomial LB(ND,XD,I)(X) is defined as the polynomial of degree 
        //    ND - 1 which is 1 at  XD(I) and 0 at the ND - 1 other abscissas.
        //
        //    Given data values YD at each of the abscissas, the value of the
        //    Lagrange interpolant may be written as
        //
        //      L(ND,XD,YD)(X) = sum ( 1 <= I <= ND ) LB(ND,XD,I)(X) * YD(I)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ND, the number of data points.
        //    ND must be at least 1.
        //
        //    Input, double XD[ND], the data points.
        //
        //    Input, double YD[ND], the data values.
        //
        //    Input, int NI, the number of interpolation points.
        //
        //    Input, double XI[NI], the interpolation points.
        //
        //    Output, double LAGRANGE_VALUE_1D[NI], the interpolated values.
        //
        {
            double[] lb;
            double[] yi;

            lb = lagrange_basis_1d ( nd, xd, ni, xi );

            yi = typeMethods.r8mat_mv_new ( ni, nd, lb, yd );
            
            return yi;
        }
    }
}