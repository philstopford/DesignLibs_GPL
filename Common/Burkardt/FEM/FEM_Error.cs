using System;

namespace Burkardt.FEM;

public static class FEM_Error
{
    public static double max_error_linear(int n, double[] x, double[] u,
        Func<double, double> exact)
//*****************************************************************************/
//
//  Purpose:
//
//    MAX_ERROR_LINEAR estimates the max error norm of a finite element solution.
//
//  Discussion:
//
//    We assume the finite element method has been used, over an interval [A,B]
//    involving N nodes, with piecewise linear elements used for the basis.
//
//    The coefficients U(1:N) have been computed, and a formula for the
//    exact solution is known.
//
//    This function estimates the max norm of the error:
//
//      MAX_NORM = Integral ( A <= X <= B ) max ( abs ( U(X) - EXACT(X) ) ) dX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Input, double X[N], the mesh points.
//
//    Input, double U[N], the finite element coefficients.
//
//    Input, function EQ = EXACT ( X ), returns the value of the exact
//    solution at the point X.
//
//    Output, double MAX_ERROR_LINEAR, the estimated max norm of the error.
//
    {
        int e;
        double eq;
        const int quad_num = 8;
        double uq;
        double xq;

        double value = 0.0;
//
//  Integrate over each interval.
//
        int e_num = n - 1;

        for (e = 0; e < e_num; e++)
        {
            double xl = x[e];
            double ul = u[e];

            int r = e + 1;
            double xr = x[r];
            double ur = u[r];

            int q;
            for (q = 0; q < quad_num; q++)
            {
                xq = ((quad_num - q) * xl
                      + q * xr)
                     / quad_num;
//
//  Use the fact that U is a linear combination of piecewise linears.
//
                uq = ((xr - xq) * ul
                      + (xq - xl) * ur)
                     / (xr - xl);

                eq = exact(xq);

                value = Math.Max(value, Math.Abs(uq - eq));
            }
        }

//
//  For completeness, check last node.
//
        xq = x[n - 1];
        uq = u[n - 1];
        eq = exact(xq);
        value = Math.Max(value, Math.Abs(uq - eq));
//
//  Integral approximation requires multiplication by interval length.
//
        value *= x[n - 1] - x[0];

        return value;
    }
        
    public static double max_error_quadratic(int n, double[] x, double[] u,
        Func<double, double> exact )
//*****************************************************************************/
//
//  Purpose:
//
//    MAX_ERROR_QUADRATIC estimates the max error norm of a finite element solution.
//
//  Discussion:
//
//    We assume the finite element method has been used, over an interval [A,B]
//    involving N nodes, with piecewise quadratic elements used for the basis.
//
//    The coefficients U(1:N) have been computed, and a formula for the
//    exact solution is known.
//
//    This function estimates the max norm of the error:
//
//      MAX_NORM = Integral ( A <= X <= B ) max ( abs ( U(X) - EXACT(X) ) ) dX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 July 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Input, double X[N], the mesh points.
//
//    Input, double U[N], the finite element coefficients.
//
//    Input, function EQ = EXACT ( X ), returns the value of the exact
//    solution at the point X.
//
//    Output, double MAX_ERROR_QUADRATIC, the estimated max norm of the error.
//
    {
        int e;
        double eq;
        const int quad_num = 8;
        double uq;
        double xq;

        double value = 0.0;
//
//  Integrate over each interval.
//
        int e_num = (n - 1) / 2;

        for (e = 0; e < e_num; e++)
        {
            int l = 2 * e;
            double xl = x[l];

            int m = 2 * e + 1;
            double xm = x[m];

            int r = 2 * e + 2;
            double xr = x[r];

            int q;
            for (q = 0; q < quad_num; q++)
            {
                xq = ((quad_num - q) * xl
                      + q * xr)
                     / quad_num;
                double vl = (xq - xm) / (xl - xm)
                            * ((xq - xr) / (xl - xr));

                double vm = (xq - xl) / (xm - xl)
                            * ((xq - xr) / (xm - xr));

                double vr = (xq - xl) / (xr - xl)
                            * ((xq - xm) / (xr - xm));

                uq = u[l] * vl + u[m] * vm + u[r] * vr;

                eq = exact(xq);

                value = Math.Max(value, Math.Abs(uq - eq));
            }
        }

//
//  For completeness, check last node.
//
        xq = x[n - 1];
        uq = u[n - 1];
        eq = exact(xq);
        value = Math.Max(value, Math.Abs(uq - eq));
//
//  Integral approximation requires multiplication by interval length.
//
        value *= x[n - 1] - x[0];

        return value;
    }
        
    public static double l1_error(int n, double[] x, double[] u, Func<double, double> exact)
//****************************************************************************80
//
//  Purpose:
//
//    L1_ERROR estimates the L1 error norm of a finite element solution.
//
//  Discussion:
//
//    We assume the finite element method has been used, over an interval [A,B]
//    involving N nodes.
//
//    The coefficients U(1:N) have been computed, and a formula for the
//    exact solution is known.
//
//    This function estimates the little l1 norm of the error:
//      L1_NORM = sum ( 1 <= I <= N ) abs ( U(i) - EXACT(X(i)) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Input, double X[N], the mesh points.
//
//    Input, double U[N], the finite element coefficients.
//
//    Input, function EQ = EXACT ( X ), returns the value of the exact
//    solution at the point X.
//
//    Output, double L1_ERROR, the estimated L2 norm of the error.
//
    {
        int i;

        double e1 = 0.0;

        for (i = 0; i < n; i++)
        {
            e1 += Math.Abs(u[i] - exact(x[i]));
        }

        e1 /= n;

        return e1;
    }

    public static double l2_error_linear(int n, double[] x, double[] u,
        Func<double, double> exact)
//****************************************************************************80
//
//  Purpose:
//
//    L2_ERROR_LINEAR estimates the L2 error norm of a finite element solution.
//
//  Discussion:
//
//    We assume the finite element method has been used, over an interval [A,B]
//    involving N nodes, with piecewise linear elements used for the basis.
//
//    The coefficients U(1:N) have been computed, and a formula for the
//    exact solution is known.
//
//    This function estimates the L2 norm of the error:
//
//      L2_NORM = Integral ( A <= X <= B ) ( U(X) - EXACT(X) )^2 dX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Input, double X[N], the mesh points.
//
//    Input, double U[N], the finite element coefficients.
//
//    Input, function EQ = EXACT ( X ), returns the value of the exact
//    solution at the point X.
//
//    Output, double L2_ERROR_LINEAR, the estimated L2 norm of the error.
//
    {
        const int QUAD_NUM = 2;

        double[] abscissa =
            {
                -0.577350269189625764509148780502,
                +0.577350269189625764509148780502
            }
            ;
        int i;
        double[] weight =
            {
                1.0, 1.0
            }
            ;

        double e2 = 0.0;
//
//  Integrate over each interval.
//
        for (i = 0; i < n - 1; i++)
        {
            double xl = x[i];
            double xr = x[i + 1];
            double ul = u[i];
            double ur = u[i + 1];

            int q;
            for (q = 0; q < QUAD_NUM; q++)
            {
                double xq = ((1.0 - abscissa[q]) * xl
                             + (1.0 + abscissa[q]) * xr)
                            / 2.0;

                double wq = weight[q] * (xr - xl) / 2.0;
//
//  Use the fact that U is a linear combination of piecewise linears.
//
                double uq = ((xr - xq) * ul
                             + (xq - xl) * ur)
                            / (xr - xl);

                double eq = exact(xq);

                e2 += wq * Math.Pow(uq - eq, 2);
            }
        }

        e2 = Math.Sqrt(e2);

        return e2;
    }

    public static double l2_error_quadratic(int n, double[] x, double[] u,
        Func<double, double> exact )
//****************************************************************************80
//
//  Purpose:
//
//    L2_ERROR_QUADRATIC estimates the L2 error norm of a finite element solution.
//
//  Discussion:
//
//    We assume the finite element method has been used, over an interval [A,B]
//    involving N nodes, with piecewise quadratic elements used for the basis.
//
//    The coefficients U(1:N) have been computed, and a formula for the
//    exact solution is known.
//
//    This function estimates the L2 norm of the error:
//
//      L2_NORM = Integral ( A <= X <= B ) ( U(X) - EXACT(X) )^2 dX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Input, double X[N], the mesh points.
//
//    Input, double U[N], the finite element coefficients.
//
//    Input, function EQ = EXACT ( X ), returns the value of the exact
//    solution at the point X.
//
//    Output, double L2_ERROR_QUADRATIC, the estimated L2 norm of the error.
//
    {
        const int QUAD_NUM = 3;

        double[] abscissa =  {
                -0.774596669241483377035853079956,
                0.000000000000000000000000000000,
                0.774596669241483377035853079956
            }
            ;
        int e;
        double[] weight =  {
                0.555555555555555555555555555556,
                0.888888888888888888888888888889,
                0.555555555555555555555555555556
            }
            ;

        double e2 = 0.0;
//
//  Integrate over element E.
//
        int e_num = (n - 1) / 2;

        for (e = 0; e < e_num; e++)
        {
            int l = 2 * e;
            int m = 2 * e + 1;
            int r = 2 * e + 2;

            double xl = x[l];
            double xm = x[m];
            double xr = x[r];

            int q;
            for (q = 0; q < QUAD_NUM; q++)
            {

                double xq = ((1.0 - abscissa[q]) * xl
                             + (1.0 + abscissa[q]) * xr)
                            / 2.0;

                double wq = weight[q] * (xr - xl) / 2.0;

                double vl = (xq - xm) / (xl - xm)
                            * ((xq - xr) / (xl - xr));

                double vm = (xq - xl) / (xm - xl)
                            * ((xq - xr) / (xm - xr));

                double vr = (xq - xl) / (xr - xl)
                            * ((xq - xm) / (xr - xm));

                double uq = u[l] * vl + u[m] * vm + u[r] * vr;
                double eq = exact(xq);

                e2 += wq * Math.Pow(uq - eq, 2);
            }
        }

        e2 = Math.Sqrt(e2);

        return e2;
    }
        
        
    public static double h1s_error_linear(int n, double[] x, double[] u,
        Func<double, double> exact_ux)
//****************************************************************************80
//
//  Purpose:
//
//    H1S_ERROR_LINEAR estimates the seminorm error of a finite element solution.
//
//  Discussion:
//
//    We assume the finite element method has been used, over an interval [A,B]
//    involving N nodes, with piecewise linear elements used for the basis.
//
//    The coefficients U(1:N) have been computed, and a formula for the
//    exact derivative is known.
//
//    This function estimates the seminorm of the error:
//
//      SEMINORM = Integral ( A <= X <= B ) ( dU(X)/dx - EXACT_UX(X) )^2 dX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Input, double X[N], the mesh points.
//
//    Input, double U[N], the finite element coefficients.
//
//    Input, function EQ = EXACT_UX ( X ), returns the value of the exact
//    derivative at the point X.
//
//    Output, double H1S_ERROR_LINEAR, the estimated seminorm of 
//    the error.
//
    {
        const int QUAD_NUM = 2;

        double[] abscissa =
        {
            -0.577350269189625764509148780502,
            +0.577350269189625764509148780502
        };
        int i;
        double[] weight = {1.0, 1.0};

        double h1s = 0.0;
//
//  Integrate over each interval.
//
        for (i = 0; i < n - 1; i++)
        {
            double xl = x[i];
            double xr = x[i + 1];
            double ul = u[i];
            double ur = u[i + 1];

            int q;
            for (q = 0; q < QUAD_NUM; q++)
            {
                double xq = ((1.0 - abscissa[q]) * xl
                             + (1.0 + abscissa[q]) * xr)
                            / 2.0;

                double wq = weight[q] * (xr - xl) / 2.0;
//
//  The piecewise linear derivative is a constant in the interval.
//
                double uxq = (ur - ul) / (xr - xl);

                double exq = exact_ux(xq);

                h1s += wq * Math.Pow(uxq - exq, 2);
            }
        }

        h1s = Math.Sqrt(h1s);

        return h1s;
    }

    public static double h1s_error_quadratic(int n, double[] x, double[] u,
        Func<double, double> exact_ux)
//****************************************************************************80
//
//  Purpose:
//
//    H1S_ERROR_QUADRATIC estimates the seminorm error of a finite element solution.
//
//  Discussion:
//
//    We assume the finite element method has been used, over an interval [A,B]
//    involving N nodes, with piecewise quadratic elements used for the basis.
//
//    The coefficients U(1:N) have been computed, and a formula for the
//    exact derivative is known.
//
//    This function estimates the seminorm of the error:
//
//      SEMINORM = Integral ( A <= X <= B ) ( dU(X)/dx - EXACT_UX(X) )^2 dX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Input, double X[N], the mesh points.
//
//    Input, double U[N], the finite element coefficients.
//
//    Input, function EQ = EXACT_UX ( X ), returns the value of the exact
//    derivative at the point X.
//
//    Output, double H1S_ERROR_QUADRATIC, the estimated seminorm of 
//    the error.
//
    {
        const int QUAD_NUM = 3;

        double[] abscissa =
        {
            -0.774596669241483377035853079956,
            0.000000000000000000000000000000,
            0.774596669241483377035853079956
        };
        int e;
        double[] weight =
        {
            0.555555555555555555555555555556,
            0.888888888888888888888888888889,
            0.555555555555555555555555555556
        };

        double h1s = 0.0;
//
//  Integrate over element E.
//
        int e_num = (n - 1) / 2;

        for (e = 0; e < e_num; e++)
        {
            int l = 2 * e;
            int m = 2 * e + 1;
            int r = 2 * e + 2;

            double xl = x[l];
            double xm = x[m];
            double xr = x[r];

            int q;
            for (q = 0; q < QUAD_NUM; q++)
            {

                double xq = ((1.0 - abscissa[q]) * xl
                             + (1.0 + abscissa[q]) * xr)
                            / 2.0;

                double wq = weight[q] * (xr - xl) / 2.0;

                double vlp = 1.0 / (xl - xm)
                             * ((xq - xr) / (xl - xr))
                             + (xq - xm) / (xl - xm)
                             * (1.0 / (xl - xr));

                double vmp = 1.0 / (xm - xl)
                             * ((xq - xr) / (xm - xr))
                             + (xq - xl) / (xm - xl)
                             * (1.0 / (xm - xr));

                double vrp = 1.0 / (xr - xl)
                             * ((xq - xm) / (xr - xm))
                             + (xq - xl) / (xr - xl)
                             * (1.0 / (xr - xm));

                double uxq = u[l] * vlp + u[m] * vmp + u[r] * vrp;

                double exq = exact_ux(xq);

                h1s += wq * Math.Pow(uxq - exq, 2);
            }
        }

        h1s = Math.Sqrt(h1s);

        return h1s;
    }        
        
        
}