using System;
using Burkardt.PolynomialNS;
using Burkardt.Types;

namespace Burkardt.Interpolation
{
    public static class Hermite
    {
        public static double hermite_basis_0(int n, double[] x, int i, double xv)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_BASIS_0 evaluates a zero-order Hermite interpolation basis function.
            //
            //  Discussion:
            //
            //    Given ND points XD, with values YD and derivative values YPD, the
            //    Hermite interpolant can be written as:
            //
            //      H(X) = sum ( 1 <= I <= ND ) YD(I)  * H0(I;X)
            //           + sum ( 1 <= I <= ND ) YPD(I) * H1(I;X)
            //
            //    where H0(I;X) is the I-th zero order Hermite interpolation basis function,
            //    and H1(I;X) is the I-th first order Hermite interpolation basis function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 May 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of abscissas.
            //
            //    Input, double X[N], the abscissas.
            //
            //    Input, int I, the index of the first-order basis function.
            //    Indices are 0-based
            //
            //    Input, double XV, the evaluation point.
            //
            //    Output, double HERMITE_BASIS_0, the value of the function.
            //
        {
            double[] factor;
            int j;
            double li;
            double lp;
            double lpp;
            double value;

            if (i < 0 || n - 1 < i)
            {
                Console.WriteLine("");
                Console.WriteLine("HERMITE_BASIS_0 - Fatal error!");
                Console.WriteLine("  I < 0 or N - 1 < I.");
                return (1);
            }

            factor = new double[n];
            //
            //  L(X) = product ( X - X(1:N) )
            //
            //  L'(X(I)).
            //
            for (j = 0; j < n; j++)
            {
                factor[j] = x[i] - x[j];
            }

            factor[i] = 1.0;

            lp = typeMethods.r8vec_product(n, factor);
            //
            //  LI(X) = L(X) / ( X - X(I) ) / L'(X(I))
            //
            for (j = 0; j < n; j++)
            {
                factor[j] = xv - x[j];
            }

            factor[i] = 1.0;

            li = typeMethods.r8vec_product(n, factor) / lp;
            //
            //  L''(X(I)).
            //
            lpp = 0.0;
            for (j = 0; j < n; j++)
            {
                factor[j] = x[i] - x[j];
            }

            factor[i] = 1.0;

            for (j = 0; j < n; j++)
            {
                if (j != i)
                {
                    factor[j] = 1.0;
                    lpp = lpp + 2.0 * typeMethods.r8vec_product(n, factor);
                    factor[j] = x[i] - x[j];
                }
            }

            value = (1.0 - (xv - x[i]) * lpp / lp) * li * li;

            return value;
        }

        public static double hermite_basis_1(int n, double[] x, int i, double xv)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_BASIS_1 evaluates a first-order Hermite interpolation basis function.
            //
            //  Discussion:
            //
            //    Given ND points XD, with values YD and derivative values YPD, the
            //    Hermite interpolant can be written as:
            //
            //      H(X) = sum ( 1 <= I <= ND ) YD(I)  * H0(I;X)
            //           + sum ( 1 <= I <= ND ) YPD(I) * H1(I;X)
            //
            //    where H0(I;X) is the I-th zero order Hermite interpolation basis function,
            //    and H1(I;X) is the I-th first order Hermite interpolation basis function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 May 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of abscissas.
            //
            //    Input, double X[N], the abscissas.
            //
            //    Input, int I, the index of the first-order basis function.
            //    Indices are 0-based
            //
            //    Input, double XV, the evaluation point.
            //
            //    Output, double VALUE, the value of the function.
            //
        {
            double bot;
            double[] factor;
            int j;
            double top;
            double value;

            if (i < 0 || n - 1 < i)
            {
                Console.WriteLine("");
                Console.WriteLine("HERMITE_BASIS_1 - Fatal error!");
                Console.WriteLine("  I < 0 or N - 1 < I.");
                return (1);
            }

            factor = new double[n];

            for (j = 0; j < n; j++)
            {
                factor[j] = xv - x[j];
            }

            factor[i] = 1.0;
            top = typeMethods.r8vec_product(n, factor);

            for (j = 0; j < n; j++)
            {
                factor[j] = x[i] - x[j];
            }

            factor[i] = 1.0;
            bot = typeMethods.r8vec_product(n, factor);

            value = (xv - x[i]) * (top / bot) * (top / bot);

            return value;
        }

        public static void hermite_demo(int n, double[] x, double[] y, double[] yp)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_DEMO computes and prints Hermite interpolant information for data.
            //
            //  Discussion:
            //
            //    Given a set of Hermite data, this routine calls HERMITE_INTERPOLANT to 
            //    determine and print the divided difference table, and then DIF_TO_R8POLY to 
            //    determine and print the coefficients of the polynomial in standard form.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 October 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of data points.
            //
            //    Input, double X[N], the abscissas.
            //
            //    Input, double Y[N], YP[N], the function and derivative
            //    values at the abscissas.
            //
        {
            double[] cd;
            int i;
            int nd;
            int ndp;
            int nv;
            double[] xd;
            double[] xdp;
            double[] xv;
            double[] yd;
            double[] ydp;
            double[] yv;
            double[] yvp;

            Console.WriteLine("");
            Console.WriteLine("HERMITE_DEMO");
            Console.WriteLine("  Compute coefficients CD of the Hermite polynomial");
            Console.WriteLine("  interpolant to given data (x,y,yp).");

            Console.WriteLine("");
            Console.WriteLine("  Data:");
            Console.WriteLine("              X           Y           Y'");
            Console.WriteLine("");
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(4)
                                       + "  " + x[i].ToString().PadLeft(10)
                                       + "  " + y[i].ToString().PadLeft(10)
                                       + "  " + yp[i].ToString().PadLeft(10) + "");
            }

            nd = 2 * n;
            xd = new double[nd];
            yd = new double[nd];

            ndp = 2 * n - 1;
            xdp = new double[ndp];
            ydp = new double[ndp];

            hermite_interpolant(n, x, y, yp, ref xd, ref yd, ref xdp, ref ydp);

            Console.WriteLine("");
            Console.WriteLine("  Difference table:");
            Console.WriteLine("              XD          YD");
            Console.WriteLine("");

            for (i = 0; i < nd; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(4)
                                       + "  " + xd[i].ToString().PadLeft(10)
                                       + "  " + yd[i].ToString().PadLeft(10) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  Difference table:");
            Console.WriteLine("              XDP          YDP");
            Console.WriteLine("");
            for (i = 0; i < nd - 1; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(4)
                                       + "  " + xdp[i].ToString().PadLeft(10)
                                       + "  " + ydp[i].ToString().PadLeft(10) + "");
            }

            cd = new double[2 * n];

            Dif.dif_to_r8poly(nd, xd, yd, ref cd);

            typeMethods.r8poly_print(nd - 1, cd, "  Hermite interpolating polynomial:");
            //
            //  Verify interpolation claim!
            //
            nv = n;
            xv = new double[nv];
            yv = new double[nv];
            yvp = new double[nv];

            for (i = 0; i < nv; i++)
            {
                xv[i] = x[i];
            }

            hermite_interpolant_value(nd, xd, yd, xdp, ydp, nv, xv, ref yv, ref yvp);

            Console.WriteLine("");
            Console.WriteLine("  Data Versus Interpolant:");
            Console.WriteLine("              X           Y           H           YP          HP");
            Console.WriteLine("");
            for (i = 0; i < nv; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(4)
                                       + "  " + xv[i].ToString().PadLeft(10)
                                       + "  " + y[i].ToString().PadLeft(10)
                                       + "  " + yv[i].ToString().PadLeft(10)
                                       + "  " + yp[i].ToString().PadLeft(10)
                                       + "  " + yvp[i].ToString().PadLeft(10) + "");
            }
        }

        public static void hermite_interpolant(int n, double[] x, double[] y, double[] yp,
                ref double[] xd, ref double[] yd, ref double[] xdp, ref double[] ydp)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_INTERPOLANT sets up a divided difference table from Hermite data.
            //
            //  Discussion:
            //
            //    The polynomial represented by the divided difference table can be
            //    evaluated by calling DIF_VALS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 October 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carl deBoor,
            //    A Practical Guide to Splines,
            //    Springer, 2001,
            //    ISBN: 0387953663,
            //    LC: QA1.A647.v27.
            //
            //  Parameters:
            //
            //    Input, int N, of items of data 
            //    ( X(I), Y(I), YP(I) ).
            //
            //    Input, double X[N], the abscissas.
            //    These values must be distinct.
            //
            //    Input, double Y[N], YP[N], the function and derivative values.
            //
            //    Output, double XD[2*N], YD[2*N], the divided difference table
            //    for the interpolant value.
            //
            //    Output, double XDP[2*N-1], YDP[2*N-1], the divided difference 
            //    table for the interpolant derivative.
            //
        {
            int i;
            int j;
            int nd;
            int ndp = 0;
            //
            //  Copy the data.
            //
            nd = 2 * n;

            for (i = 0; i < n; i++)
            {
                xd[0 + i * 2] = x[i];
                xd[1 + i * 2] = x[i];
            }

            //
            //  Carry out the first step of differencing.
            //
            yd[0] = y[0];
            for (i = 1; i < n; i++)
            {
                yd[0 + 2 * i] = (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
            }

            for (i = 0; i < n; i++)
            {
                yd[1 + 2 * i] = yp[i];
            }

            //
            //  Carry out the remaining steps in the usual way.
            //
            for (i = 2; i < nd; i++)
            {
                for (j = nd - 1; i <= j; j--)
                {
                    yd[j] = (yd[j] - yd[j - 1]) / (xd[j] - xd[j - i]);
                }
            }

            //
            //  Compute the difference table for the derivative.
            //
            Dif.dif_deriv(nd, xd, yd, ref ndp, xdp, ref ydp);

        }

        public static double[] hermite_interpolant_rule(int n, double a, double b, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_INTERPOLANT_RULE: quadrature rule for a Hermite interpolant.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 October 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of abscissas.
            //
            //    Input, double A, B, the integration limits.
            //
            //    Input, double X[N], the abscissas.
            //
            //    Output, double HERMITE_INTERPOLANT_RULE[2*N], the quadrature 
            //    coefficients, given as pairs for function and derivative values 
            //    at each abscissa.
            //
        {
            double a_value;
            double b_value;
            double[] c;
            int i;
            int k;
            int nd;
            int ndp;
            double[] w;
            double[] xd;
            double[] xdp;
            double[] y;
            double[] yd;
            double[] ydp;
            double[] yp;

            y = new double[n];
            yp = new double[n];

            nd = 2 * n;
            c = new double[nd];
            w = new double[nd];
            xd = new double[nd];
            yd = new double[nd];

            ndp = 2 * n - 1;
            xdp = new double[ndp];
            ydp = new double[ndp];

            for (i = 0; i < n; i++)
            {
                y[i] = 0.0;
                yp[i] = 0.0;
            }

            k = 0;

            for (i = 0; i < n; i++)
            {
                y[i] = 1.0;
                hermite_interpolant(n, x, y, yp, ref xd, ref yd, ref xdp, ref ydp);
                Dif.dif_to_r8poly(nd, xd, yd, ref c);
                a_value = typeMethods.r8poly_ant_val(n, c, a);
                b_value = typeMethods.r8poly_ant_val(n, c, b);
                w[k] = b_value - a_value;
                y[i] = 0.0;
                k = k + 1;

                yp[i] = 1.0;
                hermite_interpolant(n, x, y, yp, ref xd, ref yd, ref xdp, ref ydp);
                Dif.dif_to_r8poly(nd, xd, yd, ref c);
                a_value = typeMethods.r8poly_ant_val(n, c, a);
                b_value = typeMethods.r8poly_ant_val(n, c, b);
                w[k] = b_value - a_value;
                yp[i] = 0.0;
                k = k + 1;
            }

            return w;
        }

        public static void hermite_interpolant_value(int nd, double[] xd, double[] yd, double[] xdp,
                double[] ydp, int nv, double[] xv, ref double[] yv, ref double[] yvp)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_INTERPOLANT_VALUE evaluates the Hermite interpolant polynomial.
            //
            //  Discussion:
            //
            //    In fact, this function will evaluate an arbitrary polynomial that is
            //    represented by a difference table.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 October 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carl deBoor,
            //    A Practical Guide to Splines,
            //    Springer, 2001,
            //    ISBN: 0387953663,
            //    LC: QA1.A647.v27.
            //
            //  Parameters:
            //
            //    Input, int ND, the order of the difference table.
            //
            //    Input, double XD[ND], YD[ND], the difference table for the
            //    interpolant value.
            //
            //    Input, double XDP[ND-1], YDP[ND-1], the difference table for
            //    the interpolant derivative.
            //
            //    Input, int NV, the number of evaluation points.
            //
            //    Input, double XV[NV], the evaluation points.
            //
            //    Output, double YV[NV], YVP[NV], the value of the interpolant and
            //    its derivative at the evaluation points.
            //
        {
            int i;
            int j;
            int ndp;

            ndp = nd - 1;

            for (j = 0; j < nv; j++)
            {
                yv[j] = yd[nd - 1];
                for (i = nd - 2; 0 <= i; i--)
                {
                    yv[j] = yd[i] + (xv[j] - xd[i]) * yv[j];
                }

                yvp[j] = ydp[ndp - 1];
                for (i = ndp - 2; 0 <= i; i--)
                {
                    yvp[j] = ydp[i] + (xv[j] - xdp[i]) * yvp[j];
                }
            }
        }
    }
}