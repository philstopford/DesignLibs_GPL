using System;
using Burkardt.Types;

namespace Burkardt.Disk
{
    public static class QuadratureRule
    {
        public static void disk_rule(int nr, int nt, double xc, double yc, double rc, ref double[] w,
        ref double[] x, ref double[] y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISK_RULE computes a quadrature rule for the general disk.
        //
        //  Discussion:
        //
        //    The general disk is the region:
        //
        //      ( x - xc ) ^ 2 + ( y - yc ) ^ 2 <= rc ^ 2.
        //
        //    The integral I(f) is then approximated by
        //
        //      S(f) = sum ( 1 <= i <= NT * NR ) W(i) * F ( X(i), Y(i) ).
        //
        //      Area = pi * RC ^ 2
        //
        //      Q(f) = Area * S(f)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 March 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NR, the number of points in the radial rule.
        //
        //    Input, int NT, the number of angles to use.
        //
        //    Input, double XC, YC, the center of the disk.
        //
        //    Input, double RC, the radius of the disk.
        //
        //    Output, double W[NR*NT], the weights for the rule.
        //
        //    Output, double X[NR*NT], Y[NR*NT], the points for the rule.
        //
        {
            int i;
            int j;
            double[] r01;
            double[] t01;
            double[] w01;

            w01 = new double[nr];
            r01 = new double[nr];
            t01 = new double[nt];

            disk01_rule(nr, nt, ref w01, ref r01, ref t01);
            //
            //  Recompute the rule for the general circle in terms of X, Y.
            //
            for (j = 0; j < nt; j++)
            {
                for (i = 0; i < nr; i++)
                {
                    w[i + j * nr] = w01[i];
                    x[i + j * nr] = xc + rc * r01[i] * Math.Cos(t01[j]);
                    y[i + j * nr] = yc + rc * r01[i] * Math.Sin(t01[j]);
                }
            }
        }

        public static double disk01_monomial_integral(int[] e )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISK01_MONOMIAL_INTEGRAL returns monomial integrals in the unit disk in 2D.
        //
        //  Discussion:
        //
        //    The integration region is 
        //
        //      X^2 + Y^2 <= 1.
        //
        //    The monomial is F(X,Y) = X^E(1) * Y^E(2).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Academic Press, 1984, page 263.
        //
        //  Parameters:
        //
        //    Input, int E[2], the exponents of X and Y in the 
        //    monomial.  Each exponent must be nonnegative.
        //
        //    Output, double DISK01_MONOMIAL_INTEGRAL, the integral.
        //
        {
            double arg;
            int i;
            double integral;
            const double r = 1.0;
            double s;

            if (e[0] < 0 || e[1] < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("DISK01_MONOMIAL_INTEGRAL - Fatal error!");
                Console.WriteLine("  All exponents must be nonnegative.");
                Console.WriteLine("  E[0] = " + e[0] + "");
                Console.WriteLine("  E[1] = " + e[1] + "");
                return 1;
            }

            if ((e[0] % 2) == 1 || (e[1] % 2) == 1)
            {
                integral = 0.0;
            }
            else
            {
                integral = 2.0;

                for (i = 0; i < 2; i++)
                {
                    arg = 0.5 * (double) (e[i] + 1);
                    integral = integral * typeMethods.r8_gamma(arg);
                }

                arg = 0.5 * (double) (e[0] + e[1] + 2);
                integral = integral / typeMethods.r8_gamma(arg);
            }

            //
            //  Adjust the surface integral to get the volume integral.
            //
            s = e[0] + e[1] + 2;
            integral = integral * Math.Pow(r, s) / (double) (s);

            return integral;
        }

        public static void disk01_rule(int nr, int nt, ref double[] w, ref double[] r, ref double[] t )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISK01_RULE computes a quadrature rule for the unit disk.
        //
        //  Discussion:
        //
        //    The unit disk is the region:
        //
        //      x * x + y * y <= 1.
        //
        //    The integral I(f) is then approximated by
        //
        //      Q(f) = pi * sum ( 1 <= j <= NT ) sum ( 1 <= i <= NR ) 
        //        W(i) * F ( R(i) * cos(T(j)), R(i) * sin(T(j)) ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 March 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NR, the number of points in the radial rule.
        //
        //    Input, int NT, the number of angles to use.
        //
        //    Output, double W[NR], the weights for the disk rule.
        //
        //    Output, double R[NR], T[NT], the (R,Theta) points for the rule.
        //
        {
            int ir;
            int it;
            const double r8_pi = 3.141592653589793;
            double[] wr;
            double[] xr;
            //
            //  Request a Legendre rule for [-1,+1].
            //
            xr = new double[nr];
            wr = new double[nr];

            Legendre.legendre_ek_compute(nr, ref xr, ref wr);
            //
            //  Shift the rule to [0,1].
            //
            for (ir = 0; ir < nr; ir++)
            {
                xr[ir] = (xr[ir] + 1.0) / 2.0;
                wr[ir] = wr[ir] / 2.0;
            }

            //
            //  Compute the disk rule.
            //
            for (it = 0; it < nt; it++)
            {
                t[it] = 2.0 * r8_pi * (double) (it) / (double) (nt);
            }

            for (ir = 0; ir < nr; ir++)
            {
                w[ir] = wr[ir] / (double) (nt);
                r[ir] = Math.Sqrt(xr[ir]);
            }
        }
    }
}