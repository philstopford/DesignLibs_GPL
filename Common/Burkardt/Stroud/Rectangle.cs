using System;
using Burkardt.Types;

namespace Burkardt.Stroud;

public static class Rectangle
{


    public static double rectangle_3d(int settings, Func<int, double, double, double, double> func,
            double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RECTANGLE_3D approximates an integral inside a rectangular block in 3D.
        //
        //  Integration region:
        //
        //      A(1) <= X <= B(1),
        //    and
        //      A(2) <= Y <= B(2),
        //    and
        //      A(3) <= Z <= B(3).
        //
        //  Discussion:
        //
        //    An 8 point third degree formula is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Arthur Stroud,
        //    Approximate Calculation of Multiple Integrals,
        //    Prentice Hall, 1971,
        //    ISBN: 0130438936,
        //    LC: QA311.S85.
        //
        //  Parameters:
        //
        //    Input, Func< double, double, double, double > func, the name of the 
        //    user supplied function.
        //
        //    Input, double A[3], B[3], the lower and upper limits
        //    for X, Y and Z.
        //
        //    Output, double RECTANGLE_3D, the approximate integral of the function.
        //
    {
        int i;

        double sqr3 = 1.0 / Math.Sqrt(3.0);
        const double w = 1.0 / 8.0;

        double quad = 0.0;

        for (i = 1; i <= 2; i++)
        {
            double x = sqr3 * (int)Math.Pow(-1, i);
            x = 0.5 * ((1.0 - x) * b[0] + (1.0 + x) * a[0]);

            int j;
            for (j = 1; j <= 2; j++)
            {
                double y = sqr3 * (int)Math.Pow(-1, j);
                y = 0.5 * ((1.0 - y) * b[1] + (1.0 + y) * a[1]);

                int k;
                for (k = 1; k <= 2; k++)
                {
                    double z = sqr3 * (int)Math.Pow(-1, k);
                    z = 0.5 * ((1.0 - z) * b[2] + (1.0 + z) * a[2]);

                    quad += w * func(settings, x, y, z);
                }
            }
        }

        double volume = (b[0] - a[0]) * (b[1] - a[1]) * (b[2] - a[2]);
        double result = volume * quad;

        return result;
    }

    public static double rectangle_sub_2d(int setting, Func<int, double, double, double> func, double[] xval,
            double[] yval, int[] nsub, int order, double[] xtab, double[] ytab,
            double[] weight)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RECTANGLE_SUB_2D carries out a composite quadrature over a rectangle in 2D.
        //
        //  Integration region:
        //
        //      XVAL(1) <= X <= XVAL(2),
        //    and
        //      YVAL(1) <= Y <= YVAL(2).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, Func < double, double, double > func, the name of the function 
        //    to be integrated.
        //
        //    Input, double[] xval, the left and right X coordinates.
        //
        //    Input, double[] yval, the lower and upper Y coordinates.
        //
        //    Input, int[] nsub, the number of subintervals to use in the X 
        //    and Y directions.
        //
        //    Input, int ORDER, the order of the rule.
        //
        //    Input, double XTAB[ORDER], YTAB[ORDER], the abscissas.
        //
        //    Input, double WEIGHT[ORDER], the weights of the rule.
        //
        //    Output, double RECTANGLE_SUB_2D, the approximate integral of the function.
        //
    {
        double[] a = new double[2];
        double[] b = new double[2];
        int i;
        double result;

        a[0] = xval[0];
        a[1] = yval[0];
        b[0] = xval[1];
        b[1] = yval[1];

        for (i = 0; i < 2; i++)
        {
            if (!(Math.Abs(a[i] - b[i]) <= typeMethods.r8_epsilon()))
            {
                continue;
            }

            result = 0.0;
            return result;
        }

        for (i = 0; i < 2; i++)
        {
            switch (nsub[i])
            {
                case < 1:
                    Console.WriteLine("");
                    Console.WriteLine("RECTANGLE_SUB_2D - Fatal error!");
                    Console.WriteLine("  Nonpositive value of NSUB[" + i
                                                                     + "] = " + nsub[i] + "");
                    return 1;
            }
        }

        //
        //  Break up the X interval into NSUB(1) subintervals.
        //
        result = 0.0;

        for (i = 1; i <= nsub[0]; i++)
        {
            double xlo = typeMethods.r8vec_even_select(nsub[0] + 1, a[0], b[0], i);
            double xhi = typeMethods.r8vec_even_select(nsub[0] + 1, a[0], b[0], i + 1);
            //
            //  Break up the Y interval into NSUB(2) subintervals.
            //
            int j;
            for (j = 1; j <= nsub[1]; j++)
            {
                double ylo = typeMethods.r8vec_even_select(nsub[1] + 1, a[1], b[1], j);
                double yhi = typeMethods.r8vec_even_select(nsub[1] + 1, a[1], b[1], j + 1);

                double quad_sub = 0.0;
                int k;
                for (k = 0; k < order; k++)
                {
                    double x = xlo + 0.5 * (xtab[k] + 1.0) * (xhi - xlo);
                    double y = ylo + 0.5 * (ytab[k] + 1.0) * (yhi - ylo);

                    quad_sub += weight[k] * func(setting, x, y) / 4.0;
                }

                double volume_sub = (xhi - xlo) * (yhi - ylo);
                double result_sub = quad_sub * volume_sub;

                result += result_sub;
            }
        }

        return result;
    }

}