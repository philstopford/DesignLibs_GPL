using Burkardt.Types;

namespace Burkardt.ParabolaNS
{
    public static class Geometry
    {
        public static int parabola_ex(double x1, double y1, double x2, double y2, double x3,
                double y3, ref double x, ref double y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PARABOLA_EX finds the extremal point of a parabola determined by three points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 April 1999
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X1, Y1, X2, Y2, X3, Y3, the coordinates of three points
            //    on the parabola.  X1, X2 and X3 must be distinct.
            //
            //    Output, double *X, *Y, the X coordinate of the extremal point of the
            //    parabola, and the value of the parabola at that point.
            //
            //    Output, int PARABOLA_EX, error flag.
            //    0, no error.
            //    1, two of the X values are equal.
            //    2, the data lies on a straight line; there is no finite extremal
            //    point.
            //    3, the data lies on a horizontal line; every point is "extremal".
            //
        {
            double bot;

            x = 0.0;
            y = 0.0;

            if (x1 == x2 || x2 == x3 || x3 == x1)
            {
                return 1;
            }

            if (y1 == y2 && y2 == y3 && y3 == y1)
            {
                x = x1;
                y = y1;
                return 3;
            }

            bot = (x2 - x3) * y1 - (x1 - x3) * y2 + (x1 - x2) * y3;

            if (bot == 0.0)
            {
                return 2;
            }

            x = 0.5 * (
                x1 * x1 * (y3 - y2)
                + x2 * x2 * (y1 - y3)
                + x3 * x3 * (y2 - y1)) / bot;

            y = (
                    (x - x2) * (x - x3) * (x2 - x3) * y1
                    - (x - x1) * (x - x3) * (x1 - x3) * y2
                    + (x - x1) * (x - x2) * (x1 - x2) * y3) /
                ((x1 - x2) * (x2 - x3) * (x1 - x3));

            return 0;
        }

        public static int parabola_ex2(double x1, double y1, double x2, double y2, double x3,
                double y3, ref double x, ref double y, ref double a, ref double b, ref double c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PARABOLA_EX2 finds the extremal point of a parabola determined by three points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 May 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X1, Y1, X2, Y2, X3, Y3, the coordinates of three points
            //    on the parabola.  X1, X2 and X3 must be distinct.
            //
            //    Output, double *X, *Y, the X coordinate of the extremal point of the
            //    parabola, and the value of the parabola at that point.
            //
            //    Output, double *A, *B, *C, the coefficients that define the parabola:
            //    P(X) = A * X * X + B * X + C.
            //
            //    Output, int PARABOLA_EX2, error flag.
            //    0, no error.
            //    1, two of the X values are equal.
            //    2, the data lies on a straight line; there is no finite extremal
            //    point.
            //    3, the data lies on a horizontal line; any point is an "extremal point".
            //
        {
            double[] v = new double[3 * 3];
            double[] w;

            a = 0.0;
            b = 0.0;
            c = 0.0;
            x = 0.0;
            y = 0.0;

            if (x1 == x2 || x2 == x3 || x3 == x1)
            {
                return 1;
            }

            if (y1 == y2 && y2 == y3 && y3 == y1)
            {
                x = x1;
                y = y1;
                return 3;
            }

            //
            //  Set up the Vandermonde matrix.
            //
            v[0 + 0 * 3] = 1.0;
            v[1 + 0 * 3] = 1.0;
            v[2 + 0 * 3] = 1.0;

            v[0 + 1 * 3] = x1;
            v[1 + 1 * 3] = x2;
            v[2 + 1 * 3] = x3;

            v[0 + 2 * 3] = x1 * x1;
            v[1 + 2 * 3] = x2 * x2;
            v[2 + 2 * 3] = x3 * x3;
            //
            //  Get the inverse.
            //
            w = typeMethods.r8mat_inverse_3d(v);
            //
            //  Compute the parabolic coefficients.
            //
            c = w[0 + 0 * 3] * y1 + w[0 + 1 * 3] * y2 + w[0 + 2 * 3] * y3;
            b = w[1 + 0 * 3] * y1 + w[1 + 1 * 3] * y2 + w[1 + 2 * 3] * y3;
            a = w[2 + 0 * 3] * y1 + w[2 + 1 * 3] * y2 + w[2 + 2 * 3] * y3;

            //
            //  Determine the extremal point.
            //
            if (a == 0.0)
            {
                return 2;
            }

            x = -b / (2.0 * a);
            y = a * x * x + b * x + c;

            return 0;
        }

    }
}