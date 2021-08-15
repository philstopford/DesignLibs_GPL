using System;

namespace Burkardt.Stroud
{
    public static class Polygon
    {
        public static double polygon_1_2d(int n, double[] x, double[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_1_2D integrates the function 1 over a polygon in 2D.
            //
            //  Integration region:
            //
            //    The polygon bounded by the points (X(1:N), Y(1:N)).
            //
            //  Formula:
            //
            //    INTEGRAL = 0.5 * sum ( 1 <= I <= N )
            //      ( X(I) + X(I-1) ) * ( Y(I) - Y(I-1) )
            //
            //    where X(0) and Y(0) should be replaced by X(N) and Y(N).
            //
            //  Discussion:
            //
            //    The integral of 1 over a polygon is the area of the polygon.
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
            //    SF Bockman,
            //    Generalizing the Formula for Areas of Polygons to Moments,
            //    American Mathematical Society Monthly,
            //    1989, pages 131-132.
            //
            //  Parameters:
            //
            //    Input, int N, the number of vertices of the polygon.
            //    N should be at least 3 for a nonzero result.
            //
            //    Input, double X[N], Y[N], the coordinates of the vertices
            //    of the polygon.  These vertices should be given in counter-clockwise order.
            //
            //    Output, double POLYGON_1_2D, the value of the integral.
            //
        {
            int i;
            int im1;
            double result;

            result = 0.0;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("POLYGON_1_2D - Fatal error!");
                Console.WriteLine("  The number of vertices must be at least 3.");
                Console.WriteLine("  The input value of N = " + n + "");
                return (1);
            }

            for (i = 0; i < n; i++)
            {
                if (i == 0)
                {
                    im1 = n - 1;
                }
                else
                {
                    im1 = i - 1;
                }

                result = result + 0.5 * (x[i] + x[im1]) * (y[i] - y[im1]);
            }

            return result;
        }

        public static double polygon_x_2d(int n, double[] x, double[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_X_2D integrates the function X over a polygon in 2D.
            //
            //  Integration region:
            //
            //    The polygon bounded by the points (X(1:N), Y(1:N)).
            //
            //  Formula:
            //
            //    INTEGRAL = (1/6) * sum ( 1 <= I <= N )
            //      ( X(I)^2 + X(I) * X(I-1) + X(I-1)^2 ) * ( Y(I) - Y(I-1) )
            //
            //    where X(0) and Y(0) should be replaced by X(N) and Y(N).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    SF Bockman,
            //    Generalizing the Formula for Areas of Polygons to Moments,
            //    American Mathematical Society Monthly,
            //    1989, pages 131-132.
            //
            //  Parameters:
            //
            //    Input, int N, the number of vertices of the polygon.
            //    N should be at least 3 for a nonzero result.
            //
            //    Input, double X[N], Y[N], the coordinates of the vertices
            //    of the polygon.  These vertices should be given in counter-clockwise order.
            //
            //    Output, double POLYGON_X_2D, the value of the integral.
            //
        {
            int i;
            int im1;
            double result;

            result = 0.0;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("POLYGON_X_2D - Fatal error!");
                Console.WriteLine("  The number of vertices must be at least 3.");
                Console.WriteLine("  The input value of N = " + n + "");
                return (1);
            }

            for (i = 0; i < n; i++)
            {
                if (i == 0)
                {
                    im1 = n - 1;
                }
                else
                {
                    im1 = i - 1;
                }

                result = result + (x[i] * x[i] + x[i] * x[im1] + x[im1] * x[im1])
                    * (y[i] - y[im1]);
            }

            result = result / 6.0;

            return result;
        }

        public static double polygon_xx_2d(int n, double[] x, double[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_XX_2D integrates the function X*X over a polygon in 2D.
            //
            //  Integration region:
            //
            //    The polygon bounded by the points (X(1:N), Y(1:N)).
            //
            //  Formula:
            //
            //    INTEGRAL = (1/12) * sum ( 1 <= I <= N )
            //      ( X(I)^3 + X(I)^2 * X(I-1) + X(I) * X(I-1)^2 + X(I-1)^3 )
            //      * ( Y(I) - Y(I-1) )
            //
            //    where X(0) and Y(0) should be replaced by X(N) and Y(N).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    SF Bockman,
            //    Generalizing the Formula for Areas of Polygons to Moments,
            //    American Mathematical Society Monthly,
            //    Volume 96, Number 2, February 1989, pages 131-132.
            //
            //  Parameters:
            //
            //    Input, int N, the number of vertices of the polygon.
            //    N should be at least 3 for a nonzero result.
            //
            //    Input, double X[N], Y[N], the coordinates of the vertices
            //    of the polygon.  These vertices should be given in
            //    counter-clockwise order.
            //
            //    Output, double RESULT, the value of the integral.
            //
        {
            int i;
            int im1;
            double result;

            result = 0.0;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("POLYGON_XX_2D - Fatal error!");
                Console.WriteLine("  The number of vertices must be at least 3.");
                Console.WriteLine("  The input value of N = " + n + "");
                return (1);
            }

            for (i = 0; i < n; i++)
            {
                if (i == 0)
                {
                    im1 = n - 1;
                }
                else
                {
                    im1 = i - 1;
                }

                result = result + (x[i] * x[i] * x[i] + x[i] * x[i] * x[im1]
                                                      + x[i] * x[im1] * x[im1] + x[im1] * x[im1] * x[im1])
                    * (y[i] - y[im1]);
            }

            result = result / 12.0;

            return result;
        }

        public static double polygon_xy_2d(int n, double[] x, double[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_XY_2D integrates the function X*Y over a polygon in 2D.
            //
            //  Integration region:
            //
            //    The polygon bounded by the points (X(1:N), Y(1:N)).
            //
            //  Formula:
            //
            //    INTEGRAL = (1/24) * sum ( 1 <= I <= N )
            //      ( Y(I)   * ( 3 * X(I)**2 + 2 * X(I) * X(I-1) +     X(I-1)**2 )
            //      + Y(I-1) * (     X(I)**2 + 2 * X(I) * X(I-1) + 3 * X(I-1)**2 ) )
            //      * ( Y(I) - Y(I-1) )
            //
            //    where X(0) and Y(0) should be replaced by X(N) and Y(N).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    SF Bockman,
            //    Generalizing the Formula for Areas of Polygons to Moments,
            //    American Mathematical Society Monthly,
            //    Volume 96, Number 2, February 1989, pages 131-132.
            //
            //  Parameters:
            //
            //    Input, int N, the number of vertices of the polygon.
            //    N should be at least 3 for a nonzero result.
            //
            //    Input, double X[N], Y[N], the coordinates of the vertices
            //    of the polygon.  These vertices should be given in
            //    counter-clockwise order.
            //
            //    Output, double POLYGON_XY_2D, the value of the integral.
            //
        {
            int i;
            int im1;
            double result;

            result = 0.0;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("POLYGON_XY_2D - Fatal error!");
                Console.WriteLine("  The number of vertices must be at least 3.");
                Console.WriteLine("  The input value of N = " + n + "");
                return (1);
            }

            for (i = 0; i < n; i++)
            {
                if (i == 0)
                {
                    im1 = n - 1;
                }
                else
                {
                    im1 = i - 1;
                }

                result = result + (
                    y[i] * (3.0 * x[i] * x[i] + 2.0 * x[i] * x[im1] + x[im1] * x[im1])
                    + y[im1] * (x[i] * x[i] + 2.0 * x[i] * x[im1] + 3.0 * x[im1] * x[im1])
                ) * (y[i] - y[im1]);
            }

            result = result / 24.0;

            return result;
        }

        public static double polygon_y_2d(int n, double[] x, double[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_Y_2D integrates the function Y over a polygon in 2D.
            //
            //  Integration region:
            //
            //    The polygon bounded by the points (X(1:N), Y(1:N)).
            //
            //  Formula:
            //
            //    INTEGRAL = (1/6) * sum ( 1 <= I <= N )
            //      - ( Y(I)^2 + Y(I) * Y(I-1) + Y(I-1)^2 ) * ( X(I) - X(I-1) )
            //
            //    where X(0) and Y(0) should be replaced by X(N) and Y(N).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    SF Bockman,
            //    Generalizing the Formula for Areas of Polygons to Moments,
            //    American Mathematical Society Monthly,
            //    Volume 96, Number 2, February 1989, pages 131-132.
            //
            //  Parameters:
            //
            //    Input, int N, the number of vertices of the polygon.
            //    N should be at least 3 for a nonzero result.
            //
            //    Input, double X[N], Y[N], the coordinates of the vertices
            //    of the polygon.  These vertices should be given in
            //    counter-clockwise order.
            //
            //    Output, double POLYGON_Y_2D, the value of the integral.
            //
        {
            int i;
            int im1;
            double result;

            result = 0.0;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("POLYGON_Y_2D - Fatal error!");
                Console.WriteLine("  The number of vertices must be at least 3.");
                Console.WriteLine("  The input value of N = " + n + "");
                return (1);
            }

            for (i = 0; i < n; i++)
            {
                if (i == 0)
                {
                    im1 = n - 1;
                }
                else
                {
                    im1 = i - 1;
                }

                result = result - (y[i] * y[i] + y[i] * y[im1] + y[im1] * y[im1])
                    * (x[i] - x[im1]);
            }

            result = result / 6.0;

            return result;
        }

        public static double polygon_yy_2d(int n, double[] x, double[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_YY_2D integrates the function Y*Y over a polygon in 2D.
            //
            //  Integration region:
            //
            //    The polygon bounded by the points (X(1:N), Y(1:N)).
            //
            //  Formula:
            //
            //    INTEGRAL = (1/12) * sum ( 1 <= I <= N )
            //      - ( Y(I)^3 + Y(I)^2 * Y(I-1) + Y(I) * Y(I-1)^2 + Y(I-1)^3 )
            //      * ( X(I) - X(I-1) )
            //
            //    where X(0) and Y(0) should be replaced by X(N) and Y(N).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    SF Bockman,
            //    Generalizing the Formula for Areas of Polygons to Moments,
            //    American Mathematical Society Monthly,
            //    Volume 96, Number 2, February 1989, pages 131-132.
            //
            //  Parameters:
            //
            //    Input, int N, the number of vertices of the polygon.
            //    N should be at least 3 for a nonzero result.
            //
            //    Input, double X[N], Y[N], the coordinates of the vertices
            //    of the polygon.  These vertices should be given in
            //    counter-clockwise order.
            //
            //    Output, double POLYGON_YY_2D, the value of the integral.
            //
        {
            int i;
            int im1;
            double result;

            result = 0.0;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("POLYGON_YY_2D - Fatal error!");
                Console.WriteLine("  The number of polygonal vertices must be");
                Console.WriteLine("  at least 3, but the input polygon has N = " + n + "");
                return (1);
            }

            for (i = 0; i < n; i++)
            {
                if (i == 0)
                {
                    im1 = n - 1;
                }
                else
                {
                    im1 = i - 1;
                }

                result = result - (
                    y[i] * y[i] * y[i]
                    + y[i] * y[i] * y[im1]
                    + y[i] * y[im1] * y[im1]
                    + y[im1] * y[im1] * y[im1]
                ) * (x[i] - x[im1]);
            }

            result = result / 12.0;

            return result;
        }
   }
}