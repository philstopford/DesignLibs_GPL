using System;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace Burkardt.Stroud
{
    public static class Circle
    {
        public static double circle_annulus(Func<double, double, double> func, double[] center,
                double radius1, double radius2, int nr)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CIRCLE_ANNULUS approximates an integral in an annulus.
            //
            //  Integration region:
            //
            //    RADIUS1^2 <= ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= RADIUS2^2
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    William Peirce,
            //    Numerical Integration Over the Planar Annulus,
            //    Journal of the Society for Industrial and Applied Mathematics,
            //    Volume 5, Number 2, June 1957, pages 66-73.
            //
            //  Parameters:
            //
            //    Input, double FUNC ( double x, double y ), the name of the user supplied 
            //    function of two variables which is to be integrated.
            //
            //    Input, double CENTER[2], the center of the circle.
            //
            //    Input, double RADIUS1, RADIUS2, the radii of the circles.
            //
            //    Input, int NR, the order of the rule.  This quantity specifies
            //    the number of distinct radii to use.  The number of angles used will
            //    be 4*NR, for a total of 4*NR*NR points.
            //
            //    Output, double CIRCLE_ANNULUS, the approximation to the integral.
            //
        {
            double a;
            double area;
            double b;
            double c;
            double d;
            int i;
            int j;
            int nt;
            double pi = 3.141592653589793;
            double quad;
            double[] ra;
            double result;
            double[] rw;
            double t;
            double tw;
            double x;
            double y;
            //
            //  Choose radial abscissas and weights.
            //
            ra = new double[nr];
            rw = new double[nr];

            LegendreQuadrature.legendre_set(nr, ref ra, ref rw);
            a = -1.0;
            b = +1.0;
            c = radius1 * radius1;
            d = radius2 * radius2;

            QuadratureRule.rule_adjust(a, b, c, d, nr, ref ra, ref rw);

            for (i = 0; i < nr; i++)
            {
                ra[i] = Math.Sqrt(ra[i]);
            }

            for (i = 0; i < nr; i++)
            {
                rw[i] = rw[i] / (radius2 - radius1) / (radius2 + radius1);
            }

            //
            //  Set angular abscissas and weights.
            //
            nt = 4 * nr;

            tw = 1.0 / (double)(nt);
            //
            //  Approximate the integral.
            //
            quad = 0.0;
            for (i = 0; i < nt; i++)
            {
                t = 2.0 * pi * (double)(i - 1) / (double)(nt);
                for (j = 0; j < nr; j++)
                {
                    x = center[0] + ra[j] * Math.Cos(t);
                    y = center[1] + ra[j] * Math.Sin(t);
                    quad = quad + tw * rw[j] * func(x, y);
                }
            }

            area = circle_annulus_area_2d(radius1, radius2);
            result = quad * area;

            return result;
        }

        public static double circle_annulus_area_2d(double radius1, double radius2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CIRCLE_ANNULUS_AREA_2D returns the area of a circular annulus in 2D.
            //
            //  Integration region:
            //
            //    RADIUS1^2 <= ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= RADIUS2^2
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double RADIUS1, RADIUS2, the radii of the circles.
            //
            //    Output, double CIRCLE_ANNULUS_AREA_2D, the area of the annulus.
            //
        {
            double pi = 3.141592653589793;
            double value;

            value = pi * (radius1 + radius2) * (radius2 - radius1);

            return value;
        }

        public static double circle_annulus_sector(Func<double, double, double> func,
                double[] center, double radius1, double radius2, double theta1,
                double theta2, int nr)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CIRCLE_ANNULUS_SECTOR approximates an integral in a circular annulus sector.
            //
            //  Discussion:
            //
            //    A circular annulus sector comprises the area between two concentric
            //    circles and two concentric rays.
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
            //    William Peirce,
            //    Numerical Integration Over the Planar Annulus,
            //    Journal of the Society for Industrial and Applied Mathematics,
            //    Volume 5, Number 2, June 1957, pages 66-73.
            //
            //  Parameters:
            //
            //    Input, double FUNC ( double x, double y ), the name of the 
            //    user supplied function to be integrated.
            //
            //    Input, double CENTER[2], the center of the circle.
            //
            //    Input, double RADIUS1, RADIUS2, the radii of the circles.
            //
            //    Input, double THETA1, THETA2, the angles defining the sector.
            //    The sector is measured from THETA1 to THETA2.
            //
            //    Input, int NR, the order of the rule.  This quantity specifies
            //    the number of distinct radii to use.  The number of angles used will
            //    be 4*NR, for a total of 4*NR*NR points.
            //
            //    Output, double CIRCLE_ANNULUS_SECTOR, the approximation to the integral.
            //
        {
            double a;
            double area;
            double b;
            double c;
            double d;
            int i;
            int j;
            int nt;
            double quad;
            double[] ra;
            double result;
            double[] rw;
            double[] ta;
            double[] tw;
            double x;
            double y;
            //
            //  Set the radial abscissas and weights.
            //
            ra = new double[nr];
            rw = new double[nr];

            LegendreQuadrature.legendre_set(nr, ref ra, ref rw);

            a = -1.0;
            b = +1.0;
            c = radius1 * radius1;
            d = radius2 * radius2;

            QuadratureRule.rule_adjust(a, b, c, d, nr, ref ra, ref rw);

            for (i = 0; i < nr; i++)
            {
                ra[i] = Math.Sqrt(ra[i]);
            }

            for (i = 0; i < nr; i++)
            {
                rw[i] = rw[i] / (radius2 - radius1) / (radius2 + radius1);
            }

            //
            //  Pick angles evenly spaced between THETA1 and THETA2, but do not
            //  include the endpoints, and use a half interval for the first and last.
            //
            nt = 4 * nr;

            ta = typeMethods.tvec_even_bracket3(nt, theta1, theta2);
            tw = new double[nt];
            for (i = 0; i < nt; i++)
            {
                tw[i] = 1.0 / (double)(nt);
            }

            //
            //  Approximate the integral.
            //
            quad = 0.0;
            for (i = 0; i < nt; i++)
            {
                for (j = 0; j < nr; j++)
                {
                    x = center[0] + ra[j] * Math.Cos(ta[i]);
                    y = center[1] + ra[j] * Math.Sin(ta[i]);
                    quad = quad + tw[i] * rw[j] * func(x, y);
                }
            }

            area = circle_annulus_sector_area_2d(radius1, radius2, theta1, theta2);

            result = quad * area;



            return result;
        }

        public static double circle_annulus_sector_area_2d(double radius1, double radius2,
                double theta1, double theta2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CIRCLE_ANNULUS_SECTOR_AREA_2D returns the area of a circular annulus sector in 2D.
            //
            //  Discussion:
            //
            //    A circular annulus sector comprises the area between two concentric
            //    circles and two concentric rays.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double RADIUS1, RADIUS2, the radii of the circles.
            //
            //    Input, double THETA1, THETA2, the angles of the rays.
            //    Ordinarily, (THETA2-THETA1) is between 0 and 2*PI.
            //
            //    Output, double CIRCLE_ANNULUS_SECTOR_AREA_2D, the area of the
            //    circular annulus sector.
            //
        {
            double area;

            area = 0.5 * (radius1 + radius2) * (radius2 - radius1)
                   * (theta2 - theta1);

            return area;
        }

        public static double circle_area_2d(double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CIRCLE_AREA_2D returns the area of a circle in 2D.
            //
            //  Integration region:
            //
            //    ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= R * R
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the radius of the circle.
            //
            //    Output, double CIRCLE_AREA_2D, the area of the circle.
            //
        {
            double area;
            double pi = 3.141592653589793;

            area = pi * r * r;

            return area;
        }

        public static double circle_cap_area_2d(double r, double h)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CIRCLE_CAP_AREA_2D computes the area of a circle cap in 2D.
            //
            //  Discussion:
            //
            //    Draw any radius R of the circle and denote as P the point where the
            //    radius intersects the circle.  Now consider the point Q which lies
            //    on the radius and which is H units from P.  The line which is
            //    perpendicular to the radius R and passes through Q divides the
            //    circle into two pieces.  The piece including the point P is the
            //    circular cap of height (or thickness) H.
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
            //  Parameters:
            //
            //    Input, double R, the radius of the circle.
            //
            //    Input, double H, the "height" of the circle cap.  
            //
            //    Output, double CIRCLE_CAP_AREA_2D, the area of the circle cap.
            //
        {
            double area = 0;
            double pi = 3.141592653589793;
            double theta;

            if (h <= 0.0)
            {
                area = 0.0;
            }
            else if (h <= r)
            {
                theta = 2.0 * Math.Asin(Math.Sqrt(h * (2.0 * r - h)) / r);
                area = r * r * (theta - Math.Sin(theta)) / 2.0;
            }
            else if (h <= 2.0 * r)
            {
                theta = 2.0 * Math.Asin(Math.Sqrt(h * (2.0 * r - h)) / r);
                area = r * r * (pi - (theta - Math.Sin(theta)) / 2.0);
            }
            else if (2.0 * r <= h)
            {
                area = pi * r * r;
            }

            return area;
        }

        public static double circle_cum(Func<double, double, double> func, double[] center,
                double radius, int order)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CIRCLE_CUM approximates an integral on the circumference of a circle in 2D.
            //
            //  Integration region:
            //
            //    ( X - CENTER(1) )^2 + ( Y - CENTER(2) )^2 <= R * R
            //
            //  Discussion:
            //
            //    An ORDER point, (ORDER-1)-th degree formula is used, 
            //    Stroud number U2:M-1.
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
            //    Input, double FUNC ( double x, double y ), the name of the 
            //    user supplied function to be integrated.
            //
            //    Input, double CENTER[2], the coordinates of the center of 
            //    the circle.
            //
            //    Input, double RADIUS, the radius of the circle.
            //
            //    Input, int ORDER, the number of points to use.
            //
            //    Output, double CIRCLE_CUM, the approximate integral of the function.
            //
        {
            double angle;
            int i;
            double pi = 3.141592653589793;
            double quad;
            double result;
            double volume;
            double x;
            double y;

            quad = 0.0;

            for (i = 0; i < order; i++)
            {
                angle = (double)(2 * i) * pi / (double)(order);
                x = center[0] + radius * Math.Cos(angle);
                y = center[1] + radius * Math.Sin(angle);
                quad = quad + func(x, y);
            }

            quad = quad / (double)(order);

            volume = pi * radius * radius;
            result = quad * volume;

            return result;
        }

        public static void circle_rt_set(int rule, int nr, int nt, int nc, ref double[] ra,
                ref double[] rw, ref double[] ta, ref double[] tw, ref double cw)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CIRCLE_RT_SET sets an R, THETA product quadrature rule in the unit circle.
            //
            //  Discussion:
            //
            //    For a given value of RULE, here are the number of points used at the
            //    center (NC), the number of points along the radial direction (NR) and
            //    the number of points along the theta direction (NT).  The total number
            //    of points in the rule will be 
            //
            //      Total = NC + NR * NT.
            //
            //    The user, when choosing RULE, must allocate enough space in the arrays
            //    RA, RW, TA and TW for the resulting values of NR and NT.
            //
            //    RULE  NC  NR  NT  Total
            //    ----  --  --  --  -----
            //       1   1   0   0      1
            //       2   0   1   4      4
            //       3   1   1   4      5
            //       4   1   1   6      7
            //       5   1   2   4      9
            //       6   0   3   4     12
            //       7   1   2  10     21
            //       8   0   4  16     64
            //       9   0   5  20    120
            //
            //    The integral of F(X,Y) over the unit circle is approximated by
            //
            //      Integral ( X*X + Y*Y <= 1 ) F(X,Y) dx dy 
            //      = Integral ( 0 <= R <= 1, 0 <= T <= 2PI ) F(R*cos(T),R*sin(T)) r dr dt
            //      = approximately
            //        CW * F(0,0) 
            //        + sum ( 1 <= I <= NR ) Sum ( 1 <= J <= NT )
            //        RW(I) * TW(J) * F ( R(I) * Math.Cos ( TA(J) ), R(I) * Math.Sin ( TA(J) ) )
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
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, int RULE, the rule desired.
            //
            //    Input, int NR, the number of R abscissas.
            //
            //    Input, int NT, the number of Theta abscissas.
            //
            //    Input, int NC, the number of center abscissas (0 or 1 ).
            //
            //    Output, double RA[NR], RW[NR], the R abscissas and weights.
            //
            //    Output, double TA[NT], TW[NT], the THETA abscissas and weights.
            //
            //    Output, double *ZW, the weight to use for the center.
            //
        {
            double a;
            double b;
            double c;
            double d;
            int i;
            double pi = 3.141592653589793;
            double u;
            double v;
            double w;

            if (rule == 1)
            {
                cw = 1.0;
            }
            else if (rule == 2)
            {
                ra[0] = 0.5;
                rw[0] = 1.0;

                for (i = 0; i < nt; i++)
                {
                    ta[i] = (double)(2 * i) * pi / (double)(nt);
                }

                for (i = 0; i < nt; i++)
                {
                    tw[i] = 1.0 / (double)(nt);
                }

                cw = 0.0;
            }
            else if (rule == 3)
            {
                ra[0] = 1.0;
                rw[0] = 1.0;

                for (i = 0; i < nt; i++)
                {
                    ta[i] = (double)(2 * i) * pi / (double)(nt);
                }

                for (i = 0; i < nt; i++)
                {
                    tw[i] = 0.125;
                }

                cw = 0.5;
            }
            else if (rule == 4)
            {
                ra[0] = Math.Sqrt(2.0 / 3.0);
                rw[0] = 1.0;

                for (i = 0; i < nt; i++)
                {
                    ta[i] = (double)(2 * i) * pi / (double)(nt);
                }

                for (i = 0; i < nt; i++)
                {
                    tw[i] = 0.125;
                }

                cw = 0.25;
            }
            else if (rule == 5)
            {
                a = 1.0;
                b = Math.Sqrt(2.0) / 2.0;
                u = 1.0 / 6.0;
                v = 4.0 / 6.0;

                ra[0] = a;
                ra[1] = b;
                rw[0] = u;
                rw[1] = v;

                for (i = 0; i < nt; i++)
                {
                    ta[i] = (double)(2 * i) * pi / (double)(nt);
                }

                for (i = 0; i < nt; i++)
                {
                    tw[i] = 1.0 / (double)(nt);
                }

                cw = 4.0 / 24.0;
            }
            else if (rule == 6)
            {
                a = Math.Sqrt(3.0) / 2.0;
                b = Math.Sqrt((27.0 - 3.0 * Math.Sqrt(29.0)) / 52.0);
                c = Math.Sqrt((27.0 + 3.0 * Math.Sqrt(29.0)) / 52.0);

                u = 8.0 / 27.0;
                v = (551.0 + 41.0 * Math.Sqrt(29.0)) / 1566.0;
                w = (551.0 - 41.0 * Math.Sqrt(29.0)) / 1566.0;

                ra[0] = a;
                ra[1] = b;
                ra[2] = c;
                rw[0] = u;
                rw[1] = v;
                rw[2] = w;

                for (i = 0; i < nt; i++)
                {
                    ta[i] = (double)(2 * i) * pi / (double)(nt);
                }

                for (i = 0; i < nt; i++)
                {
                    tw[i] = 1.0 / (double)(nt);
                }

                cw = 0.0;
            }
            else if (rule == 7)
            {
                a = Math.Sqrt((6.0 - Math.Sqrt(6.0)) / 10.0);
                b = Math.Sqrt((6.0 + Math.Sqrt(6.0)) / 10.0);
                u = (16.0 + Math.Sqrt(6.0)) / 36.0;
                v = (16.0 - Math.Sqrt(6.0)) / 36.0;

                ra[0] = a;
                ra[1] = b;
                rw[0] = u;
                rw[1] = v;

                for (i = 0; i < nt; i++)
                {
                    ta[i] = (double)(2 * i) * pi / (double)(nt);
                }

                for (i = 0; i < nt; i++)
                {
                    tw[i] = 1.0 / (double)(nt);
                }

                cw = 1.0 / 9.0;
            }
            else if (rule == 8)
            {
                LegendreQuadrature.legendre_set(nr, ref ra, ref rw);
                a = -1.0;
                b = +1.0;
                c = 0.0;
                d = +1.0;
                QuadratureRule.rule_adjust(a, b, c, d, nr, ref ra, ref rw);

                for (i = 0; i < nr; i++)
                {
                    ra[i] = Math.Sqrt(ra[i]);
                }

                for (i = 0; i < nt; i++)
                {
                    ta[i] = (double)(2 * i) * pi / (double)(nt);
                }

                for (i = 0; i < nt; i++)
                {
                    tw[i] = 1.0 / (double)(nt);
                }

                cw = 0.0;
            }
            else if (rule == 9)
            {
                LegendreQuadrature.legendre_set(nr, ref ra, ref rw);
                a = -1.0;
                b = +1.0;
                c = 0.0;
                d = +1.0;
                QuadratureRule.rule_adjust(a, b, c, d, nr, ref ra, ref rw);

                for (i = 0; i < nr; i++)
                {
                    ra[i] = Math.Sqrt(ra[i]);
                }

                for (i = 0; i < nt; i++)
                {
                    ta[i] = (double)(2 * i) * pi / (double)(nt);
                }

                for (i = 0; i < nt; i++)
                {
                    tw[i] = 1.0 / (double)(nt);
                }

                cw = 0.0;
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("CIRCLE_RT_SET - Fatal error!");
                Console.WriteLine("  There is no rule of index " + rule + "");
            }
        }

        public static void circle_rt_size(int rule, ref int nr, ref int nt, ref int nc)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CIRCLE_RT_SIZE sizes an R, THETA product quadrature rule in the unit circle.
            //
            //  Discussion:
            //
            //    For a given value of RULE, here are the number of points used at the
            //    center (NC), the number of points along the radial direction (NR) and
            //    the number of points along the theta direction (NT).  The total number
            //    of points in the rule will be 
            //
            //      Total = NC + NR * NT.
            //
            //    The user, when choosing RULE, must allocate enough space in the arrays
            //    RA, RW, TA and TW for the resulting values of NR and NT.
            //
            //    RULE  NC  NR  NT  Total
            //    ----  --  --  --  -----
            //       1   1   0   0      1
            //       2   0   1   4      4
            //       3   1   1   4      5
            //       4   1   1   6      7
            //       5   1   2   4      9
            //       6   0   3   4     12
            //       7   1   2  10     21
            //       8   0   4  16     64
            //       9   0   5  20    120
            //
            //    The integral of F(X,Y) over the unit circle is approximated by
            //
            //      Integral ( X*X + Y*Y <= 1 ) F(X,Y) dx dy 
            //      = Integral ( 0 <= R <= 1, 0 <= T <= 2PI ) F(R*cos(T),R*sin(T)) r dr dt
            //      = approximately
            //        ZW * F(0,0) 
            //        + sum ( 1 <= I <= NR ) Sum ( 1 <= J <= NT )
            //        RW(I) * TW(J) * F ( R(I) * Math.Cos ( TA(J) ), R(I) * Math.Sin ( TA(J) ) )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, int RULE, the rule desired.
            //
            //    Output, int *NR, the number of R abscissas.
            //    
            //    Output, int *NT, the number of Theta abscissas.
            //
            //    Output, int *NC, the number of center abscissas (0 or 1).
            //
        {
            if (rule == 1)
            {
                nr = 0;
                nt = 0;
                nc = 1;
            }
            else if (rule == 2)
            {
                nr = 1;
                nt = 4;
                nc = 0;
            }
            else if (rule == 3)
            {
                nr = 1;
                nt = 4;
                nc = 1;
            }
            else if (rule == 4)
            {
                nr = 1;
                nt = 6;
                nc = 1;
            }
            else if (rule == 5)
            {
                nr = 2;
                nt = 4;
                nc = 1;
            }
            else if (rule == 6)
            {
                nr = 3;
                nt = 4;
                nc = 0;
            }
            else if (rule == 7)
            {
                nr = 2;
                nt = 10;
                nc = 1;
            }
            else if (rule == 8)
            {
                nr = 4;
                nt = 16;
                nc = 0;
            }
            else if (rule == 9)
            {
                nr = 5;
                nt = 20;
                nc = 0;
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("CIRCLE_RT_SIZE - Fatal error!");
                Console.WriteLine("  There is no rule of index " + rule + "");
            }
        }

        public static double circle_rt_sum(Func<double, double, double> func, double[] center,
                double radius, int nr, double[] ra, double[] rw, int nt, double[] ta,
                double[] tw, double zw)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CIRCLE_RT_SUM applies an R, THETA product quadrature rule inside a circle.
            //
            //  Integration region:
            //
            //    (X-CENTER(1))^2 + (Y-CENTER(2))^2 <= RADIUS^2.
            //
            //  Discussion:
            //
            //    The product rule is assumed to be have the form:
            //
            //      Integral_Approx = ZW * F(CENTER(1),CENTER(2)) +
            //        sum ( 1 <= IR <= NR ) Sum ( 1 <= IT <= NT )
            //        RW(IR) * TW(IT) * F ( CENTER(1) + R(IR) * RADIUS * Cos ( TA(IT) ),
            //                              CENTER(2) + R(IR) * RADIUS * Sin ( TA(IT) ) )
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
            //  Parameters:
            //
            //    Input, double FUNC ( double x, double y ), the name of the 
            //    user supplied function to be integrated.
            //
            //    Input, double CENTER[2], the center of the circle.
            //
            //    Input, double RADIUS, the radius of the circle.
            //
            //    Input, int NR, the number of R abscissas.
            //
            //    Input, double RA[NR], RW[NR], the R abscissas and weights.
            //
            //    Input, int NT, the number of Theta abscissas.
            //
            //    Input, double TA[NT], TW[NT], the THETA abscissas and weights.
            //
            //    Input, double ZW, the weight to use for the center.
            //
            //    Output, double CIRCLE_RT_SUM, the approximate integral of the function.
            //
        {
            int ir;
            int it;
            double quad;
            double rct;
            double result;
            double rst;
            double volume;
            double x;
            double y;

            quad = 0.0;

            if (zw != 0.0)
            {
                x = center[0];
                y = center[1];
                quad = quad + zw * func(x, y);
            }

            for (it = 0; it < nt; it++)
            {
                rct = radius * Math.Cos(ta[it]);
                rst = radius * Math.Sin(ta[it]);
                for (ir = 0; ir < nr; ir++)
                {
                    x = center[0] + ra[ir] * rct;
                    y = center[1] + ra[ir] * rst;
                    quad = quad + tw[it] * rw[ir] * func(x, y);
                }
            }

            volume = circle_area_2d(radius);
            result = quad * volume;

            return result;
        }

        public static double circle_sector(Func<double, double, double> func, double[] center,
                double radius, double theta1, double theta2, int nr)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CIRCLE_SECTOR approximates an integral in a circular sector.
            //
            //  Discussion:
            //
            //    A sector is contained within a circular arc and the lines joining each
            //    endpoint of the arc to the center of the circle.
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
            //  Parameters:
            //
            //    Input, double FUNC ( double x, double y ), the name of the 
            //    user supplied function to be integrated.
            //
            //    Input, double CENTER[2], the center of the circle.
            //
            //    Input, double RADIUS, the radius of the circle.
            //
            //    Input, double THETA1, THETA2, the angles defining the sector.
            //    The sector is measured from THETA1 to THETA2.
            //
            //    Input, int NR, the number of radial values used in the approximation
            //    of the integral.  NR must be at least 1.  Higher values improve the
            //    accuracy of the integration, at the cost of more function evaluations.
            //
            //    Output, double CIRCLE_SECTOR, the approximation to the integral.
            //
        {
            double a;
            double area;
            double b;
            double c;
            double d;
            int i;
            int j;
            int nt;
            double quad;
            double[] ra;
            double result;
            double[] rw;
            double[] ta;
            double[] tw;
            double x;
            double y;
            //
            //  Set the radial abscissas and weights.
            //
            ra = new double[nr];
            rw = new double[nr];

            LegendreQuadrature.legendre_set(nr, ref ra, ref rw);

            a = -1.0;
            b = +1.0;
            c = 0.0;
            d = radius * radius;

            QuadratureRule.rule_adjust(a, b, c, d, nr, ref ra, ref rw);

            for (i = 0; i < nr; i++)
            {
                ra[i] = Math.Sqrt(ra[i]);
            }

            for (i = 0; i < nr; i++)
            {
                rw[i] = rw[i] / radius / radius;
            }

            //
            //  Pick angles evenly spaced between THETA1 and THETA2, but do not
            //  include the endpoints, and use a half interval for the first and last.
            //
            nt = 4 * nr;

            ta = typeMethods.tvec_even_bracket3(nt, theta1, theta2);

            tw = new double[nt];
            for (i = 0; i < nt; i++)
            {
                tw[i] = 1.0 / (double)(nt);
            }

            //
            //  Approximate the integral.
            //
            quad = 0.0;
            for (i = 0; i < nr; i++)
            {
                for (j = 0; j < nt; j++)
                {
                    x = center[0] + ra[i] * Math.Cos(ta[j]);
                    y = center[1] + ra[i] * Math.Sin(ta[j]);
                    quad = quad + rw[i] * tw[j] * func(x, y);
                }
            }

            area = circle_sector_area_2d(radius, theta1, theta2);
            result = quad * area;

            return result;
        }

        public static double circle_sector_area_2d(double r, double theta1, double theta2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CIRCLE_SECTOR_AREA_2D returns the area of a circular sector in 2D.
            //
            //  Discussion:
            //
            //    A sector is contained within a circular arc and the lines joining each
            //    endpoint of the arc to the center of the circle.
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
            //  Parameters:
            //
            //    Input, double R, the radius of the circle.
            //
            //    Input, double THETA1, THETA2, the angles of the rays
            //    that delimit the sector.
            //
            //    Output, double CIRCLE_SECTOR_AREA_2D, the area of the sector.
            //
        {
            double value;

            value = 0.5 * r * r * (theta2 - theta1);

            return value;
        }

        public static double circle_triangle_area_2d(double r, double theta1, double theta2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CIRCLE_TRIANGLE_AREA_2D returns the area of a circle triangle in 2D.
            //
            //  Discussion:
            //
            //    A circle triangle is formed by drawing a circular arc, and considering
            //    the triangle formed by the endpoints of the arc plus the center of
            //    the circle.
            //
            //    The normal situation is that 0 < ( THETA2 - THETA1 ) < PI.  Outside
            //    this range, the triangle can actually have NEGATIVE area.
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
            //  Parameters:
            //
            //    Input, double R, the radius of the circle.
            //
            //    Input, double THETA1, THETA2, the angles of the rays that
            //    delimit the arc.
            //
            //    Output, double CIRCLE_TRIANGLE_AREA_2D, the (signed) area
            //    of the triangle.
            //
        {
            double value;

            value = 0.5 * r * r * Math.Sin(theta2 - theta1);

            return value;
        }

        public static void circle_xy_set(int rule, int order, double[] xtab, double[] ytab,
                double[] weight)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CIRCLE_XY_SET sets an XY quadrature rule inside the unit circle in 2D.
            //
            //  Integration region:
            //
            //    X*X + Y*Y <= 1.0.
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
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //    Frank Lether,
            //    A Generalized Product Rule for the Circle,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 8, Number 2, June 1971, pages 249-253.
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, int RULE, the rule desired.
            //      1, 1 point 1-st degree;
            //      2, 4 point 3-rd degree, Stroud S2:3-1;
            //      3, 4 point 3-rd degree, Lether #1;
            //      4, 4 point 3-rd degree, Stroud S2:3-2;
            //      5, 5 point 3-rd degree;
            //      6, 7 point 5-th degree;
            //      7, 9 point 5-th degree;
            //      8, 9 point 5-th degree, Lether #2;
            //      9, 12 point 7-th degree;
            //     10, 16 point 7-th degree, Lether #3;
            //     11, 21 point 9-th degree, Stroud S2:9-3;
            //     12, 25 point 9-th degree, Lether #4 (after correcting error);
            //     13, 64 point 15-th degree Gauss product rule.
            //
            //    Input, int ORDER, the order of the desired rule.
            //
            //    Output, double XTAB[ORDER], YTAB[ORDER], the abscissas of 
            //    the rule.
            //
            //    Output, double WEIGHT[ORDER], the ORDER weights of the rule.
            //
        {
            double a;
            double b;
            double c;
            double d;
            double e;
            double f;
            double g;
            double h;
            int i;
            int j;
            int k;
            int nr;
            double pi = 3.141592653589793;
            double r;
            double[] ra;
            double[] rw;
            double s;
            double u;
            double v;
            double w;
            double w1;
            double w2;
            double w3;
            double w4;
            double w5;
            double w6;
            double w7;
            double w8;
            double w9;
            double z;

            if (rule == 1)
            {
                xtab[0] = 0.0;
                ytab[0] = 0.0;
                weight[0] = 1.0;
            }
            else if (rule == 2)
            {
                a = 0.5;
                b = 0.25;
                z = 0.0;

                xtab[0] = a;
                xtab[1] = -a;
                xtab[2] = z;
                xtab[3] = z;

                ytab[0] = z;
                ytab[1] = z;
                ytab[2] = a;
                ytab[3] = -a;

                weight[0] = b;
                weight[1] = b;
                weight[2] = b;
                weight[3] = b;
            }
            else if (rule == 3)
            {
                a = 0.5;
                b = 0.25;

                xtab[0] = a;
                xtab[1] = -a;
                xtab[2] = -a;
                xtab[3] = a;

                ytab[0] = a;
                ytab[1] = a;
                ytab[2] = -a;
                ytab[3] = -a;

                weight[0] = b;
                weight[1] = b;
                weight[2] = b;
                weight[3] = b;
            }
            else if (rule == 4)
            {
                a = Math.Sqrt(2.0) / 2.0;
                b = 0.25;

                xtab[0] = a;
                xtab[1] = -a;
                xtab[2] = -a;
                xtab[3] = a;

                ytab[0] = a;
                ytab[1] = a;
                ytab[2] = -a;
                ytab[3] = -a;

                weight[0] = b;
                weight[1] = b;
                weight[2] = b;
                weight[3] = b;
            }
            else if (rule == 5)
            {
                a = 1.0;
                b = 0.5;
                c = 0.125;
                z = 0.0;

                xtab[0] = z;
                xtab[1] = a;
                xtab[2] = z;
                xtab[3] = -a;
                xtab[4] = z;

                ytab[0] = z;
                ytab[1] = z;
                ytab[2] = a;
                ytab[3] = z;
                ytab[4] = -a;

                weight[0] = b;
                weight[1] = c;
                weight[2] = c;
                weight[3] = c;
                weight[4] = c;
            }
            else if (rule == 6)
            {
                a = Math.Sqrt(2.0 / 3.0);
                b = Math.Sqrt(1.0 / 6.0);
                c = Math.Sqrt(2.0) / 2.0;
                d = 0.125;
                e = 0.25;
                z = 0.0;

                xtab[0] = z;
                xtab[1] = a;
                xtab[2] = -a;
                xtab[3] = b;
                xtab[4] = -b;
                xtab[5] = b;
                xtab[6] = -b;

                ytab[0] = z;
                ytab[1] = z;
                ytab[2] = z;
                ytab[3] = c;
                ytab[4] = c;
                ytab[5] = -c;
                ytab[6] = -c;

                weight[0] = e;
                weight[1] = d;
                weight[2] = d;
                weight[3] = d;
                weight[4] = d;
                weight[5] = d;
                weight[6] = d;
            }
            else if (rule == 7)
            {
                a = 0.5;
                b = 1.0;
                c = 4.0 / 24.0;
                d = 1.0 / 24.0;
                z = 0.0;

                xtab[0] = z;
                xtab[1] = b;
                xtab[2] = -b;
                xtab[3] = z;
                xtab[4] = z;
                xtab[5] = a;
                xtab[6] = -a;
                xtab[7] = -a;
                xtab[8] = a;

                ytab[0] = z;
                ytab[1] = z;
                ytab[2] = z;
                ytab[3] = b;
                ytab[4] = -b;
                ytab[5] = a;
                ytab[6] = a;
                ytab[7] = -a;
                ytab[8] = -a;

                weight[0] = c;
                weight[1] = d;
                weight[2] = d;
                weight[3] = d;
                weight[4] = d;
                weight[5] = c;
                weight[6] = c;
                weight[7] = c;
                weight[8] = c;
            }
            else if (rule == 8)
            {
                a = Math.Sqrt(2.0) / 2.0;
                b = Math.Sqrt(3.0 / 5.0);
                c = Math.Sqrt(3.0 / 10.0);

                w1 = 16.0 / 72.0;
                w2 = 8.0 / 72.0;
                w3 = 10.0 / 72.0;
                w4 = 5.0 / 72.0;

                z = 0.0;

                xtab[0] = z;
                xtab[1] = a;
                xtab[2] = -a;
                xtab[3] = z;
                xtab[4] = z;
                xtab[5] = a;
                xtab[6] = a;
                xtab[7] = -a;
                xtab[8] = -a;

                ytab[0] = z;
                ytab[1] = z;
                ytab[2] = z;
                ytab[3] = b;
                ytab[4] = -b;
                ytab[5] = c;
                ytab[6] = -c;
                ytab[7] = c;
                ytab[8] = -c;

                weight[0] = w1;
                weight[1] = w2;
                weight[2] = w2;
                weight[3] = w3;
                weight[4] = w3;
                weight[5] = w4;
                weight[6] = w4;
                weight[7] = w4;
                weight[8] = w4;
            }
            else if (rule == 9)
            {
                a = Math.Sqrt(3.0) / 2.0;
                b = Math.Sqrt((27.0 - 3.0 * Math.Sqrt(29.0)) / 104.0);
                c = Math.Sqrt((27.0 + 3.0 * Math.Sqrt(29.0)) / 104.0);
                u = 2.0 / 27.0;
                v = (551.0 + 41.0 * Math.Sqrt(29.0)) / 6264.0;
                w = (551.0 - 41.0 * Math.Sqrt(29.0)) / 6264.0;
                z = 0.0;

                xtab[0] = a;
                xtab[1] = -a;
                xtab[2] = z;
                xtab[3] = z;
                xtab[4] = b;
                xtab[5] = -b;
                xtab[6] = b;
                xtab[7] = -b;
                xtab[8] = c;
                xtab[9] = c;
                xtab[10] = -c;
                xtab[11] = -c;

                ytab[0] = z;
                ytab[1] = z;
                ytab[2] = a;
                ytab[3] = -a;
                ytab[4] = b;
                ytab[5] = b;
                ytab[6] = -b;
                ytab[7] = -b;
                ytab[8] = c;
                ytab[9] = -c;
                ytab[10] = c;
                ytab[11] = -c;

                weight[0] = u;
                weight[1] = u;
                weight[2] = u;
                weight[3] = u;
                weight[4] = v;
                weight[5] = v;
                weight[6] = v;
                weight[7] = v;
                weight[8] = w;
                weight[9] = w;
                weight[10] = w;
                weight[11] = w;
            }
            else if (rule == 10)
            {
                a = Math.Sqrt((3.0 - Math.Sqrt(5.0)) / 8.0);
                b = Math.Sqrt((15.0 + 3.0 * Math.Sqrt(5.0)
                               - 2.0 * Math.Sqrt(30.0) - 2.0 * Math.Sqrt(6.0)) / 56.0);
                c = Math.Sqrt((15.0 + 3.0 * Math.Sqrt(5.0)
                                    + 2.0 * Math.Sqrt(30.0) + 2.0 * Math.Sqrt(6.0)) / 56.0);
                d = Math.Sqrt((3.0 + Math.Sqrt(5.0)) / 8.0);
                e = Math.Sqrt((15.0 - 3.0 * Math.Sqrt(5.0)
                                    - 2.0 * Math.Sqrt(30.0) + 2.0 * Math.Sqrt(6.0)) / 56.0);
                f = Math.Sqrt((15.0 - 3.0 * Math.Sqrt(5.0)
                    + 2.0 * Math.Sqrt(30.0) - 2.0 * Math.Sqrt(6.0)) / 56.0);
                w1 = (90.0 + 5.0 * Math.Sqrt(30.0) + 18.0 * Math.Sqrt(5.0)
                      + 5.0 * Math.Sqrt(6.0)) / 1440.0;
                w2 = (90.0 - 5.0 * Math.Sqrt(30.0) + 18.0 * Math.Sqrt(5.0)
                      - 5.0 * Math.Sqrt(6.0)) / 1440.0;
                w3 = (90.0 + 5.0 * Math.Sqrt(30.0) - 18.0 * Math.Sqrt(5.0)
                                                   - 5.0 * Math.Sqrt(6.0)) / 1440.0;
                w4 = (90.0 - 5.0 * Math.Sqrt(30.0) - 18.0 * Math.Sqrt(5.0)
                      + 5.0 * Math.Sqrt(6.0)) / 1440.0;

                xtab[0] = a;
                xtab[1] = a;
                xtab[2] = -a;
                xtab[3] = -a;
                xtab[4] = a;
                xtab[5] = a;
                xtab[6] = -a;
                xtab[7] = -a;
                xtab[8] = d;
                xtab[9] = d;
                xtab[10] = -d;
                xtab[11] = -d;
                xtab[12] = d;
                xtab[13] = d;
                xtab[14] = -d;
                xtab[15] = -d;

                ytab[0] = b;
                ytab[1] = -b;
                ytab[2] = b;
                ytab[3] = -b;
                ytab[4] = c;
                ytab[5] = -c;
                ytab[6] = c;
                ytab[7] = -c;
                ytab[8] = e;
                ytab[9] = -e;
                ytab[10] = e;
                ytab[11] = -e;
                ytab[12] = f;
                ytab[13] = -f;
                ytab[14] = f;
                ytab[15] = -f;

                weight[0] = w1;
                weight[1] = w1;
                weight[2] = w1;
                weight[3] = w1;
                weight[4] = w2;
                weight[5] = w2;
                weight[6] = w2;
                weight[7] = w2;
                weight[8] = w3;
                weight[9] = w3;
                weight[10] = w3;
                weight[11] = w3;
                weight[12] = w4;
                weight[13] = w4;
                weight[14] = w4;
                weight[15] = w4;
            }
            else if (rule == 11)
            {
                xtab[0] = 0.0;
                ytab[0] = 0.0;
                weight[0] = 1.0 / 9.0;

                for (i = 1; i < 11; i++)
                {
                    weight[i] = (16.0 + Math.Sqrt(6.0)) / 360.0;
                }

                for (i = 11; i < 21; i++)
                {
                    weight[i] = (16.0 - Math.Sqrt(6.0)) / 360.0;
                }

                r = Math.Sqrt((6.0 - Math.Sqrt(6.0)) / 10.0);

                for (i = 0; i < 10; i++)
                {
                    a = 2.0 * pi * (double)(i) / 10.0;
                    xtab[i + 1] = r * Math.Cos(a);
                    ytab[i + 1] = r * Math.Sin(a);
                }

                r = Math.Sqrt((6.0 + Math.Sqrt(6.0)) / 10.0);

                for (i = 0; i < 10; i++)
                {
                    a = 2.0 * pi * (double)(i) / 10.0;
                    xtab[i + 11] = r * Math.Cos(a);
                    ytab[i + 11] = r * Math.Sin(a);
                }
            }
            //
            //  There was apparently a misprint in the Lether paper.  The quantity
            //  which here reads "322" was printed there as "332".
            //
            else if (rule == 12)
            {
                a = 0.5;
                b = Math.Sqrt(3.0) / 2.0;
                c = Math.Sqrt((35.0 + 2.0 * Math.Sqrt(70.0)) / 252.0);
                d = Math.Sqrt((35.0 - 2.0 * Math.Sqrt(70.0)) / 252.0);
                e = Math.Sqrt((35.0 + 2.0 * Math.Sqrt(70.0)) / 84.0);
                f = Math.Sqrt((35.0 - 2.0 * Math.Sqrt(70.0)) / 84.0);
                g = Math.Sqrt((35.0 + 2.0 * Math.Sqrt(70.0)) / 63.0);
                h = Math.Sqrt((35.0 - 2.0 * Math.Sqrt(70.0)) / 63.0);

                w1 = 64.0 / 675.0;
                w2 = 16.0 / 225.0;
                w3 = 16.0 / 675.0;
                w4 = (322.0 - 13.0 * Math.Sqrt(70.0)) / 21600.0;
                w5 = (322.0 + 13.0 * Math.Sqrt(70.0)) / 21600.0;
                w6 = (322.0 - 13.0 * Math.Sqrt(70.0)) / 7200.0;
                w7 = (322.0 + 13.0 * Math.Sqrt(70.0)) / 7200.0;
                w8 = (322.0 - 13.0 * Math.Sqrt(70.0)) / 5400.0;
                w9 = (322.0 + 13.0 * Math.Sqrt(70.0)) / 5400.0;
                z = 0.0;

                xtab[0] = z;
                xtab[1] = a;
                xtab[2] = -a;
                xtab[3] = b;
                xtab[4] = -b;
                xtab[5] = b;
                xtab[6] = b;
                xtab[7] = -b;
                xtab[8] = -b;
                xtab[9] = b;
                xtab[10] = b;
                xtab[11] = -b;
                xtab[12] = -b;
                xtab[13] = a;
                xtab[14] = a;
                xtab[15] = -a;
                xtab[16] = -a;
                xtab[17] = a;
                xtab[18] = a;
                xtab[19] = -a;
                xtab[20] = -a;
                xtab[21] = z;
                xtab[22] = z;
                xtab[23] = z;
                xtab[24] = z;

                ytab[0] = z;
                ytab[1] = z;
                ytab[2] = z;
                ytab[3] = z;
                ytab[4] = z;
                ytab[5] = c;
                ytab[6] = -c;
                ytab[7] = c;
                ytab[8] = -c;
                ytab[9] = d;
                ytab[10] = -d;
                ytab[11] = d;
                ytab[12] = -d;
                ytab[13] = e;
                ytab[14] = -e;
                ytab[15] = e;
                ytab[16] = -e;
                ytab[17] = f;
                ytab[18] = -f;
                ytab[19] = f;
                ytab[20] = -f;
                ytab[21] = g;
                ytab[22] = -g;
                ytab[23] = h;
                ytab[24] = -h;

                weight[0] = w1;
                weight[1] = w2;
                weight[2] = w2;
                weight[3] = w3;
                weight[4] = w3;
                weight[5] = w4;
                weight[6] = w4;
                weight[7] = w4;
                weight[8] = w4;
                weight[9] = w5;
                weight[10] = w5;
                weight[11] = w5;
                weight[12] = w5;
                weight[13] = w6;
                weight[14] = w6;
                weight[15] = w6;
                weight[16] = w6;
                weight[17] = w7;
                weight[18] = w7;
                weight[19] = w7;
                weight[20] = w7;
                weight[21] = w8;
                weight[22] = w8;
                weight[23] = w9;
                weight[24] = w9;
            }
            else if (rule == 13)
            {
                nr = 4;
                ra = new double[nr];
                rw = new double[nr];

                LegendreQuadrature.legendre_set(nr, ref ra, ref rw);

                a = -1.0;
                b = +1.0;
                c = 0.0;
                d = +1.0;
                QuadratureRule.rule_adjust(a, b, c, d, nr, ref ra, ref rw);

                for (i = 0; i < nr; i++)
                {
                    ra[i] = Math.Sqrt(ra[i]);
                }

                i = 0;
                for (j = 0; j < 16; j++)
                {
                    c = Math.Cos(pi * (double)(j) / 8.0);
                    s = Math.Sin(pi * (double)(j) / 8.0);

                    for (k = 0; k < nr; k++)
                    {
                        xtab[i] = c * ra[k];
                        ytab[i] = s * ra[k];
                        weight[i] = rw[k] / 16.0;
                        i = i + 1;
                    }
                }
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("CIRCLE_XY_SET - Fatal error!");
                Console.WriteLine("  There is no rule of index " + rule + "");
            }

        }

        public static int circle_xy_size(int rule)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CIRCLE_XY_SIZE sizes an XY quadrature rule inside the unit circle in 2D.
            //
            //  Integration region:
            //
            //    X*X + Y*Y <= 1.0.
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
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //    Frank Lether,
            //    A Generalized Product Rule for the Circle,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 8, Number 2, June 1971, pages 249-253.
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, int RULE, the rule desired.
            //      1, 1 point 1-st degree;
            //      2, 4 point 3-rd degree, Stroud S2:3-1;
            //      3, 4 point 3-rd degree, Lether #1;
            //      4, 4 point 3-rd degree, Stroud S2:3-2;
            //      5, 5 point 3-rd degree;
            //      6, 7 point 5-th degree;
            //      7, 9 point 5-th degree;
            //      8, 9 point 5-th degree, Lether #2;
            //      9, 12 point 7-th degree;
            //     10, 16 point 7-th degree, Lether #3;
            //     11, 21 point 9-th degree, Stroud S2:9-3;
            //     12, 25 point 9-th degree, Lether #4 (after correcting error);
            //     13, 64 point 15-th degree Gauss product rule.
            //
            //    Output, int CIRCLE_XY_SIZE, the order of the desired rule.
            //
        {
            int order;

            if (rule == 1)
            {
                order = 1;
            }
            else if (rule == 2)
            {
                order = 4;
            }
            else if (rule == 3)
            {
                order = 4;
            }
            else if (rule == 4)
            {
                order = 4;
            }
            else if (rule == 5)
            {
                order = 5;
            }
            else if (rule == 6)
            {
                order = 7;
            }
            else if (rule == 7)
            {
                order = 9;
            }
            else if (rule == 8)
            {
                order = 9;
            }
            else if (rule == 9)
            {
                order = 12;
            }
            else if (rule == 10)
            {
                order = 16;
            }
            else if (rule == 11)
            {
                order = 21;
            }
            else if (rule == 12)
            {
                order = 25;
            }
            else if (rule == 13)
            {
                order = 64;
            }
            else
            {
                order = -1;
                Console.WriteLine("");
                Console.WriteLine("CIRCLE_XY_SIZE - Fatal error!");
                Console.WriteLine("  There is no rule of index " + rule + "");
                return (1);
            }

            return order;
        }

        public static double circle_xy_sum(Func<double, double, double> func, double[] center,
                double r, int order, double[] xtab, double[] ytab, double[] weight)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CIRCLE_XY_SUM applies an XY quadrature rule inside a circle in 2D.
            //
            //  Integration region:
            //
            //    (X-CENTER(1))^2 + (Y-CENTER(2))^2 <= R * R.
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
            //  Parameters:
            //
            //    Input, double FUNC ( double x, double y ), the name of the 
            //    user supplied function to be integrated.
            //
            //    Input, double CENTER[2], the coordinates of the center of 
            //    the circle.
            //
            //    Input, double R, the radius of the circle.
            //
            //    Input, int ORDER, the order of the rule.  The rule is
            //    assumed to be defined on the unit circle.
            //
            //    Input, double XTAB[ORDER], YTAB[ORDER], the XY
            //    coordinates of the abscissas of the quadrature rule for the unit circle.
            //
            //    Input, double WEIGHT[ORDER], the weights of the rule.
            //
            //    Output, double RESULT, the approximate integral of the function.
            //
        {
            int i;
            double quad;
            double result;
            double volume;
            double x;
            double y;

            quad = 0.0;
            for (i = 0; i < order; i++)
            {
                x = center[0] + r * xtab[i];
                y = center[1] + r * ytab[i];
                quad = quad + weight[i] * func(x, y);
            }

            volume = circle_area_2d(r);
            result = quad * volume;

            return result;
        }

    }
}