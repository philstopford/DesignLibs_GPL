using System;

namespace Burkardt.Stroud;

public static class Hexagon
{
    public static double hexagon_area_2d(double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HEXAGON_AREA_2D returns the area of a regular hexagon in 2D.
        //
        //  Discussion:
        //
        //    The formula for the area only requires the radius, and does
        //    not depend on the location of the center, or the orientation
        //    of the hexagon.
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
        //    Input, double R, the radius of the hexagon.
        //
        //    Output, double HEXAGON_AREA_2D, the area of the hexagon.
        //
    {
        double value = r * r * hexagon_unit_area_2d();

        return value;
    }

    public static double hexagon_sum(int setting, Func<int, double, double, double> func, double[] center,
            double r, int order, double[] xtab, double[] ytab, double[] weight)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HEXAGON_SUM applies a quadrature rule inside a hexagon in 2D.
        //
        //  Discussion:
        //
        //    The input quadrature rule is assumed to be defined for a unit hexagon.
        //
        //    The input quadrature rule may be defined by calling HEXAGON_UNIT_SET.
        //
        //  Integration region:
        //
        //    The definition is given in terms of THETA, the angle in degrees of the
        //    vector (X-CENTER(1),Y-CENTER(2)).  The following six conditions apply,
        //    respectively, between the bracketing values of THETA of 0, 60, 120, 
        //    180, 240, 300, and 360.
        //
        //      0 <= Y-CENTER(2) <= -SQRT(3) * (X-CENTER(1)) + R * SQRT(3)
        //      0 <= Y-CENTER(2) <=                     R * SQRT(3)/2
        //      0 <= Y-CENTER(2) <=  SQRT(3) * (X-CENTER(1)) + R * SQRT(3) 
        //      -SQRT(3) * (X-CENTER(1)) - R * SQRT(3)	<= Y-CENTER(2) <= 0
        //                        - R * SQRT(3)/2 <= Y-CENTER(2) <= 0
        //       SQRT(3) * (X-CENTER(1)) - R * SQRT(3)   <= Y-CENTER(2) <= 0
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
        //    Input, Func < double, double, double > func, the name of the 
        //    user supplied function of two variables which is to be integrated.
        //
        //    Input, double[] center, the center of the hexagon.
        //
        //    Input, double R, the radius of the hexagon.
        //
        //    Input, int ORDER, the order of the rule.
        //
        //    Input, double XTAB[ORDER], YTAB[ORDER], the abscissas.
        //
        //    Input, double WEIGHT[ORDER], the weights of the rule.
        //
        //    Output, double RESULT, the approximate integral of the function.
        //
    {
        int i;

        double quad = 0.0;
        for (i = 0; i < order; i++)
        {
            double x = center[0] + r * xtab[i];
            double y = center[1] + r * ytab[i];
            quad += weight[i] * func(setting, x, y);
        }

        double volume = hexagon_area_2d(r);
        double result = quad * volume;

        return result;
    }

    public static double hexagon_unit_area_2d()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HEXAGON_UNIT_AREA_2D returns the area of the unit regular hexagon in 2D.
        //
        //  Integration region:
        //
        //    The definition is given in terms of THETA, the angle in degrees of the
        //    vector (X,Y).  The following six conditions apply, respectively,
        //    between the bracketing values of THETA of 0, 60, 120, 180, 240,
        //    300, and 360.
        //
        //                              0 <= Y <= -SQRT(3) * X + SQRT(3)
        //                              0 <= Y <=                 SQRT(3)/2
        //                              0 <= Y <=  SQRT(3) * X + SQRT(3)
        //      - SQRT(3) * X - SQRT(3)   <= Y <= 0
        //                    - SQRT(3)/2 <= Y <= 0
        //        SQRT(3) * X - SQRT(3)   <= Y <= 0
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
        //    Output, double HEXAGON_UNIT_AREA_2D, the area of the hexagon.
        //
    {
        double value = 3.0 * Math.Sqrt(3.0) / 2.0;

        return value;
    }

    public static void hexagon_unit_set(int rule, int order, ref double[] xtab, ref double[] ytab,
            ref double[] weight)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HEXAGON_UNIT_SET sets a quadrature rule inside the unit hexagon in 2D.
        //
        //  Integration region:
        //
        //    The definition is given in terms of THETA, the angle in degrees of the
        //    vector (X,Y).  The following six conditions apply, respectively,
        //    between the bracketing values of THETA of 0, 60, 120, 180, 240,
        //    300, and 360.
        //
        //                              0 <= Y <= -SQRT(3) * X + SQRT(3)
        //                              0 <= Y <=                 SQRT(3)/2
        //                              0 <= Y <=  SQRT(3) * X + SQRT(3)
        //       -SQRT(3) * X - SQRT(3)   <= Y <= 0
        //                    - SQRT(3)/2 <= Y <= 0
        //        SQRT(3) * X - SQRT(3)   <= Y <= 0
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
        //      1, 1 point,  degree 1;
        //      2, 4 points, degree 3;
        //      3, 7 points, degree 3;
        //      4, 7 points, degree 5;
        //
        //    Input, int ORDER, the order of the desired rule.
        //
        //    Output, double XTAB[ORDER], YTAB[ORDER], the abscissas of the rule.
        //
        //    Output, double WEIGHT[ORDER], the weights of the rule.
        //
    {
        double a;
        double b;
        double c;
        double d;
        double e;
        double z;

        switch (rule)
        {
            case 1:
                xtab[0] = 0.0;
                ytab[0] = 0.0;
                weight[0] = 1.0;
                break;
            //
            //  Stroud rule H2:3-1.
            //
            case 2:
                a = Math.Sqrt(5.0 / 12.0);
                b = 1.0 / 4.0;
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
                break;
            //
            //  Stroud rule H2:3-2.
            //
            case 3:
                a = Math.Sqrt(3.0) / 2.0;
                b = 0.5;
                c = 1.0;
                d = 5.0 / 72.0;
                e = 42.0 / 72.0;
                z = 0.0;

                xtab[0] = z;
                xtab[1] = c;
                xtab[2] = -c;
                xtab[3] = b;
                xtab[4] = -b;
                xtab[5] = b;
                xtab[6] = -b;

                ytab[0] = z;
                ytab[1] = z;
                ytab[2] = z;
                ytab[3] = a;
                ytab[4] = a;
                ytab[5] = -a;
                ytab[6] = -a;

                weight[0] = e;
                weight[1] = d;
                weight[2] = d;
                weight[3] = d;
                weight[4] = d;
                weight[5] = d;
                weight[6] = d;
                break;
            //
            //  Stroud rule H2:5-1.
            //
            case 4:
                a = Math.Sqrt(14.0) / 5.0;
                b = Math.Sqrt(14.0) / 10.0;
                c = Math.Sqrt(42.0) / 10.0;
                d = 125.0 / 1008.0;
                e = 258.0 / 1008.0;
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
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("HEXAGON_UNIT_SET - Fatal error!");
                Console.WriteLine("  Illegal input value of RULE = " + rule + "");
                break;
        }

    }

    public static int hexagon_unit_size(int rule)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HEXAGON_UNIT_SIZE sizes a quadrature rule inside the unit hexagon in 2D.
        //
        //  Integration region:
        //
        //    The definition is given in terms of THETA, the angle in degrees of the
        //    vector (X,Y).  The following six conditions apply, respectively,
        //    between the bracketing values of THETA of 0, 60, 120, 180, 240,
        //    300, and 360.
        //
        //                              0 <= Y <= -SQRT(3) * X + SQRT(3)
        //                              0 <= Y <=                 SQRT(3)/2
        //                              0 <= Y <=  SQRT(3) * X + SQRT(3)
        //       -SQRT(3) * X - SQRT(3)   <= Y <= 0
        //                    - SQRT(3)/2 <= Y <= 0
        //        SQRT(3) * X - SQRT(3)   <= Y <= 0
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
        //      1, 1 point,  degree 1;
        //      2, 4 points, degree 3;
        //      3, 7 points, degree 3;
        //      4, 7 points, degree 5;
        //
        //    Output, int HEXAGON_UNIT_SIZE, the order of the desired rule.
        //    If RULE is not legal, then ORDER is returned as -1.
        //
    {
        int order;

        switch (rule)
        {
            case 1:
                order = 1;
                break;
            //
            //  Stroud rule H2:3-1.
            //
            case 2:
                order = 4;
                break;
            //
            //  Stroud rule H2:3-2.
            //
            case 3:
            //
            //  Stroud rule H2:5-1.
            //
            case 4:
                order = 7;
                break;
            default:
                order = -1;
                break;
        }

        return order;
    }

}