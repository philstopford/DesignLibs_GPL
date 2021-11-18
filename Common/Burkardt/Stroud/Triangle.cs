using System;
using Burkardt.Quadrature;

namespace Burkardt.Stroud;

public static class Triangle
{
    public static void triangle_rule_adjust(double[] xval, double[] yval, int order,
            double[] xtab, double[] ytab, double[] weight, ref double[] xtab2, ref double[] ytab2,
            ref double[] weight2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_RULE_ADJUST adjusts a unit quadrature rule to an arbitrary triangle.
        //
        //  Integration region:
        //
        //      (X,Y) = ALPHA * (X1,Y1) + BETA * (X2,Y2) + ( 1 - ALPHA - BETA ) * (X3,Y3)
        //    and
        //      0 <= ALPHA <= 1 - BETA
        //    and
        //      0 <= BETA <= 1 - ALPHA
        //
        //  Discussion:
        //
        //    This routine accepts as input abscissas and weights appropriate for
        //    quadrature in the unit triangle, and returns abscissas and weights
        //    appropriate for quadrature in a given triangle.
        //
        //    Once this routine has been called, an integral over the given triangle
        //    can be approximated as:
        //
        //      QUAD = sum ( 1 <= I <= ORDER ) WTAB2(I) * FUNC ( XTAB2(I), YTAB2(I) )
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
        //    Input, double XVAL[3], YVAL[3], the coordinates of the nodes.
        //
        //    Input, int ORDER, the order of the rule.
        //
        //    Input, double XTAB[ORDER], YTAB[ORDER], the abscissas for
        //    the unit triangle.
        //
        //    Input, double WEIGHT[ORDER], the weights for the unit triangle.
        //
        //    Output, double XTAB2[ORDER], YTAB2[ORDER], the adjusted
        //    abscissas.
        //
        //    Output, double WEIGHT2[ORDER], the adjusted weights.
        //
    {
        int i;

        double volume = triangle_volume(xval, yval);

        for (i = 0; i < order; i++)
        {
            xtab2[i] = xtab[i] * xval[0]
                       + ytab[i] * xval[1]
                       + (1.0 - xtab[i] - ytab[i]) * xval[2];

            ytab2[i] = xtab[i] * yval[0]
                       + ytab[i] * yval[1]
                       + (1.0 - xtab[i] - ytab[i]) * yval[2];

            weight2[i] = weight[i] * 2.0 * volume;
        }
    }

    public static double triangle_sub(int setting, Func<int, double, double, double> func, double[] xval,
            double[] yval, int nsub, int order, double[] xtab, double[] ytab,
            double[] weight)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_SUB carries out quadrature over subdivisions of a triangular region.
        //
        //  Integration region:
        //
        //      (X,Y) =       ALPHA          * ( XVAL[0], YVAL[0] )
        //            +               BETA   * ( XVAL[1], YVAL[1] )
        //            + ( 1 - ALPHA - BETA ) * ( XVAL[2], YVAL[2] )
        //    and
        //      0 <= ALPHA <= 1 - BETA
        //    and
        //      0 <= BETA <= 1 - ALPHA
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
        //  Parameters:
        //
        //    Input, Func < double, double, double > func, the name of the user supplied
        //    function to be integrated.
        //
        //    Input, double XVAL[3], YVAL[3], the coordinates of the triangle vertices.
        //
        //    Input, int NSUB, the number of subdivisions of each side of the
        //    input triangle to be made.  NSUB = 1 means no subdivisions are made.
        //    NSUB = 3 means that each side of the triangle is subdivided into
        //    three portions, and that the original triangle is subdivided into
        //    NSUB * NSUB triangles.  NSUB must be at least 1.
        //
        //    Input, int ORDER, the order of the rule.
        //
        //    Input, double XTAB[ORDER], YTAB[ORDER], the abscissas.
        //
        //    Input, double WEIGHT[ORDER], the weights of the rule.
        //
        //    Output, double TRIANGLE_SUB, the approximate integral of the function.
        //
    {
        int i;
        //
        //  Initialize RESULT, the approximate integral.
        //
        double result = 0.0;
        switch (nsub)
        {
            //
            //  NSUB must be positive.
            //
            case <= 0:
                return result;
        }

        //
        //  Initialize QUAD, the quadrature sum.
        //
        double quad = 0.0;
        //
        //  The sub-triangles can be grouped into NSUB strips.
        //
        for (i = 1; i <= nsub; i++)
        {
            double temp1 = 0.0;
            double temp2 = i / (double)nsub;

            double x2 = xval[1] + temp1 * (xval[2] - xval[1])
                                + temp2 * (xval[0] - xval[1]);

            double y2 = yval[1] + temp1 * (yval[2] - yval[1])
                                + temp2 * (yval[0] - yval[1]);

            temp1 = 0.0;
            temp2 = (i - 1) / (double)nsub;

            double x3 = xval[1] + temp1 * (xval[2] - xval[1])
                                + temp2 * (xval[0] - xval[1]);

            double y3 = yval[1] + temp1 * (yval[2] - yval[1])
                                + temp2 * (yval[0] - yval[1]);
            //
            //  There are 2*I-1 triangles in strip number I.
            //  The next triangle in the strip shares two nodes with the previous one.
            //  Compute its corners, (X1,Y1), (X2,Y2), (X3,Y3).
            //
            int j;
            for (j = 1; j <= 2 * i - 1; j++)
            {
                double x1 = x2;
                double y1 = y2;
                x2 = x3;
                y2 = y3;
                temp1 = (double)(j + 1) / 2 / nsub;
                temp2 = (double)(i - 1 - j / 2) / nsub;

                x3 = xval[1] + temp1 * (xval[2] - xval[1])
                             + temp2 * (xval[0] - xval[1]);

                y3 = yval[1] + temp1 * (yval[2] - yval[1])
                             + temp2 * (yval[0] - yval[1]);
                //
                //  Now integrate over the triangle, mapping the points ( XTAB(K), YTAB(K) )
                //  into the triangle.
                //
                int k;
                for (k = 0; k < order; k++)
                {
                    double x = x2 + xtab[k] * (x3 - x2) + ytab[k] * (x1 - x2);
                    double y = y2 + xtab[k] * (y3 - y2) + ytab[k] * (y1 - y2);
                    quad += weight[k] * func(setting, x, y);
                }
            }
        }

        double volume = triangle_volume(xval, yval) / (nsub * nsub);
        result = quad * volume;

        return result;
    }

    public static double triangle_sum(int setting, Func<int, double, double, double> func,
            double[] xval, double[] yval, int order, double[] xtab, double[] ytab,
            double[] weight)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_SUM carries out a unit quadrature rule in an arbitrary triangle.
        //
        //  Integration region:
        //
        //      (X,Y) =       ALPHA          * (X1,Y1) 
        //            +               BETA   * (X2,Y2) 
        //            + ( 1 - ALPHA - BETA ) * (X3,Y3)
        //    and
        //      0 <= ALPHA <= 1 - BETA
        //    and
        //      0 <= BETA <= 1 - ALPHA
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
        //  Parameters:
        //
        //    Input, Func< double, double, double, double > func, the name of the 
        //    user supplied function to be integrated.
        //
        //    Input, double XVAL[3], YVAL[3], the coordinates of the nodes.
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
            double x = xtab[i] * xval[0]
                       + ytab[i] * xval[1]
                       + (1.0 - xtab[i] - ytab[i]) * xval[2];

            double y = xtab[i] * yval[0]
                       + ytab[i] * yval[1]
                       + (1.0 - xtab[i] - ytab[i]) * yval[2];

            quad += weight[i] * func(setting, x, y);
        }

        double volume = triangle_volume(xval, yval);
        double result = quad * volume;

        return result;
    }

    public static double triangle_sum_adjusted(Func<double, double, double> func,
            int order, double[] xtab, double[] ytab, double[] weight)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_SUM_ADJUSTED carries out an adjusted quadrature rule in a triangle.
        //
        //  Integration region:
        //
        //      (X,Y) =       ALPHA          * (X1,Y1) 
        //                          + BETA   * (X2,Y2) 
        //            + ( 1 - ALPHA - BETA ) * (X3,Y3)
        //    and
        //      0 <= ALPHA <= 1 - BETA
        //    and
        //      0 <= BETA <= 1 - ALPHA
        //
        //  Discussion:
        //
        //    It is assumed that a quadrature rule approprate for the unit triangle
        //    was generated, and then adjusted to a particular triangle by calling
        //    TRIANGLE_RULE_ADJUST.
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
        //    Input, Func < double, double, double > func, the name of the 
        //    user supplied function to be integrated.
        //
        //    Input, int ORDER, the order of the rule.
        //
        //    Input, double XTAB[ORDER], YTAB[ORDER], the abscissas.
        //
        //    Input, double WEIGHT[ORDER], the weights of the rule.
        //
        //    Output, double TRIANGLE_SUM_ADJUSTED, the approximate integral 
        //    of the function.
        //
    {
        int i;

        double result = 0.0;

        for (i = 0; i < order; i++)
        {
            result += weight[i] * func(xtab[i], ytab[i]);
        }

        return result;
    }

    public static void triangle_unit_product_set(int rule, int order, ref double[] xtab,
            ref double[] ytab, ref double[] weight)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_PRODUCT_SET sets a product rule on the unit triangle.
        //
        //  Discussion:
        //
        //    For a given order of accuracy, a product rule on a triangle usually
        //    uses more points than necessary.  That is, there is usually a rule
        //    of the same order that uses fewer points.
        //
        //    However, one advantage of product rules is that a rule of any
        //    desired order can be generated automatically.
        //   
        //    The integration region is:
        //
        //      0 <= X,
        //    and
        //      0 <= Y, 
        //    and
        //      X + Y <= 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 September 2010
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
        //    Input, int RULE, the order of the 1D rule.
        //
        //    Input, int ORDER, the order of the rule.
        //
        //    Output, double XTAB[ORDER], YTAB[ORDER], the abscissas.
        //
        //    Output, double WEIGHT[ORDER], the weights of the rule.
        //
    {
        int j;

        double[] weight0 = new double[rule];
        double[] weight1 = new double[rule];
        double[] xtab0 = new double[rule];
        double[] xtab1 = new double[rule];

        const double a = -1.0;
        const double b = +1.0;
        const double c = 0.0;
        const double d = +1.0;

        int order0 = rule;
        LegendreQuadrature.legendre_set(order0, ref xtab0, ref weight0);
        QuadratureRule.rule_adjust(a, b, c, d, order0, ref xtab0, ref weight0);

        int order1 = rule;
        LegendreQuadrature.legendre_set_x1(order1, ref xtab1, ref weight1);
        QuadratureRule.rule_adjust(a, b, c, d, order1, ref xtab1, ref weight1);

        int k = 0;
        for (j = 0; j < order1; j++)
        {
            int i;
            for (i = 0; i < order0; i++)
            {
                xtab[k] = 1.0 - xtab1[j];
                ytab[k] = xtab0[i] * xtab1[j];
                weight[k] = weight0[i] * weight1[j];
                k += 1;
            }
        }
    }

    public static int triangle_unit_product_size(int rule)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_PRODUCT_SIZE sizes a product rule on the unit triangle.
        //
        //  Discussion:
        //
        //    The integration region is:
        //
        //      0 <= X,
        //    and
        //      0 <= Y, 
        //    and
        //      X + Y <= 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 April 2008
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
        //    Input, int RULE, the order of the 1D rule.
        //
        //    Input, int TRIANGLE_UNIT_PRODUCT_SIZE, the order of the rule.
        //
    {
        int order = rule * rule;

        return order;
    }

    public static void triangle_unit_set(int rule, int order, ref double[] xtab, ref double[] ytab,
            ref double[] weight)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_SET sets a quadrature rule in the unit triangle.
        //
        //  Discussion:
        //
        //    The user is responsible for determining the value of ORDER,
        //    and appropriately dimensioning the arrays XTAB, YTAB and
        //    WEIGHT so that they can accommodate the data.
        //
        //    The value of ORDER for each rule can be found by invoking
        //    the function TRIANGLE_RULE_SIZE.
        //
        //  Integration region:
        //
        //      0 <= X,
        //    and
        //      0 <= Y, 
        //    and
        //      X + Y <= 1.
        //
        //  Graph:
        //
        //      ^
        //    1 | *
        //      | |.
        //    Y | | .
        //      | |  .
        //    0 | *---*
        //      +------->
        //        0 X 1
        //
        //   The rules are accessed by an index number, RULE.  The indices,
        //   and the descriptions of the corresponding rules, are:
        //
        //     1, ORDER =  1, precision 1, Zienkiewicz #1.
        //     2, ORDER =  2, precision 1, (the "vertex rule").
        //     3, ORDER =  3, precision 2, Strang and Fix formula #1.
        //     4, ORDER =  3, precision 2, Strang and Fix formula #2,
        //                                 Zienkiewicz #2.
        //     5, ORDER =  4, precision 3, Strang and Fix formula #3,
        //                                 Zienkiewicz #3.
        //     6, ORDER =  6, precision 3, Strang and Fix formula #4.
        //     7, ORDER =  6, precision 3, Stroud formula T2:3-1.
        //     8, ORDER =  6, precision 4, Strang and Fix formula #5.
        //     9, ORDER =  7, precision 4, Strang and Fix formula #6.
        //    10, ORDER =  7, precision 5, Strang and Fix formula #7,
        //                                 Stroud formula T2:5-1, 
        //                                 Zienkiewicz #4, 
        //                                 Schwarz Table 2.2.
        //    11, ORDER =  9, precision 6, Strang and Fix formula #8.
        //    12, ORDER = 12, precision 6, Strang and Fix formula #9.
        //    13, ORDER = 13, precision 7, Strang and Fix formula #10.
        //        Note that there is a typographical error in Strang and Fix
        //        which lists the value of the XSI(3) component of the
        //        last generator point as 0.4869... when it should be 0.04869...
        //    14, ORDER =  7, precision 3.
        //    15, ORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
        //    16, ORDER = 64, precision 15, triangular product Gauss rule.
        //    17, ORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
        //    18, ORDER = 19, precision 9, from TRIEX, ACM TOMS #612.
        //    19, ORDER = 28, precision 11, from TRIEX, ACM TOMS #612.
        //    20, ORDER = 37, precision 13, from ACM TOMS #706.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Jarle Berntsen, Terje Espelid,
        //    Algorithm 706,
        //    DCUTRI: an algorithm for adaptive cubature over a collection of triangles, 
        //    ACM Transactions on Mathematical Software,
        //    Volume 18, Number 3, September 1992, pages 329-342.
        //
        //    Elise deDoncker, Ian Robinson,
        //    Algorithm 612:
        //    Integration over a Triangle Using Nonlinear Extrapolation,
        //    ACM Transactions on Mathematical Software,
        //    Volume 10, Number 1, March 1984, pages 17-22.
        //
        //    Dirk Laurie,
        //    Algorithm 584,
        //    CUBTRI, Automatic Cubature Over a Triangle,
        //    ACM Transactions on Mathematical Software,
        //    Volume 8, Number 2, 1982, pages 210-218.
        //
        //    James Lyness, Dennis Jespersen,
        //    Moderate Degree Symmetric Quadrature Rules for the Triangle,
        //    Journal of the Institute of Mathematics and its Applications,
        //    Volume 15, Number 1, February 1975, pages 19-32.
        //
        //    Hans Rudolf Schwarz,
        //    Finite Element Methods,
        //    Academic Press, 1988,
        //    ISBN: 0126330107,
        //    LC: TA347.F5.S3313.
        //
        //    Gilbert Strang, George Fix,
        //    An Analysis of the Finite Element Method,
        //    Cambridge, 1973,
        //    ISBN: 096140888X,
        //    LC: TA335.S77.
        //
        //    Arthur Stroud,
        //    Approximate Calculation of Multiple Integrals,
        //    Prentice Hall, 1971,
        //    ISBN: 0130438936,
        //    LC: QA311.S85.
        //
        //    Olgierd Zienkiewicz,
        //    The Finite Element Method,
        //    Sixth Edition,
        //    Butterworth-Heinemann, 2005,
        //    ISBN: 0750663200,
        //    LC: TA640.2.Z54
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //
        //    Input, int ORDER, the order of the rule.
        //
        //    Output, double XTAB[ORDER], YTAB[ORDER], the abscissas.
        //
        //    Output, double WEIGHT[ORDER], the weights of the rule.
        //
    {
        double a;
        double b;
        double c;
        double d;
        double e;
        double f;
        double g;
        int i;
        int j;
        int k;
        int order2;
        double p;
        double q;
        double r;
        double s;
        double t;
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
        double[] weight1 = new double[8];
        double[] weight2 = new double[8];
        double[] xtab1 = new double[8];
        double[] xtab2 = new double[8];

        switch (rule)
        {
            //
            //  1 point, precision 1.
            //
            case 1:
                xtab[0] = 0.33333333333333333333;

                ytab[0] = 0.33333333333333333333;

                weight[0] = 1.00000000000000000000;
                break;
            //
            //  3 points, precision 1, the "vertex rule".
            //
            case 2:
                xtab[0] = 1.00000000000000000000;
                xtab[1] = 0.00000000000000000000;
                xtab[2] = 0.00000000000000000000;

                ytab[0] = 0.00000000000000000000;
                ytab[1] = 1.00000000000000000000;
                ytab[2] = 0.00000000000000000000;

                weight[0] = 0.33333333333333333333;
                weight[1] = 0.33333333333333333333;
                weight[2] = 0.33333333333333333333;
                break;
            //
            //  3 points, precision 2, Strang and Fix formula #1.
            //
            case 3:
                xtab[0] = 0.66666666666666666667;
                xtab[1] = 0.16666666666666666667;
                xtab[2] = 0.16666666666666666667;

                ytab[0] = 0.16666666666666666667;
                ytab[1] = 0.66666666666666666667;
                ytab[2] = 0.16666666666666666667;

                weight[0] = 0.33333333333333333333;
                weight[1] = 0.33333333333333333333;
                weight[2] = 0.33333333333333333333;
                break;
            //
            //  3 points, precision 2, Strang and Fix formula #2.
            //
            case 4:
                xtab[0] = 0.50000000000000000000;
                xtab[1] = 0.50000000000000000000;
                xtab[2] = 0.00000000000000000000;

                ytab[0] = 0.00000000000000000000;
                ytab[1] = 0.50000000000000000000;
                ytab[2] = 0.50000000000000000000;

                weight[0] = 0.33333333333333333333;
                weight[1] = 0.33333333333333333333;
                weight[2] = 0.33333333333333333333;
                break;
            //
            //  4 points, precision 3, Strang and Fix formula #3.
            //
            case 5:
                a = 6.0 / 30.0;
                b = 10.0 / 30.0;
                c = 18.0 / 30.0;

                d = 25.0 / 48.0;
                e = -27.0 / 48.0;

                xtab[0] = b;
                xtab[1] = c;
                xtab[2] = a;
                xtab[3] = a;

                ytab[0] = b;
                ytab[1] = a;
                ytab[2] = c;
                ytab[3] = a;

                weight[0] = e;
                weight[1] = d;
                weight[2] = d;
                weight[3] = d;
                break;
            //
            //  6 points, precision 3, Strang and Fix formula #4.
            //
            case 6:
                a = 0.659027622374092;
                b = 0.231933368553031;
                c = 0.109039009072877;

                xtab[0] = a;
                xtab[1] = a;
                xtab[2] = b;
                xtab[3] = b;
                xtab[4] = c;
                xtab[5] = c;

                ytab[0] = b;
                ytab[1] = c;
                ytab[2] = a;
                ytab[3] = c;
                ytab[4] = a;
                ytab[5] = b;

                weight[0] = 0.16666666666666666667;
                weight[1] = 0.16666666666666666667;
                weight[2] = 0.16666666666666666667;
                weight[3] = 0.16666666666666666667;
                weight[4] = 0.16666666666666666667;
                weight[5] = 0.16666666666666666667;
                break;
            //
            //  6 points, precision 3, Stroud T2:3-1.
            //
            case 7:
                a = 0.0;
                b = 0.5;
                c = 2.0 / 3.0;
                d = 1.0 / 6.0;
                v = 1.0 / 30.0;
                w = 3.0 / 10.0;

                xtab[0] = a;
                xtab[1] = b;
                xtab[2] = c;
                xtab[3] = c;
                xtab[4] = d;
                xtab[5] = d;

                ytab[0] = b;
                ytab[1] = a;
                ytab[2] = b;
                ytab[3] = d;
                ytab[4] = c;
                ytab[5] = d;

                weight[0] = v;
                weight[1] = v;
                weight[2] = v;
                weight[3] = w;
                weight[4] = w;
                weight[5] = w;
                break;
            //
            //  6 points, precision 4, Strang and Fix, formula #5.
            //
            case 8:
                a = 0.816847572980459;
                b = 0.091576213509771;
                c = 0.108103018168070;
                d = 0.445948490915965;
                v = 0.109951743655322;
                w = 0.223381589678011;

                xtab[0] = a;
                xtab[1] = b;
                xtab[2] = b;
                xtab[3] = c;
                xtab[4] = d;
                xtab[5] = d;

                ytab[0] = b;
                ytab[1] = a;
                ytab[2] = b;
                ytab[3] = d;
                ytab[4] = c;
                ytab[5] = d;

                weight[0] = v;
                weight[1] = v;
                weight[2] = v;
                weight[3] = w;
                weight[4] = w;
                weight[5] = w;
                break;
            //
            //  7 points, precision 4, Strang and Fix formula #6.
            //
            case 9:
                a = 1.0 / 3.0;
                c = 0.736712498968435;
                d = 0.237932366472434;
                e = 0.025355134551932;
                v = 0.375000000000000;
                w = 0.104166666666667;

                xtab[0] = a;
                xtab[1] = c;
                xtab[2] = c;
                xtab[3] = d;
                xtab[4] = d;
                xtab[5] = e;
                xtab[6] = e;

                ytab[0] = a;
                ytab[1] = d;
                ytab[2] = e;
                ytab[3] = c;
                ytab[4] = e;
                ytab[5] = c;
                ytab[6] = d;

                weight[0] = v;
                weight[1] = w;
                weight[2] = w;
                weight[3] = w;
                weight[4] = w;
                weight[5] = w;
                weight[6] = w;
                break;
            //
            //  7 points, precision 5, Strang and Fix formula #7, Stroud T2:5-1
            //
            case 10:
                a = 1.0 / 3.0;
                b = (9.0 + 2.0 * Math.Sqrt(15.0)) / 21.0;
                c = (6.0 - Math.Sqrt(15.0)) / 21.0;
                d = (9.0 - 2.0 * Math.Sqrt(15.0)) / 21.0;
                e = (6.0 + Math.Sqrt(15.0)) / 21.0;
                u = 0.225;
                v = (155.0 - Math.Sqrt(15.0)) / 1200.0;
                w = (155.0 + Math.Sqrt(15.0)) / 1200.0;

                xtab[0] = a;
                xtab[1] = b;
                xtab[2] = c;
                xtab[3] = c;
                xtab[4] = d;
                xtab[5] = e;
                xtab[6] = e;

                ytab[0] = a;
                ytab[1] = c;
                ytab[2] = b;
                ytab[3] = c;
                ytab[4] = e;
                ytab[5] = d;
                ytab[6] = e;

                weight[0] = u;
                weight[1] = v;
                weight[2] = v;
                weight[3] = v;
                weight[4] = w;
                weight[5] = w;
                weight[6] = w;
                break;
            //
            //  9 points, precision 6, Strang and Fix formula #8.
            //
            case 11:
                a = 0.124949503233232;
                b = 0.437525248383384;
                c = 0.797112651860071;
                d = 0.165409927389841;
                e = 0.037477420750088;

                u = 0.205950504760887;
                v = 0.063691414286223;

                xtab[0] = a;
                xtab[1] = b;
                xtab[2] = b;
                xtab[3] = c;
                xtab[4] = c;
                xtab[5] = d;
                xtab[6] = d;
                xtab[7] = e;
                xtab[8] = e;

                ytab[0] = b;
                ytab[1] = a;
                ytab[2] = b;
                ytab[3] = d;
                ytab[4] = e;
                ytab[5] = c;
                ytab[6] = e;
                ytab[7] = c;
                ytab[8] = d;

                weight[0] = u;
                weight[1] = u;
                weight[2] = u;
                weight[3] = v;
                weight[4] = v;
                weight[5] = v;
                weight[6] = v;
                weight[7] = v;
                weight[8] = v;
                break;
            //
            //  12 points, precision 6, Strang and Fix, formula #9.
            //
            case 12:
                a = 0.873821971016996;
                b = 0.063089014491502;
                c = 0.501426509658179;
                d = 0.249286745170910;
                e = 0.636502499121399;
                f = 0.310352451033785;
                g = 0.053145049844816;

                u = 0.050844906370207;
                v = 0.116786275726379;
                w = 0.082851075618374;

                xtab[0] = a;
                xtab[1] = b;
                xtab[2] = b;
                xtab[3] = c;
                xtab[4] = d;
                xtab[5] = d;
                xtab[6] = e;
                xtab[7] = e;
                xtab[8] = f;
                xtab[9] = f;
                xtab[10] = g;
                xtab[11] = g;

                ytab[0] = b;
                ytab[1] = a;
                ytab[2] = b;
                ytab[3] = d;
                ytab[4] = c;
                ytab[5] = d;
                ytab[6] = f;
                ytab[7] = g;
                ytab[8] = e;
                ytab[9] = g;
                ytab[10] = e;
                ytab[11] = f;

                weight[0] = u;
                weight[1] = u;
                weight[2] = u;
                weight[3] = v;
                weight[4] = v;
                weight[5] = v;
                weight[6] = w;
                weight[7] = w;
                weight[8] = w;
                weight[9] = w;
                weight[10] = w;
                weight[11] = w;
                break;
            //
            //  13 points, precision 7, Strang and Fix, formula #10.
            //
            //  Note that there is a typographical error in Strang and Fix
            //  which lists the value of the XSI[2] component of the
            //  last generator point as 0.4869... when it should be 0.04869...
            //
            case 13:
                double h = 1.0 / 3.0;
                a = 0.479308067841923;
                b = 0.260345966079038;
                c = 0.869739794195568;
                d = 0.065130102902216;
                e = 0.638444188569809;
                f = 0.312865496004875;
                g = 0.048690315425316;

                w = -0.149570044467670;
                t = 0.175615257433204;
                u = 0.053347235608839;
                v = 0.077113760890257;

                xtab[0] = h;
                xtab[1] = a;
                xtab[2] = b;
                xtab[3] = b;
                xtab[4] = c;
                xtab[5] = d;
                xtab[6] = d;
                xtab[7] = e;
                xtab[8] = e;
                xtab[9] = f;
                xtab[10] = f;
                xtab[11] = g;
                xtab[12] = g;

                ytab[0] = h;
                ytab[1] = b;
                ytab[2] = a;
                ytab[3] = b;
                ytab[4] = d;
                ytab[5] = c;
                ytab[6] = d;
                ytab[7] = f;
                ytab[8] = g;
                ytab[9] = e;
                ytab[10] = g;
                ytab[11] = e;
                ytab[12] = f;

                weight[0] = w;
                weight[1] = t;
                weight[2] = t;
                weight[3] = t;
                weight[4] = u;
                weight[5] = u;
                weight[6] = u;
                weight[7] = v;
                weight[8] = v;
                weight[9] = v;
                weight[10] = v;
                weight[11] = v;
                weight[12] = v;
                break;
            //
            //  7 points, precision 3.
            //
            case 14:
                a = 1.0 / 3.0;
                b = 1.0;
                c = 0.5;
                double z = 0.0;

                u = 27.0 / 60.0;
                v = 3.0 / 60.0;
                w = 8.0 / 60.0;

                xtab[0] = a;
                xtab[1] = b;
                xtab[2] = z;
                xtab[3] = z;
                xtab[4] = z;
                xtab[5] = c;
                xtab[6] = c;

                ytab[0] = a;
                ytab[1] = z;
                ytab[2] = b;
                ytab[3] = z;
                ytab[4] = c;
                ytab[5] = z;
                ytab[6] = c;

                weight[0] = u;
                weight[1] = v;
                weight[2] = v;
                weight[3] = v;
                weight[4] = w;
                weight[5] = w;
                weight[6] = w;
                break;
            //
            //  16 points, precision 5, Stroud T2:7-1.
            //
            case 15:
            {
                //
                //  Legendre rule of order 4.
                //
                order2 = 4;

                xtab[0] = -0.861136311594052575223946488893;
                xtab[1] = -0.339981043584856264802665759103;
                xtab[2] = 0.339981043584856264802665759103;
                xtab[3] = 0.861136311594052575223946488893;

                weight1[0] = 0.347854845137453857373063949222;
                weight1[1] = 0.652145154862546142626936050778;
                weight1[2] = 0.652145154862546142626936050778;
                weight1[3] = 0.347854845137453857373063949222;

                for (i = 0; i < order2; i++)
                {
                    xtab1[i] = 0.5 * (xtab1[i] + 1.0);
                }

                weight2[0] = 0.1355069134;
                weight2[1] = 0.2034645680;
                weight2[2] = 0.1298475476;
                weight2[3] = 0.0311809709;

                xtab2[0] = 0.0571041961;
                xtab2[1] = 0.2768430136;
                xtab2[2] = 0.5835904324;
                xtab2[3] = 0.8602401357;

                k = 0;
                for (i = 0; i < order2; i++)
                {
                    for (j = 0; j < order2; j++)
                    {
                        xtab[k] = xtab2[j];
                        ytab[k] = xtab1[i] * (1.0 - xtab2[j]);
                        weight[k] = weight1[i] * weight2[j];
                        k += 1;
                    }
                }

                break;
            }
            //
            //  64 points, precision 15.
            //
            case 16:
            {
                //
                //  Legendre rule of order 8.
                //
                order2 = 8;

                xtab1[0] = -0.960289856497536231683560868569;
                xtab1[1] = -0.796666477413626739591553936476;
                xtab1[2] = -0.525532409916328985817739049189;
                xtab1[3] = -0.183434642495649804939476142360;
                xtab1[4] = 0.183434642495649804939476142360;
                xtab1[5] = 0.525532409916328985817739049189;
                xtab1[6] = 0.796666477413626739591553936476;
                xtab1[7] = 0.960289856497536231683560868569;

                weight1[0] = 0.101228536290376259152531354310;
                weight1[1] = 0.222381034453374470544355994426;
                weight1[2] = 0.313706645877887287337962201987;
                weight1[3] = 0.362683783378361982965150449277;
                weight1[4] = 0.362683783378361982965150449277;
                weight1[5] = 0.313706645877887287337962201987;
                weight1[6] = 0.222381034453374470544355994426;
                weight1[7] = 0.101228536290376259152531354310;

                weight2[0] = 0.00329519144;
                weight2[1] = 0.01784290266;
                weight2[2] = 0.04543931950;
                weight2[3] = 0.07919959949;
                weight2[4] = 0.10604735944;
                weight2[5] = 0.11250579947;
                weight2[6] = 0.09111902364;
                weight2[7] = 0.04455080436;

                xtab2[0] = 0.04463395529;
                xtab2[1] = 0.14436625704;
                xtab2[2] = 0.28682475714;
                xtab2[3] = 0.45481331520;
                xtab2[4] = 0.62806783542;
                xtab2[5] = 0.78569152060;
                xtab2[6] = 0.90867639210;
                xtab2[7] = 0.98222008485;

                k = 0;
                for (i = 0; i < order2; i++)
                {
                    for (j = 0; j < order2; j++)
                    {
                        xtab[k] = 1.0 - xtab2[j];
                        ytab[k] = 0.5 * (1.0 + xtab1[i]) * xtab2[j];
                        weight[k] = weight1[i] * weight2[j];
                        k += 1;
                    }
                }

                break;
            }
            //
            //  19 points, precision 8, from CUBTRI.
            //
            case 17:
                a = 1.0 / 3.0;
                b = (9.0 + 2.0 * Math.Sqrt(15.0)) / 21.0;
                c = (6.0 - Math.Sqrt(15.0)) / 21.0;
                d = (9.0 - 2.0 * Math.Sqrt(15.0)) / 21.0;
                e = (6.0 + Math.Sqrt(15.0)) / 21.0;
                f = (40.0 - 10.0 * Math.Sqrt(15.0)
                     + 10.0 * Math.Sqrt(7.0) + 2.0 * Math.Sqrt(105.0)) / 90.0;
                g = (25.0 + 5.0 * Math.Sqrt(15.0)
                     - 5.0 * Math.Sqrt(7.0) - Math.Sqrt(105.0)) / 90.0;
                p = (40.0 + 10.0 * Math.Sqrt(15.0)
                          + 10.0 * Math.Sqrt(7.0) - 2.0 * Math.Sqrt(105.0)) / 90.0;
                q = (25.0 - 5.0 * Math.Sqrt(15.0)
                          - 5.0 * Math.Sqrt(7.0) + Math.Sqrt(105.0)) / 90.0;
                r = (40.0 + 10.0 * Math.Sqrt(7.0)) / 90.0;
                s = (25.0 + 5.0 * Math.Sqrt(15.0) - 5.0 * Math.Sqrt(7.0)
                                                  - Math.Sqrt(105.0)) / 90.0;
                t = (25.0 - 5.0 * Math.Sqrt(15.0) - 5.0 * Math.Sqrt(7.0)
                     + Math.Sqrt(105.0)) / 90.0;

                w1 = (7137.0 - 1800.0 * Math.Sqrt(7.0)) / 62720.0;
                w2 = -9301697.0 / 4695040.0 - 13517313.0 * Math.Sqrt(15.0)
                    / 23475200.0 + 764885.0 * Math.Sqrt(7.0) / 939008.0
                                 + 198763.0 * Math.Sqrt(105.0) / 939008.0;
                w2 /= 3.0;
                w3 = -9301697.0 / 4695040.0 + 13517313.0 * Math.Sqrt(15.0)
                                            / 23475200.0
                                            + 764885.0 * Math.Sqrt(7.0) / 939008.0
                     - 198763.0 * Math.Sqrt(105.0) / 939008.0;
                w3 /= 3.0;
                w4 = (102791225.0 - 23876225.0 * Math.Sqrt(15.0)
                                  - 34500875.0 * Math.Sqrt(7.0)
                      + 9914825.0 * Math.Sqrt(105.0)) / 59157504.0;
                w4 /= 3.0;
                w5 = (102791225.0 + 23876225.0 * Math.Sqrt(15.0)
                      - 34500875.0 * Math.Sqrt(7.0)
                      - 9914825 * Math.Sqrt(105.0)) / 59157504.0;
                w5 /= 3.0;
                w6 = (11075.0 - 3500.0 * Math.Sqrt(7.0)) / 8064.0;
                w6 /= 6.0;

                xtab[0] = a;
                xtab[1] = b;
                xtab[2] = c;
                xtab[3] = c;
                xtab[4] = d;
                xtab[5] = e;
                xtab[6] = e;
                xtab[7] = f;
                xtab[8] = g;
                xtab[9] = g;
                xtab[10] = p;
                xtab[11] = q;
                xtab[12] = q;
                xtab[13] = r;
                xtab[14] = r;
                xtab[15] = s;
                xtab[16] = s;
                xtab[17] = t;
                xtab[18] = t;

                ytab[0] = a;
                ytab[1] = c;
                ytab[2] = b;
                ytab[3] = c;
                ytab[4] = e;
                ytab[5] = d;
                ytab[6] = e;
                ytab[7] = g;
                ytab[8] = f;
                ytab[9] = g;
                ytab[10] = q;
                ytab[11] = p;
                ytab[12] = q;
                ytab[13] = s;
                ytab[14] = t;
                ytab[15] = r;
                ytab[16] = t;
                ytab[17] = r;
                ytab[18] = s;

                weight[0] = w1;
                weight[1] = w2;
                weight[2] = w2;
                weight[3] = w2;
                weight[4] = w3;
                weight[5] = w3;
                weight[6] = w3;
                weight[7] = w4;
                weight[8] = w4;
                weight[9] = w4;
                weight[10] = w5;
                weight[11] = w5;
                weight[12] = w5;
                weight[13] = w6;
                weight[14] = w6;
                weight[15] = w6;
                weight[16] = w6;
                weight[17] = w6;
                weight[18] = w6;
                break;
            //
            //  19 points, precision 9.
            //  Lyness and Jesperson.
            //
            case 18:
                a = 1.0 / 3.0;
                b = 0.02063496160252593;
                c = 0.4896825191987370;
                d = 0.1258208170141290;
                e = 0.4370895914929355;
                f = 0.6235929287619356;
                g = 0.1882035356190322;
                r = 0.9105409732110941;
                s = 0.04472951339445297;
                t = 0.7411985987844980;
                u = 0.03683841205473626;
                v = 0.22196298916076574;

                w1 = 0.09713579628279610;
                w2 = 0.03133470022713983;
                w3 = 0.07782754100477543;
                w4 = 0.07964773892720910;
                w5 = 0.02557767565869810;
                w6 = 0.04328353937728940;

                xtab[0] = a;
                xtab[1] = b;
                xtab[2] = c;
                xtab[3] = c;
                xtab[4] = d;
                xtab[5] = e;
                xtab[6] = e;
                xtab[7] = f;
                xtab[8] = g;
                xtab[9] = g;
                xtab[10] = r;
                xtab[11] = s;
                xtab[12] = s;
                xtab[13] = t;
                xtab[14] = t;
                xtab[15] = u;
                xtab[16] = u;
                xtab[17] = v;
                xtab[18] = v;

                ytab[0] = a;
                ytab[1] = c;
                ytab[2] = b;
                ytab[3] = c;
                ytab[4] = e;
                ytab[5] = d;
                ytab[6] = e;
                ytab[7] = g;
                ytab[8] = f;
                ytab[9] = s;
                ytab[10] = r;
                ytab[11] = s;
                ytab[12] = u;
                ytab[13] = v;
                ytab[14] = t;
                ytab[15] = t;
                ytab[16] = v;
                ytab[17] = t;
                ytab[18] = u;

                weight[0] = w1;
                weight[1] = w2;
                weight[2] = w2;
                weight[3] = w2;
                weight[4] = w3;
                weight[5] = w3;
                weight[6] = w3;
                weight[7] = w4;
                weight[8] = w4;
                weight[9] = w4;
                weight[10] = w5;
                weight[11] = w5;
                weight[12] = w5;
                weight[13] = w6;
                weight[14] = w6;
                weight[15] = w6;
                weight[16] = w6;
                weight[17] = w6;
                weight[18] = w6;
                break;
            //
            //  28 points, precision 11.
            //  Lyness and Jesperson.
            //
            case 19:
                a = 1.0 / 3.0;
                b = 0.9480217181434233;
                c = 0.02598914092828833;
                d = 0.8114249947041546;
                e = 0.09428750264792270;
                f = 0.01072644996557060;
                g = 0.4946367750172147;
                p = 0.5853132347709715;
                q = 0.2073433826145142;
                r = 0.1221843885990187;
                s = 0.4389078057004907;
                t = 0.6779376548825902;
                u = 0.04484167758913055;
                v = 0.27722066752827925;
                w = 0.8588702812826364;
                double x = 0.0;
                double y = 0.1411297187173636;

                w1 = 0.08797730116222190;
                w2 = 0.008744311553736190;
                w3 = 0.03808157199393533;
                w4 = 0.01885544805613125;
                w5 = 0.07215969754474100;
                w6 = 0.06932913870553720;
                w7 = 0.04105631542928860;
                w8 = 0.007362383783300573;

                xtab[0] = a;
                xtab[1] = b;
                xtab[2] = c;
                xtab[3] = c;
                xtab[4] = d;
                xtab[5] = e;
                xtab[6] = e;
                xtab[7] = f;
                xtab[8] = g;
                xtab[9] = g;
                xtab[10] = p;
                xtab[11] = q;
                xtab[12] = q;
                xtab[13] = r;
                xtab[14] = s;
                xtab[15] = s;
                xtab[16] = t;
                xtab[17] = t;
                xtab[18] = u;
                xtab[19] = u;
                xtab[20] = v;
                xtab[21] = v;
                xtab[22] = w;
                xtab[23] = w;
                xtab[24] = x;
                xtab[25] = x;
                xtab[26] = y;
                xtab[27] = y;

                ytab[0] = a;
                ytab[1] = c;
                ytab[2] = b;
                ytab[3] = c;
                ytab[4] = e;
                ytab[5] = d;
                ytab[6] = e;
                ytab[7] = g;
                ytab[8] = f;
                ytab[9] = g;
                ytab[10] = q;
                ytab[11] = p;
                ytab[12] = q;
                ytab[13] = s;
                ytab[14] = r;
                ytab[15] = s;
                ytab[16] = u;
                ytab[17] = v;
                ytab[18] = t;
                ytab[19] = v;
                ytab[20] = t;
                ytab[21] = u;
                ytab[22] = x;
                ytab[23] = y;
                ytab[24] = w;
                ytab[25] = y;
                ytab[26] = w;
                ytab[27] = x;

                weight[0] = w1;
                weight[1] = w2;
                weight[2] = w2;
                weight[3] = w2;
                weight[4] = w3;
                weight[5] = w3;
                weight[6] = w3;
                weight[7] = w4;
                weight[8] = w4;
                weight[9] = w4;
                weight[10] = w5;
                weight[11] = w5;
                weight[12] = w5;
                weight[13] = w6;
                weight[14] = w6;
                weight[15] = w6;
                weight[16] = w7;
                weight[17] = w7;
                weight[18] = w7;
                weight[19] = w7;
                weight[20] = w7;
                weight[21] = w7;
                weight[22] = w8;
                weight[23] = w8;
                weight[24] = w8;
                weight[25] = w8;
                weight[26] = w8;
                weight[27] = w8;
                break;
            //
            //  37 points, precision 13.
            //
            case 20:
                a = 1.0 / 3.0;
                b = 0.950275662924105565450352089520;
                c = 0.024862168537947217274823955239;
                d = 0.171614914923835347556304795551;
                e = 0.414192542538082326221847602214;
                f = 0.539412243677190440263092985511;
                g = 0.230293878161404779868453507244;

                w1 = 0.051739766065744133555179145422;
                w2 = 0.008007799555564801597804123460;
                w3 = 0.046868898981821644823226732071;
                w4 = 0.046590940183976487960361770070;
                w5 = 0.031016943313796381407646220131;
                w6 = 0.010791612736631273623178240136;
                w7 = 0.032195534242431618819414482205;
                w8 = 0.015445834210701583817692900053;
                double w9 = 0.017822989923178661888748319485;
                double wx = 0.037038683681384627918546472190;

                xtab[0] = a;
                xtab[1] = b;
                xtab[2] = c;
                xtab[3] = c;
                xtab[4] = d;
                xtab[5] = e;
                xtab[6] = e;
                xtab[7] = f;
                xtab[8] = g;
                xtab[9] = g;

                ytab[0] = a;
                ytab[1] = c;
                ytab[2] = b;
                ytab[3] = c;
                ytab[4] = e;
                ytab[5] = d;
                ytab[6] = e;
                ytab[7] = g;
                ytab[8] = f;
                ytab[9] = g;

                weight[0] = w1;
                weight[1] = w2;
                weight[2] = w2;
                weight[3] = w2;
                weight[4] = w3;
                weight[5] = w3;
                weight[6] = w3;
                weight[7] = w4;
                weight[8] = w4;
                weight[9] = w4;
                weight[10] = w5;
                weight[11] = w5;
                weight[12] = w5;
                weight[13] = w6;
                weight[14] = w6;
                weight[15] = w6;
                weight[16] = w7;
                weight[17] = w7;
                weight[18] = w7;
                weight[19] = w8;
                weight[20] = w8;
                weight[21] = w8;
                weight[22] = w8;
                weight[23] = w8;
                weight[24] = w8;
                weight[25] = w9;
                weight[26] = w9;
                weight[27] = w9;
                weight[28] = w9;
                weight[29] = w9;
                weight[30] = w9;
                weight[31] = wx;
                weight[32] = wx;
                weight[33] = wx;
                weight[34] = wx;
                weight[35] = wx;
                weight[36] = wx;

                a = 0.772160036676532561750285570113;
                b = 0.113919981661733719124857214943;

                xtab[10] = a;
                ytab[10] = b;

                xtab[11] = b;
                ytab[11] = a;

                xtab[12] = b;
                ytab[12] = b;

                a = 0.009085399949835353883572964740;
                b = 0.495457300025082323058213517632;

                xtab[13] = a;
                ytab[13] = b;

                xtab[14] = b;
                ytab[14] = a;

                xtab[15] = b;
                ytab[15] = b;

                a = 0.062277290305886993497083640527;
                b = 0.468861354847056503251458179727;

                xtab[16] = a;
                ytab[16] = b;

                xtab[17] = b;
                ytab[17] = a;

                xtab[18] = b;
                ytab[18] = b;

                a = 0.022076289653624405142446876931;
                b = 0.851306504174348550389457672223;
                c = 0.126617206172027096933163647918263;

                xtab[19] = a;
                ytab[19] = b;

                xtab[20] = a;
                ytab[20] = c;

                xtab[21] = b;
                ytab[21] = a;

                xtab[22] = b;
                ytab[22] = c;

                xtab[23] = c;
                ytab[23] = a;

                xtab[24] = c;
                ytab[24] = b;

                a = 0.018620522802520968955913511549;
                b = 0.689441970728591295496647976487;
                c = 0.291937506468887771754472382212953;

                xtab[25] = a;
                ytab[25] = b;

                xtab[26] = a;
                ytab[26] = c;

                xtab[27] = b;
                ytab[27] = a;

                xtab[28] = b;
                ytab[28] = c;

                xtab[29] = c;
                ytab[29] = a;

                xtab[30] = c;
                ytab[30] = b;

                a = 0.096506481292159228736516560903;
                b = 0.635867859433872768286976979827;
                c = 0.267625659273967961282458816185681;

                xtab[31] = a;
                ytab[31] = b;

                xtab[32] = a;
                ytab[32] = c;

                xtab[33] = b;
                ytab[33] = a;

                xtab[34] = b;
                ytab[34] = c;

                xtab[35] = c;
                ytab[35] = a;

                xtab[36] = c;
                ytab[36] = b;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("TRIANGLE_UNIT_SET - Fatal error!");
                Console.WriteLine("  Illegal value of RULE = " + rule + "");
                break;
        }
    }

    public static int triangle_unit_size(int rule)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_SIZE returns the "size" of a unit triangle quadrature rule.
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
        //    Jarle Berntsen, Terje Espelid,
        //    Algorithm 706,
        //    DCUTRI: an algorithm for adaptive cubature over a collection of triangles, 
        //    ACM Transactions on Mathematical Software,
        //    Volume 18, Number 3, September 1992, pages 329-342.
        //
        //    Elise deDoncker, Ian Robinson,
        //    Algorithm 612:
        //    Integration over a Triangle Using Nonlinear Extrapolation,
        //    ACM Transactions on Mathematical Software,
        //    Volume 10, Number 1, March 1984, pages 17-22.
        //
        //    DP Laurie,
        //    Algorithm 584,
        //    CUBTRI, Automatic Cubature Over a Triangle,
        //    ACM Transactions on Mathematical Software,
        //    Volume 8, Number 2, 1982, pages 210-218.
        //
        //    James Lyness, Dennis Jespersen,
        //    Moderate Degree Symmetric Quadrature Rules for the Triangle,
        //    Journal of the Institute of Mathematics and its Applications,
        //    Volume 15, Number 1, February 1975, pages 19-32.
        //
        //    Hans Rudolf Schwarz,
        //    Methode der Finiten Elemente,
        //    Teubner Studienbuecher, 1980,
        //    ISBN: 3-519-02349-0.
        //
        //    Gilbert Strang, George Fix,
        //    An Analysis of the Finite Element Method,
        //    Prentice Hall, 1973, page 184,
        //    ISBN: 096140888X,
        //    LC: TA335.S77.
        //
        //    Arthur Stroud,
        //    Approximate Calculation of Multiple Integrals,
        //    Prentice Hall, 1971,
        //    ISBN: 0130438936,
        //    LC: QA311.S85.
        //
        //    Olgierd Zienkiewicz,
        //    The Finite Element Method,
        //    Sixth Edition,
        //    Butterworth-Heinemann, 2005,
        //    ISBN: 0750663200,
        //    TA640.2.Z54
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //     1, ORDER =  1, precision 1, Zienkiewicz #1.
        //     2, ORDER =  2, precision 1, (the "vertex rule").
        //     3, ORDER =  3, precision 2, Strang and Fix formula #1.
        //     4, ORDER =  3, precision 2, Strang and Fix formula #2, Zienkiewicz #2.
        //     5, ORDER =  4, precision 3, Strang and Fix formula #3, Zienkiewicz #3.
        //     6, ORDER =  6, precision 3, Strang and Fix formula #4.
        //     7, ORDER =  6, precision 3, Stroud formula T2:3-1.
        //     8, ORDER =  6, precision 4, Strang and Fix formula #5.
        //     9, ORDER =  7, precision 4, Strang and Fix formula #6.
        //    10, ORDER =  7, precision 5, Strang and Fix formula #7,
        //        Stroud formula T2:5-1, Zienkiewicz #4, Schwarz Table 2.2.
        //    11, ORDER =  9, precision 6, Strang and Fix formula #8.
        //    12, ORDER = 12, precision 6, Strang and Fix formula #9.
        //    13, ORDER = 13, precision 7, Strang and Fix formula #10.
        //    14, ORDER =  7, precision ?.
        //    15, ORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
        //    16, ORDER = 64, precision 15, triangular product Gauss rule.
        //    17, ORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
        //    18, ORDER = 19, precision 9, from TRIEX, Lyness and Jespersen.
        //    19, ORDER = 28, precision 11, from TRIEX, Lyness and Jespersen.
        //    20, ORDER = 37, precision 13, from ACM TOMS #706.
        //
        //    Output, int TRIANGLE_UNIT_SIZE, the order of the rule.
        //
    {
        int size;

        switch (rule)
        {
            case 1:
                size = 1;
                break;
            case 2:
            case 3:
            case 4:
                size = 3;
                break;
            case 5:
                size = 4;
                break;
            case 6:
            case 7:
            case 8:
                size = 6;
                break;
            case 9:
            case 10:
                size = 7;
                break;
            case 11:
                size = 9;
                break;
            case 12:
                size = 12;
                break;
            case 13:
                size = 13;
                break;
            case 14:
                size = 7;
                break;
            case 15:
                size = 16;
                break;
            case 16:
                size = 64;
                break;
            case 17:
            case 18:
                size = 19;
                break;
            case 19:
                size = 28;
                break;
            case 20:
                size = 37;
                break;
            default:
                size = -1;
                break;
        }

        return size;
    }

    public static double triangle_unit_sum(int setting, Func<int, double, double, double> func, int order,
            double[] xtab, double[] ytab, double[] weight)

        //****************************************************************************80
        //
        //// TRIANGLE_UNIT_SUM carries out a quadrature rule in the unit triangle.
        //
        //  Integration region:
        //
        //      0 <= X,
        //    and
        //      0 <= Y, 
        //    and
        //      X + Y <= 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, Func < double, double, double > func, the name of the
        //    user supplied function to be integrated.
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
            quad += weight[i] * func(setting, xtab[i], ytab[i]);
        }

        double volume = triangle_unit_volume();
        double result = quad * volume;

        return result;
    }

    public static double triangle_unit_volume()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_VOLUME returns the "volume" of the unit triangle in 2D.
        //
        //  Integration region:
        //
        //      0 <= X,
        //    and
        //      0 <= Y, 
        //    and
        //      X + Y <= 1.
        //
        //  Discussion:
        //
        //    The "volume" of a triangle is usually called its area.
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
        //    Output, double TRIANGLE_UNIT_VOLUME, the volume of the unit
        //    triangle.
        //
    {
        const double volume = 1.0 / 2.0;

        return volume;
    }

    public static double triangle_volume(double[] x, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_VOLUME returns the "volume" of a triangle in 2D.
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
        //    Input, double X[3], Y[3], the vertices of the triangle.
        //
        //    Output, double TRIANGLE_VOLUME, the volume of the triangle.
        //
    {
        double value = 0.5 * Math.Abs(
            x[0] * (y[1] - y[2]) +
            x[1] * (y[2] - y[0]) +
            x[2] * (y[0] - y[1]));

        return value;
    }


}