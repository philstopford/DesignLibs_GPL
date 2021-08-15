using System;
using Burkardt.Quadrature;

namespace Burkardt.Stroud
{
    public static class Square
    {
        public static double square_sum(Func<double, double, double> func, double[] center,
                double r, int order, double[] xtab, double[] ytab, double[] weight)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SQUARE_SUM carries out a quadrature rule over a square.
            //
            //  Integration region:
            //
            //      abs ( X - CENTER(1) ) <= R 
            //    and
            //      abs ( Y - CENTER(2) ) <= R
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
            //    Input, Func < double, double, double > func, the name of the function 
            //    to be integrated.  
            //
            //    Input, double[] center, the center of the square.
            //
            //    Input, double R, the radius of the square.
            //
            //    Input, int ORDER, the order of the rule.
            //
            //    Input, double XTAB[ORDER], YTAB[ORDER], the abscissas of 
            //    the rule.
            //
            //    Input, double WEIGHT[ORDER], the weights of the rule.
            //
            //    Output, double SQUARE_SUM, the approximate integral of the function.
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
                quad = quad + 0.25 * weight[i] * func(x, y);
            }

            volume = 4.0 * r * r;
            result = quad * volume;

            return result;
        }

        public static void square_unit_set(int rule, int order, double[] xtab, double[] ytab,
                double[] weight)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SQUARE_UNIT_SET sets quadrature weights and abscissas in the unit square.
            //
            //  Discussion;
            //
            //    To get the value of ORDER associated with a given rule, 
            //    call SQUARE_UNIT_SIZE first.
            //
            //  Integration region:
            //
            //      -1 <= X <= 1,
            //    and
            //      -1 <= Y <= 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
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
            //  Parameters:
            //
            //    Input, int RULE, the rule number.
            //    1, order 1, degree 1 rule.
            //    2, order 4, degree 3, rule.
            //    3, order 9, degree 5 rule.
            //    4, order 12 degree 7 rule, Stroud number C2:7-1.
            //    5, order 13 degree 7 rule, Stroud number C2:7-3.
            //    6, order 64 degree 15 product rule.
            //
            //    Input, int ORDER, the order of the rule.
            //
            //    Output, double XTAB[ORDER], YTAB[ORDER], the abscissas.
            //
            //    Output, double WEIGHT[ORDER], the weights.
            //
        {
            double a;
            double c;
            int i;
            int j;
            int k;
            int order2 = 8;
            double r;
            double s;
            double t;
            double w1;
            double w2;
            double w3;
            double[] weight2;
            double[] xtab2;
            double z;

            if (rule == 1)
            {
                weight[0] = 4.0;

                xtab[0] = 0.0;
                ytab[0] = 0.0;
            }
            else if (rule == 2)
            {
                a = 1.0;
                s = 1.0 / Math.Sqrt(3.0);

                xtab[0] = -s;
                xtab[1] = s;
                xtab[2] = -s;
                xtab[3] = s;

                ytab[0] = -s;
                ytab[1] = -s;
                ytab[2] = s;
                ytab[3] = s;

                weight[0] = a;
                weight[1] = a;
                weight[2] = a;
                weight[3] = a;
            }
            else if (rule == 3)
            {
                s = Math.Sqrt(0.6);
                z = 0.0;
                w1 = 64.0 / 81.0;
                w2 = 25.0 / 81.0;
                w3 = 40.0 / 81.0;

                xtab[0] = z;
                xtab[1] = -s;
                xtab[2] = s;
                xtab[3] = -s;
                xtab[4] = s;
                xtab[5] = z;
                xtab[6] = -s;
                xtab[7] = s;
                xtab[8] = z;

                ytab[0] = z;
                ytab[1] = -s;
                ytab[2] = -s;
                ytab[3] = s;
                ytab[4] = s;
                ytab[5] = -s;
                ytab[6] = z;
                ytab[7] = z;
                ytab[8] = s;

                weight[0] = w1;
                weight[1] = w2;
                weight[2] = w2;
                weight[3] = w2;
                weight[4] = w2;
                weight[5] = w3;
                weight[6] = w3;
                weight[7] = w3;
                weight[8] = w3;
            }
            else if (rule == 4)
            {
                r = Math.Sqrt(6.0 / 7.0);
                c = 3.0 * Math.Sqrt(583.0);
                s = Math.Sqrt((114.0 - c) / 287.0);
                t = Math.Sqrt((114.0 + c) / 287.0);
                w1 = 4.0 * 49.0 / 810.0;
                w2 = 4.0 * (178981.0 + 923.0 * c) / 1888920.0;
                w3 = 4.0 * (178981.0 - 923.0 * c) / 1888920.0;
                z = 0.0;

                xtab[0] = r;
                xtab[1] = z;
                xtab[2] = -r;
                xtab[3] = z;
                xtab[4] = s;
                xtab[5] = -s;
                xtab[6] = -s;
                xtab[7] = s;
                xtab[8] = t;
                xtab[9] = -t;
                xtab[10] = -t;
                xtab[11] = t;

                ytab[0] = z;
                ytab[1] = r;
                ytab[2] = z;
                ytab[3] = -r;
                ytab[4] = s;
                ytab[5] = s;
                ytab[6] = -s;
                ytab[7] = -s;
                ytab[8] = t;
                ytab[9] = t;
                ytab[10] = -t;
                ytab[11] = -t;

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
            }
            else if (rule == 5)
            {
                r = Math.Sqrt(12.0 / 35.0);
                c = 3.0 * Math.Sqrt(186.0);
                s = Math.Sqrt((93.0 + c) / 155.0);
                t = Math.Sqrt((93.0 - c) / 155.0);
                w1 = 8.0 / 162.0;
                w2 = 98.0 / 162.0;
                w3 = 31.0 / 162.0;
                z = 0.0;

                xtab[0] = z;
                xtab[1] = r;
                xtab[2] = -r;
                xtab[3] = z;
                xtab[4] = z;
                xtab[5] = s;
                xtab[6] = s;
                xtab[7] = -s;
                xtab[8] = -s;
                xtab[9] = t;
                xtab[10] = t;
                xtab[11] = -t;
                xtab[12] = -t;

                ytab[0] = z;
                ytab[1] = z;
                ytab[2] = z;
                ytab[3] = r;
                ytab[4] = -r;
                ytab[5] = t;
                ytab[6] = -t;
                ytab[7] = t;
                ytab[8] = -t;
                ytab[9] = s;
                ytab[10] = -s;
                ytab[11] = s;
                ytab[12] = -s;

                weight[0] = w1;
                weight[1] = w2;
                weight[2] = w2;
                weight[3] = w2;
                weight[4] = w2;
                weight[5] = w3;
                weight[6] = w3;
                weight[7] = w3;
                weight[8] = w3;
                weight[9] = w3;
                weight[10] = w3;
                weight[11] = w3;
                weight[12] = w3;
            }
            else if (rule == 6)
            {
                xtab2 = new double[order2];
                weight2 = new double[order2];

                LegendreQuadrature.legendre_set(order2, ref xtab2, ref weight2);

                k = 0;

                for (i = 0; i < order2; i++)
                {
                    for (j = 0; j < order2; j++)
                    {
                        xtab[k] = xtab2[i];
                        ytab[k] = xtab2[j];
                        weight[k] = weight2[i] * weight2[j];
                        k = k + 1;
                    }
                }
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("SQUARE_UNIT_SET - Fatal error!");
                Console.WriteLine("  Illegal value of RULE = " + rule + "");
            }
        }

        public static int square_unit_size(int rule)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SQUARE_UNIT_SIZE sizes a quadrature rule in the unit square.
            //
            //  Integration region:
            //
            //      -1 <= X <= 1,
            //    and
            //      -1 <= Y <= 1.
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
            //  Reference:
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
            //  Parameters:
            //
            //    Input, int RULE, the rule number.
            //    1, a 1 point 1st degree rule.
            //    2, a 4 point 3rd degree rule.
            //    3, a 9 point 5th degree rule.
            //    4, a 12 point 7-th degree rule, Stroud number C2:7-1.
            //    5, a 13 point 7-th degree rule, Stroud number C2:7-3.
            //    6, a 64 point 15-th degree product rule.
            //
            //    Output, int SQUARE_UNIT_SIZE, the order of the rule.
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
                order = 9;
            }
            else if (rule == 4)
            {
                order = 12;
            }
            else if (rule == 5)
            {
                order = 13;
            }
            else if (rule == 6)
            {
                order = 64;
            }
            else
            {
                order = -1;
            }

            return order;
        }

        public static double square_unit_sum(Func<double, double, double> func, int order,
                double[] xtab, double[] ytab, double[] weight)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SQUARE_UNIT_SUM carries out a quadrature rule over the unit square.
            //
            //  Integration region:
            //
            //      -1 <= X <= 1, 
            //    and
            //      -1 <= Y <= 1.
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
            //    Input, Func < double, double, double > func, the name of the function
            //    to be integrated. 
            //
            //    Input, int ORDER, the order of the rule.
            //
            //    Input, double XTAB[ORDER], YTAB[ORDER], the abscissas.
            //
            //    Input, double WEIGHT[ORDER], the weights.
            //
            //    Output, double SQUARE_UNIT_SUM, the approximate integral of the function.
            //
        {
            int i;
            double quad;
            double result;
            double volume;

            quad = 0.0;
            for (i = 0; i < order; i++)
            {
                quad = quad + weight[i] * func(xtab[i], ytab[i]) / 4.0;
            }

            volume = 1.0;
            result = quad * volume;

            return result;
        }


    }
}