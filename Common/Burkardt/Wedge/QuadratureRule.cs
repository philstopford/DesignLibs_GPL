using System;
using Burkardt.Types;

namespace Burkardt.Wedge
{
    public static class QuadratureRule
    {
        public static double wedge_integral ( int[] expon )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEDGE_INTEGRAL: monomial integral in a unit wedge.
        //
        //  Discussion:
        //
        //    This routine returns the integral of
        //
        //      product ( 1 <= I <= 3 ) X(I)^EXPON(I)
        //
        //    over the unit wedge.
        //
        //    The integration region is:
        //
        //      0 <= X
        //      0 <= Y
        //      X + Y <= 1
        //      -1 <= Z <= 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 March 2008
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
        //    Input, int EXPON[3], the exponents.
        //
        //    Output, double WEDGE_INTEGRAL, the integral of the monomial.
        //
        {
            int i;
            int k;
            double value;
            //
            //  The first computation ends with VALUE = 1.0;
            //
            value = 1.0;

            k = expon[0];

            for ( i = 1; i <= expon[1]; i++ )
            {
                k = k + 1;
                value = value * ( double ) ( i ) / ( double ) ( k );
            }

            k = k + 1;
            value = value / ( double ) ( k );

            k = k + 1;
            value = value / ( double ) ( k );
            //
            //  Now account for integration in Z.
            //
            if ( expon[2] == - 1 )
            {
                Console.WriteLine("");
                Console.WriteLine("WEDGE_INTEGRAL - Fatal error!");
                Console.WriteLine("  EXPON[2] = -1 is not a legal input.");
                return ( 1 );
            }
            else if ( ( expon[2] % 2 ) == 1 )
            {
                value = 0.0;
            }
            else
            {
                value = value * 2.0 / ( double ) ( expon[2] + 1 );
            }

            return value;
        }
        
        public static void wedge_rule(int line_order, int triangle_order, ref double[] w, ref double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    WEDGE_RULE returns a quadrature rule for the unit wedge.
            //
            //  Discussion:
            //
            //    It is usually sensible to take LINE_ORDER and TRIG_ORDER so that
            //    the line and triangle rules are roughly the same precision.  For that
            //    criterion, we recommend the following combinations:
            //
            //      TRIANGLE_ORDER  LINE_ORDER  Precision
            //      --------------  ----------  ---------
            //          1               1       1 x 1
            //          3               2       2 x 3
            //         -3               2       2 x 3
            //          6               3       4 x 5
            //         -6               2       3 x 3
            //          7               3       5 x 5
            //         12               4       6 x 7
            //
            //    The integration region is:
            //
            //      0 <= X
            //      0 <= Y
            //      X + Y <= 1
            //      -1 <= Z <= 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 April 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carlos Felippa,
            //    A compendium of FEM integration formulas for symbolic work,
            //    Engineering Computation,
            //    Volume 21, Number 8, 2004, pages 867-890.
            //
            //  Parameters:
            //
            //    Input, int LINE_ORDER, the index of the line rule.
            //    The index of the rule is equal to the order of the rule.
            //    1 <= LINE_ORDER <= 5.
            //
            //    Input, int TRIANGLE_ORDER, the indes of the triangle rule.
            //    The index of the rule is 1, 3, -3, 6, -6, 7 or 12.
            //
            //    Output, double W[LINE_ORDER*abs(TRIANGLE_ORDER)], the weights.
            //
            //    Output, double XYZ[3*LINE_ORDER*abs(TRIANGLE_ORDER)], the abscissas.
            //
        {
            int i;
            int j;
            int k;
            double[] line_w;
            double[] line_x;
            double[] triangle_w;
            double[] triangle_xy;

            line_w = new double[line_order];
            line_x = new double[line_order];

            if (line_order == 1)
            {
                line_o01(ref line_w, ref line_x);
            }
            else if (line_order == 2)
            {
                line_o02(ref line_w, ref line_x);
            }
            else if (line_order == 3)
            {
                line_o03(ref line_w, ref line_x);
            }
            else if (line_order == 4)
            {
                line_o04(ref line_w, ref line_x);
            }
            else if (line_order == 5)
            {
                line_o05(ref line_w, ref line_x);
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("WEDGE_RULE - Fatal error!");
                Console.WriteLine("  Illegal value of LINE_ORDER.");
                return;
            }

            triangle_w = new double[Math.Abs(triangle_order)];
            triangle_xy = new double[2 * Math.Abs(triangle_order)];

            if (triangle_order == 1)
            {
                triangle_o01(ref triangle_w, ref triangle_xy);
            }
            else if (triangle_order == 3)
            {
                triangle_o03(ref triangle_w, ref triangle_xy);
            }
            else if (triangle_order == -3)
            {
                triangle_o03b(ref triangle_w, ref triangle_xy);
            }
            else if (triangle_order == 6)
            {
                triangle_o06(ref triangle_w, ref triangle_xy);
            }
            else if (triangle_order == -6)
            {
                triangle_o06b(ref triangle_w, ref triangle_xy);
            }
            else if (triangle_order == 7)
            {
                triangle_o07(ref triangle_w, ref triangle_xy);
            }
            else if (triangle_order == 12)
            {
                triangle_o12(ref triangle_w, ref triangle_xy);
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("WEDGE_RULE - Fatal error!");
                Console.WriteLine("  Illegal value of TRIANGLE_ORDER.");
                return;
            }

            k = 0;
            for (i = 0; i < line_order; i++)
            {
                for (j = 0; j < Math.Abs(triangle_order); j++)
                {
                    w[k] = line_w[i] * triangle_w[j];
                    xyz[0 + k * 3] = triangle_xy[0 + j * 2];
                    xyz[1 + k * 3] = triangle_xy[1 + j * 2];
                    xyz[2 + k * 3] = line_x[i];
                    k = k + 1;
                }
            }
        }

        public static void line_o01(ref double[] w, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_O01 returns a 1 point quadrature rule for the unit line.
        //
        //  Discussion:
        //
        //    The integration region is:
        //
        //    - 1.0 <= X <= 1.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 April 2009
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
        //    Output, double W[1], the weights.
        //
        //    Output, double X[1], the abscissas.
        //
        {
            int order = 1;
            double[] w_save =  {
                1.0
            }
            ;
            double[] x_save =  {
                0.0
            }
            ;

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(order, x_save, ref x);
        }

        public static void line_o02(ref double[] w, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_O02 returns a 2 point quadrature rule for the unit line.
        //
        //  Discussion:
        //
        //    The integration region is:
        //
        //    - 1.0 <= X <= 1.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 April 2009
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
        //    Output, double W[2], the weights.
        //
        //    Output, double X[2], the abscissas.
        //
        {
            int order = 2;
            double[] w_save =  {
                0.5,
                0.5
            }
            ;
            double[] x_save =  {
                -0.57735026918962576451,
                0.57735026918962576451
            }
            ;

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(order, x_save, ref x);
        }

        public static void line_o03(ref double[] w, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_O03 returns a 3 point quadrature rule for the unit line.
        //
        //  Discussion:
        //
        //    The integration region is:
        //
        //    - 1.0 <= X <= 1.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 April 2009
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
        //    Output, double W[3], the weights.
        //
        //    Output, double X[3], the abscissas.
        //
        {
            int order = 3;
            double[] w_save =  {
                0.27777777777777777777,
                0.44444444444444444444,
                0.27777777777777777777
            }
            ;
            double[] x_save =  {
                -0.77459666924148337704,
                0.00000000000000000000,
                0.77459666924148337704
            }
            ;

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(order, x_save, ref x);
        }

        public static void line_o04(ref double[] w, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_O04 returns a 4 point quadrature rule for the unit line.
        //
        //  Discussion:
        //
        //    The integration region is:
        //
        //    - 1.0 <= X <= 1.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 April 2009
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
        //    Output, double W[4], the weights.
        //
        //    Output, double X[4], the abscissas.
        //
        {
            int order = 4;
            double[] w_save =  {
                0.173927422568727,
                0.326072577431273,
                0.326072577431273,
                0.173927422568727
            }
            ;
            double[] x_save =  {
                -0.86113631159405257522,
                -0.33998104358485626480,
                0.33998104358485626480,
                0.86113631159405257522
            }
            ;

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(order, x_save, ref x);
        }

        public static void line_o05(ref double[] w, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_O05 returns a 5 point quadrature rule for the unit line.
        //
        //  Discussion:
        //
        //    The integration region is:
        //
        //    - 1.0 <= X <= 1.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 April 2009
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
        //    Output, double W[5], the weights.
        //
        //    Output, double X[5], the abscissas.
        //
        {
            int order = 5;
            double[] w_save =  {
                0.118463442528095,
                0.239314335249683,
                0.284444444444444,
                0.239314335249683,
                0.118463442528095
            }
            ;
            double[] x_save =  {
                -0.90617984593866399280,
                -0.53846931010568309104,
                0.00000000000000000000,
                0.53846931010568309104,
                0.90617984593866399280
            }
            ;

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(order, x_save, ref x);

        }

        public static void triangle_o01(ref double[] w, ref double[] xy)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_O01 returns a 1 point quadrature rule for the unit triangle.
            //
            //  Discussion:
            //
            //    This rule is precise for monomials through degree 1.
            //
            //    The integration region is:
            //
            //      0 <= X
            //      0 <= Y
            //      X + Y <= 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 April 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carlos Felippa,
            //    A compendium of FEM integration formulas for symbolic work,
            //    Engineering Computation,
            //    Volume 21, Number 8, 2004, pages 867-890.
            //
            //  Parameters:
            //
            //    Output, double W[1], the weights.
            //
            //    Output, double XY[2*1], the abscissas.
            //
        {
            int order = 1;

            double[] w_save =  {
                1.0
            }
            ;
            double[] xy_save =  {
                0.33333333333333333333, 0.33333333333333333333
            }
            ;

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(2 * order, xy_save, ref xy);
        }

        public static void triangle_o03(ref double[] w, ref double[] xy)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_O03 returns a 3 point quadrature rule for the unit triangle.
            //
            //  Discussion:
            //
            //    This rule is precise for monomials through degree 2.
            //
            //    The integration region is:
            //
            //      0 <= X
            //      0 <= Y
            //      X + Y <= 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 April 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carlos Felippa,
            //    A compendium of FEM integration formulas for symbolic work,
            //    Engineering Computation,
            //    Volume 21, Number 8, 2004, pages 867-890.
            //
            //  Parameters:
            //
            //    Output, double W[3], the weights.
            //
            //    Output, double XY[2*3], the abscissas.
            //
        {
            int order = 3;

            double[] w_save =  {
                0.33333333333333333333,
                0.33333333333333333333,
                0.33333333333333333333
            }
            ;
            double[] xy_save =  {
                0.66666666666666666667, 0.16666666666666666667,
                0.16666666666666666667, 0.66666666666666666667,
                0.16666666666666666667, 0.16666666666666666667
            }
            ;

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(2 * order, xy_save, ref xy);
        }

        public static void triangle_o03b(ref double[] w, ref double[] xy)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_O03B returns a 3 point quadrature rule for the unit triangle.
            //
            //  Discussion:
            //
            //    This rule is precise for monomials through degree 2.
            //
            //    The integration region is:
            //
            //      0 <= X
            //      0 <= Y
            //      X + Y <= 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 April 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carlos Felippa,
            //    A compendium of FEM integration formulas for symbolic work,
            //    Engineering Computation,
            //    Volume 21, Number 8, 2004, pages 867-890.
            //
            //  Parameters:
            //
            //    Output, double W[3], the weights.
            //
            //    Output, double XY[2*3], the abscissas.
            //
        {
            int order = 3;

            double[] w_save =  {
                0.33333333333333333333,
                0.33333333333333333333,
                0.33333333333333333333
            }
            ;
            double[] xy_save =  {
                0.0, 0.5,
                0.5, 0.0,
                0.5, 0.5
            }
            ;

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(2 * order, xy_save, ref xy);
        }

        public static void triangle_o06(ref double[] w, ref double[] xy)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_O06 returns a 6 point quadrature rule for the unit triangle.
            //
            //  Discussion:
            //
            //    This rule is precise for monomials through degree 4.
            //
            //    The integration region is:
            //
            //      0 <= X
            //      0 <= Y
            //      X + Y <= 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 April 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carlos Felippa,
            //    A compendium of FEM integration formulas for symbolic work,
            //    Engineering Computation,
            //    Volume 21, Number 8, 2004, pages 867-890.
            //
            //  Parameters:
            //
            //    Output, double W[6], the weights.
            //
            //    Output, double XY[2*6], the abscissas.
            //
        {
            int order = 6;

            double[] w_save =
            {
                0.22338158967801146570,
                0.22338158967801146570,
                0.22338158967801146570,
                0.10995174365532186764,
                0.10995174365532186764,
                0.10995174365532186764
            };
            double[] xy_save =
            {
                0.10810301816807022736, 0.44594849091596488632,
                0.44594849091596488632, 0.10810301816807022736,
                0.44594849091596488632, 0.44594849091596488632,
                0.81684757298045851308, 0.091576213509770743460,
                0.091576213509770743460, 0.81684757298045851308,
                0.091576213509770743460, 0.091576213509770743460
            };

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(2 * order, xy_save, ref xy);
        }

        public static void triangle_o06b(ref double[] w, ref double[] xy)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_O06B returns a 6 point quadrature rule for the unit triangle.
            //
            //  Discussion:
            //
            //    This rule is precise for monomials through degree 3.
            //
            //    The integration region is:
            //
            //      0 <= X
            //      0 <= Y
            //      X + Y <= 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 April 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carlos Felippa,
            //    A compendium of FEM integration formulas for symbolic work,
            //    Engineering Computation,
            //    Volume 21, Number 8, 2004, pages 867-890.
            //
            //  Parameters:
            //
            //    Output, double W[6], the weights.
            //
            //    Output, double XY[2*6], the abscissas.
            //
        {
            int order = 6;

            double[] w_save =
            {
                0.30000000000000000000,
                0.30000000000000000000,
                0.30000000000000000000,
                0.033333333333333333333,
                0.033333333333333333333,
                0.033333333333333333333
            };
            double[] xy_save =
            {
                0.66666666666666666667, 0.16666666666666666667,
                0.16666666666666666667, 0.66666666666666666667,
                0.16666666666666666667, 0.16666666666666666667,
                0.0, 0.5,
                0.5, 0.0,
                0.5, 0.5
            };

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(2 * order, xy_save, ref xy);
        }

        public static void triangle_o07(ref double[] w, ref double[] xy)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_O07 returns a 7 point quadrature rule for the unit triangle.
            //
            //  Discussion:
            //
            //    This rule is precise for monomials through degree 5.
            //
            //    The integration region is:
            //
            //      0 <= X
            //      0 <= Y
            //      X + Y <= 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 April 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carlos Felippa,
            //    A compendium of FEM integration formulas for symbolic work,
            //    Engineering Computation,
            //    Volume 21, Number 8, 2004, pages 867-890.
            //
            //  Parameters:
            //
            //    Output, double W[7], the weights.
            //
            //    Output, double XY[2*7], the abscissas.
            //
        {
            int order = 7;

            double[] w_save =
            {
                0.12593918054482715260,
                0.12593918054482715260,
                0.12593918054482715260,
                0.13239415278850618074,
                0.13239415278850618074,
                0.13239415278850618074,
                0.22500000000000000000
            };
            double[] xy_save =
            {
                0.79742698535308732240, 0.10128650732345633880,
                0.10128650732345633880, 0.79742698535308732240,
                0.10128650732345633880, 0.10128650732345633880,
                0.059715871789769820459, 0.47014206410511508977,
                0.47014206410511508977, 0.059715871789769820459,
                0.47014206410511508977, 0.47014206410511508977,
                0.33333333333333333333, 0.33333333333333333333
            };

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(2 * order, xy_save, ref xy);
        }

        public static void triangle_o12(ref double[] w, ref double[] xy)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_O12 returns a 12 point quadrature rule for the unit triangle.
            //
            //  Discussion:
            //
            //    This rule is precise for monomials through degree 6.
            //
            //    The integration region is:
            //
            //      0 <= X
            //      0 <= Y
            //      X + Y <= 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 April 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carlos Felippa,
            //    A compendium of FEM integration formulas for symbolic work,
            //    Engineering Computation,
            //    Volume 21, Number 8, 2004, pages 867-890.
            //
            //  Parameters:
            //
            //    Output, double W[12], the weights.
            //
            //    Output, double XY[2*12], the abscissas.
            //
        {
            int order = 12;

            double[] w_save =
            {
                0.050844906370206816921,
                0.050844906370206816921,
                0.050844906370206816921,
                0.11678627572637936603,
                0.11678627572637936603,
                0.11678627572637936603,
                0.082851075618373575194,
                0.082851075618373575194,
                0.082851075618373575194,
                0.082851075618373575194,
                0.082851075618373575194,
                0.082851075618373575194
            };
            double[] xy_save =
            {
                0.87382197101699554332, 0.063089014491502228340,
                0.063089014491502228340, 0.87382197101699554332,
                0.063089014491502228340, 0.063089014491502228340,
                0.50142650965817915742, 0.24928674517091042129,
                0.24928674517091042129, 0.50142650965817915742,
                0.24928674517091042129, 0.24928674517091042129,
                0.053145049844816947353, 0.31035245103378440542,
                0.31035245103378440542, 0.053145049844816947353,
                0.053145049844816947353, 0.63650249912139864723,
                0.31035245103378440542, 0.63650249912139864723,
                0.63650249912139864723, 0.053145049844816947353,
                0.63650249912139864723, 0.31035245103378440542
            };

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(2 * order, xy_save, ref xy);

        }


        public static double wedge_volume ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    WEDGE_VOLUME: volume of a unit wedge.
            //
            //  Discussion:
            //
            //    The integration region is:
            //
            //      0 <= X
            //      0 <= Y
            //      X + Y <= 1
            //      -1 <= Z <= 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 April 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double WEDGE_VOLUME, the volume.
            //
        {
            double value = 1.0;

            return value;
        }
    }
}