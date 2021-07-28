using Burkardt.Types;

namespace Burkardt.TriangleNS
{
    public static class QuadratureRule
    {
        public static double triangle_unit_monomial(int[] expon )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_MONOMIAL integrates a monomial over the unit triangle.
        //
        //  Discussion:
        //
        //    This routine integrates a monomial of the form
        //
        //      product ( 1 <= dim <= 2 ) x(dim)^expon(dim)
        //
        //    where the exponents are nonnegative integers.  Note that
        //    if the combination 0^0 is encountered, it should be treated
        //    as 1.
        //
        //    Integral ( over unit triangle ) x^m y^n dx dy = m! * n! / ( m + n + 2 )!
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
        //    18 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int EXPON[2], the exponents.
        //
        //    Output, double TRIANGLE_UNIT_MONOMIAL, the integral of the monomial.
        //
        {
            int i;
            int k;
            double value;
            //
            //  The first computation ends with VALUE = 1.0;
            //
            value = 1.0;

            // k = 0;
            //
            //  The first loop simply computes 1 so we short circuit it!
            //
            // for ( i = 1; i <= expon[0]; i++ )
            // {
            //   k = k + 1;
            //   value = value * ( double ) ( i ) / ( double ) ( k );
            // }

            k = expon[0];

            for (i = 1; i <= expon[1]; i++)
            {
                k = k + 1;
                value = value * (double) (i) / (double) (k);
            }

            k = k + 1;
            value = value / (double) (k);

            k = k + 1;
            value = value / (double) (k);

            return value;
        }

        public static void triangle_unit_o01(ref double[] w, ref double[] xy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_O01 returns a 1 point quadrature rule for the unit triangle.
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

        public static void triangle_unit_o03(ref double[] w, ref double[] xy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_O03 returns a 3 point quadrature rule for the unit triangle.
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

        public static void triangle_unit_o03b(ref double[] w, ref double[] xy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_O03B returns a 3 point quadrature rule for the unit triangle.
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

        public static void triangle_unit_o06(ref double[] w, ref double[] xy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_O06 returns a 6 point quadrature rule for the unit triangle.
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

            double[] w_save =  {
                0.22338158967801146570,
                0.22338158967801146570,
                0.22338158967801146570,
                0.10995174365532186764,
                0.10995174365532186764,
                0.10995174365532186764
            }
            ;
            double[] xy_save =  {
                0.10810301816807022736, 0.44594849091596488632,
                0.44594849091596488632, 0.10810301816807022736,
                0.44594849091596488632, 0.44594849091596488632,
                0.81684757298045851308, 0.091576213509770743460,
                0.091576213509770743460, 0.81684757298045851308,
                0.091576213509770743460, 0.091576213509770743460
            }
            ;

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(2 * order, xy_save, ref xy);
        }

        public static void triangle_unit_o06b(ref double[] w, ref double[] xy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_O06B returns a 6 point quadrature rule for the unit triangle.
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

            double[] w_save =  {
                0.30000000000000000000,
                0.30000000000000000000,
                0.30000000000000000000,
                0.033333333333333333333,
                0.033333333333333333333,
                0.033333333333333333333
            }
            ;
            double[] xy_save =  {
                0.66666666666666666667, 0.16666666666666666667,
                0.16666666666666666667, 0.66666666666666666667,
                0.16666666666666666667, 0.16666666666666666667,
                0.0, 0.5,
                0.5, 0.0,
                0.5, 0.5
            }
            ;

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(2 * order, xy_save, ref xy);
        }

        public static void triangle_unit_o07(ref double[] w, ref double[] xy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_O07 returns a 7 point quadrature rule for the unit triangle.
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

            double[] w_save =  {
                0.12593918054482715260,
                0.12593918054482715260,
                0.12593918054482715260,
                0.13239415278850618074,
                0.13239415278850618074,
                0.13239415278850618074,
                0.22500000000000000000
            }
            ;
            double[] xy_save =  {
                0.79742698535308732240, 0.10128650732345633880,
                0.10128650732345633880, 0.79742698535308732240,
                0.10128650732345633880, 0.10128650732345633880,
                0.059715871789769820459, 0.47014206410511508977,
                0.47014206410511508977, 0.059715871789769820459,
                0.47014206410511508977, 0.47014206410511508977,
                0.33333333333333333333, 0.33333333333333333333
            }
            ;

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(2 * order, xy_save, ref xy);
        }

        public static void triangle_unit_o12(ref double[] w, ref double[] xy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_O12 returns a 12 point quadrature rule for the unit triangle.
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

            double[] w_save =  {
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
            }
            ;
            double[] xy_save =  {
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
            }
            ;

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(2 * order, xy_save, ref xy);
        }

        public static double triangle_unit_volume()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_UNIT_VOLUME returns the "volume" of the unit triangle in 2D.
            //
            //  Discussion:
            //
            //    The "volume" of a triangle is usually called its area.
            //
            //    The integration region is:
            //
            //      0 <= X,
            //      0 <= Y, 
            //      X + Y <= 1.
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
            //    Output, double TRIANGLE_UNIT_VOLUME, the volume.
            //
        {
            double volume = 1.0 / 2.0;

            return volume;
        }
    }
}