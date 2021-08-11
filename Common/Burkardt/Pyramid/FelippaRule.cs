using Burkardt.Types;

namespace Burkardt.Pyramid
{
    public static class FelippaRule
    {
        public static double pyramid_unit_monomial(int[] expon)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRAMID__UNIT_MONOMIAL: monomial integral in a unit pyramid.
            //
            //  Discussion:
            //
            //    This function returns the value of the integral of X^ALPHA Y^BETA Z^GAMMA
            //    over the unit pyramid.
            //
            //    The integration region is:
            //
            //    - ( 1 - Z ) <= X <= 1 - Z
            //    - ( 1 - Z ) <= Y <= 1 - Z
            //              0 <= Z <= 1.
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
            //    Output, double PYRAMID__UNIT_MONOMIAL, the volume of the pyramid.
            //
        {
            int i;
            int i_hi;
            double value;

            value = 0.0;

            if ((expon[0] % 2) == 0 && (expon[1] % 2) == 0)
            {
                i_hi = 2 + expon[0] + expon[1];

                for (i = 0; i <= i_hi; i++)
                {
                    value = value + typeMethods.r8_mop(i) * typeMethods.r8_choose(i_hi, i)
                        / (double)(i + expon[2] + 1);
                }

                value = value
                    * 2.0 / (double)(expon[0] + 1)
                    * 2.0 / (double)(expon[1] + 1);
            }

            return value;
        }

        public static void pyramid_unit_o01(ref double[] w, ref double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRAMID__UNIT_O01 returns a 1 point quadrature rule for the unit pyramid.
            //
            //  Discussion:
            //
            //    The integration region is:
            //
            //     - ( 1 - Z ) <= X <= 1 - Z
            //     - ( 1 - Z ) <= Y <= 1 - Z
            //               0 <= Z <= 1.
            //
            //    When Z is zero, the integration region is a square lying in the (X,Y) 
            //    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
            //    radius of the square diminishes, and when Z reaches 1, the square has 
            //    contracted to the single point (0,0,1).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 April 2009
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
            //    Output, double XYZ[3*1], the abscissas.
            //
        {
            int order = 1;

            double[] w_save = { 1.0 };

            double[] xyz_save = { 0.00, 0.00, 0.25 };

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(3 * order, xyz_save, ref xyz);

        }

        public static void pyramid_unit_o05(ref double[] w, ref double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRAMID__UNIT_O05 returns a 5 point quadrature rule for the unit pyramid.
            //
            //  Discussion:
            //
            //    The integration region is:
            //
            //     - ( 1 - Z ) <= X <= 1 - Z
            //     - ( 1 - Z ) <= Y <= 1 - Z
            //               0 <= Z <= 1.
            //
            //    When Z is zero, the integration region is a square lying in the (X,Y) 
            //    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
            //    radius of the square diminishes, and when Z reaches 1, the square has 
            //    contracted to the single point (0,0,1).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 April 2009
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
            //    Output, double W[5], the weights.
            //
            //    Output, double XYZ[3*5], the abscissas.
            //
        {
            int order = 5;

            double[] w_save =
            {
                0.21093750000000000000,
                0.21093750000000000000,
                0.21093750000000000000,
                0.21093750000000000000,
                0.15625000000000000000
            };

            double[] xyz_save =
            {
                -0.48686449556014765641, -0.48686449556014765641, 0.16666666666666666667,
                0.48686449556014765641, -0.48686449556014765641, 0.16666666666666666667,
                0.48686449556014765641, 0.48686449556014765641, 0.16666666666666666667,
                -0.48686449556014765641, 0.48686449556014765641, 0.16666666666666666667,
                0.00000000000000000000, 0.00000000000000000000, 0.70000000000000000000
            };

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(3 * order, xyz_save, ref xyz);

        }

        public static void pyramid_unit_o06(ref double[] w, ref double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRAMID__UNIT_O06 returns a 6 point quadrature rule for the unit pyramid.
            //
            //  Discussion:
            //
            //    The integration region is:
            //
            //     - ( 1 - Z ) <= X <= 1 - Z
            //     - ( 1 - Z ) <= Y <= 1 - Z
            //               0 <= Z <= 1.
            //
            //    When Z is zero, the integration region is a square lying in the (X,Y) 
            //    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
            //    radius of the square diminishes, and when Z reaches 1, the square has 
            //    contracted to the single point (0,0,1).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 April 2009
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
            //    Output, double XYZ[3*6], the abscissas.
            //
        {
            int order = 6;

            double[] w_save =
            {
                0.21000000000000000000,
                0.21000000000000000000,
                0.21000000000000000000,
                0.21000000000000000000,
                0.06000000000000000000,
                0.10000000000000000000
            };

            double[] xyz_save =
            {
                -0.48795003647426658968, -0.48795003647426658968, 0.16666666666666666667,
                0.48795003647426658968, -0.48795003647426658968, 0.16666666666666666667,
                0.48795003647426658968, 0.48795003647426658968, 0.16666666666666666667,
                -0.48795003647426658968, 0.48795003647426658968, 0.16666666666666666667,
                0.00000000000000000000, 0.00000000000000000000, 0.58333333333333333333,
                0.00000000000000000000, 0.00000000000000000000, 0.75000000000000000000
            };

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(3 * order, xyz_save, ref xyz);

        }

        public static void pyramid_unit_o08(ref double[] w, ref double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRAMID__UNIT_O08 returns an 8 point quadrature rule for the unit pyramid.
            //
            //  Discussion:
            //
            //    The integration region is:
            //
            //     - ( 1 - Z ) <= X <= 1 - Z
            //     - ( 1 - Z ) <= Y <= 1 - Z
            //               0 <= Z <= 1.
            //
            //    When Z is zero, the integration region is a square lying in the (X,Y) 
            //    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
            //    radius of the square diminishes, and when Z reaches 1, the square has 
            //    contracted to the single point (0,0,1).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 April 2009
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
            //    Output, double W[8], the weights.
            //
            //    Output, double XYZ[3*8], the abscissas.
            //
        {
            int order = 8;

            double[] w_save =
            {
                0.075589411559869072938,
                0.075589411559869072938,
                0.075589411559869072938,
                0.075589411559869072938,
                0.17441058844013092706,
                0.17441058844013092706,
                0.17441058844013092706,
                0.17441058844013092706
            };

            double[] xyz_save =
            {
                -0.26318405556971359557, -0.26318405556971359557, 0.54415184401122528880,
                0.26318405556971359557, -0.26318405556971359557, 0.54415184401122528880,
                0.26318405556971359557, 0.26318405556971359557, 0.54415184401122528880,
                -0.26318405556971359557, 0.26318405556971359557, 0.54415184401122528880,
                -0.50661630334978742377, -0.50661630334978742377, 0.12251482265544137787,
                0.50661630334978742377, -0.50661630334978742377, 0.12251482265544137787,
                0.50661630334978742377, 0.50661630334978742377, 0.12251482265544137787,
                -0.50661630334978742377, 0.50661630334978742377, 0.12251482265544137787
            };

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(3 * order, xyz_save, ref xyz);

        }

        public static void pyramid_unit_o08b(ref double[] w, ref double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRAMID__UNIT_O08B returns an 8 point quadrature rule for the unit pyramid.
            //
            //  Discussion:
            //
            //    The integration region is:
            //
            //     - ( 1 - Z ) <= X <= 1 - Z
            //     - ( 1 - Z ) <= Y <= 1 - Z
            //               0 <= Z <= 1.
            //
            //    When Z is zero, the integration region is a square lying in the (X,Y) 
            //    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
            //    radius of the square diminishes, and when Z reaches 1, the square has 
            //    contracted to the single point (0,0,1).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 April 2009
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
            //    Output, double W[8], the weights.
            //
            //    Output, double XYZ[3*8], the abscissas.
            //
        {
            int order = 8;

            double[] w_save =
            {
                0.16438287736328777572,
                0.16438287736328777572,
                0.16438287736328777572,
                0.16438287736328777572,
                0.085617122636712224276,
                0.085617122636712224276,
                0.085617122636712224276,
                0.085617122636712224276
            };

            double[] xyz_save =
            {
                -0.51197009372656270107, -0.51197009372656270107, 0.11024490204163285720,
                0.51197009372656270107, -0.51197009372656270107, 0.11024490204163285720,
                0.51197009372656270107, 0.51197009372656270107, 0.11024490204163285720,
                -0.51197009372656270107, 0.51197009372656270107, 0.11024490204163285720,
                -0.28415447557052037456, -0.28415447557052037456, 0.518326526529795714229,
                0.28415447557052037456, -0.28415447557052037456, 0.518326526529795714229,
                0.28415447557052037456, 0.28415447557052037456, 0.518326526529795714229,
                -0.28415447557052037456, 0.28415447557052037456, 0.518326526529795714229
            };

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(3 * order, xyz_save, ref xyz);

        }

        public static void pyramid_unit_o09(ref double[] w, ref double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRAMID__UNIT_O09 returns a 9 point quadrature rule for the unit pyramid.
            ///
            //  Discussion:
            //
            //    The integration region is:
            //
            //     - ( 1 - Z ) <= X <= 1 - Z
            //     - ( 1 - Z ) <= Y <= 1 - Z
            //               0 <= Z <= 1.
            //
            //    When Z is zero, the integration region is a square lying in the (X,Y) 
            //    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
            //    radius of the square diminishes, and when Z reaches 1, the square has 
            //    contracted to the single point (0,0,1).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 April 2009
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
            //    Output, double W[9], the weights.
            //
            //    Output, double XYZ[3*9], the abscissas.
            //
        {
            int order = 9;

            double[] w_save =
            {
                0.13073389672275944791,
                0.13073389672275944791,
                0.13073389672275944791,
                0.13073389672275944791,
                0.10989110327724055209,
                0.10989110327724055209,
                0.10989110327724055209,
                0.10989110327724055209,
                0.03750000000000000000
            };

            double[] xyz_save =
            {
                -0.52966422253852215131, -0.52966422253852215131, 0.08176876558246862335,
                0.52966422253852215131, -0.52966422253852215131, 0.08176876558246862335,
                0.52966422253852215131, 0.52966422253852215131, 0.08176876558246862335,
                -0.52966422253852215131, 0.52966422253852215131, 0.08176876558246862335,
                -0.34819753825720418039, -0.34819753825720418039, 0.400374091560388519511,
                0.34819753825720418039, -0.34819753825720418039, 0.400374091560388519511,
                0.34819753825720418039, 0.34819753825720418039, 0.400374091560388519511,
                -0.34819753825720418039, 0.34819753825720418039, 0.400374091560388519511,
                0.00000000000000000000, 0.00000000000000000000, 0.83333333333333333333
            };

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(3 * order, xyz_save, ref xyz);
        }

        public static void pyramid_unit_o13(ref double[] w, ref double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRAMID__UNIT_O13 returns a 13 point quadrature rule for the unit pyramid.
            //
            //  Discussion:
            //
            //    The integration region is:
            //
            //     - ( 1 - Z ) <= X <= 1 - Z
            //     - ( 1 - Z ) <= Y <= 1 - Z
            //               0 <= Z <= 1.
            //
            //    When Z is zero, the integration region is a square lying in the (X,Y) 
            //    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
            //    radius of the square diminishes, and when Z reaches 1, the square has 
            //    contracted to the single point (0,0,1).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 April 2009
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
            //    Output, double W[13], the weights.
            //
            //    Output, double XYZ[3*13], the abscissas.
            //
        {
            int order = 13;

            double[] w_save =
            {
                0.063061594202898550725,
                0.063061594202898550725,
                0.063061594202898550725,
                0.063061594202898550725,
                0.042101946815575556199,
                0.042101946815575556199,
                0.042101946815575556199,
                0.042101946815575556199,
                0.13172030707666776585,
                0.13172030707666776585,
                0.13172030707666776585,
                0.13172030707666776585,
                0.05246460761943250889
            };

            double[] xyz_save =
            {
                -0.38510399211870384331, -0.38510399211870384331, 0.428571428571428571429,
                0.38510399211870384331, -0.38510399211870384331, 0.428571428571428571429,
                0.38510399211870384331, 0.38510399211870384331, 0.428571428571428571429,
                -0.38510399211870384331, 0.38510399211870384331, 0.428571428571428571429,
                -0.40345831960728204766, 0.00000000000000000000, 0.33928571428571428571,
                0.40345831960728204766, 0.00000000000000000000, 0.33928571428571428571,
                0.00000000000000000000, -0.40345831960728204766, 0.33928571428571428571,
                0.00000000000000000000, 0.40345831960728204766, 0.33928571428571428571,
                -0.53157877436961973359, -0.53157877436961973359, 0.08496732026143790850,
                0.53157877436961973359, -0.53157877436961973359, 0.08496732026143790850,
                0.53157877436961973359, 0.53157877436961973359, 0.08496732026143790850,
                -0.53157877436961973359, 0.53157877436961973359, 0.08496732026143790850,
                0.00000000000000000000, 0.00000000000000000000, 0.76219701803768503595
            };

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(3 * order, xyz_save, ref xyz);
        }

        public static void pyramid_unit_o18(ref double[] w, ref double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRAMID__UNIT_O18 returns an 18 point quadrature rule for the unit pyramid.
            //
            //  Discussion:
            //
            //    The integration region is:
            //
            //     - ( 1 - Z ) <= X <= 1 - Z
            //     - ( 1 - Z ) <= Y <= 1 - Z
            //               0 <= Z <= 1.
            //
            //    When Z is zero, the integration region is a square lying in the (X,Y) 
            //    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
            //    radius of the square diminishes, and when Z reaches 1, the square has 
            //    contracted to the single point (0,0,1).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 April 2009
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
            //    Output, double W[18], the weights.
            //
            //    Output, double XYZ[3*18], the abscissas.
            //
        {
            int order = 18;

            double[] w_save =
            {
                0.023330065296255886709,
                0.037328104474009418735,
                0.023330065296255886709,
                0.037328104474009418735,
                0.059724967158415069975,
                0.037328104474009418735,
                0.023330065296255886709,
                0.037328104474009418735,
                0.023330065296255886709,
                0.05383042853090460712,
                0.08612868564944737139,
                0.05383042853090460712,
                0.08612868564944737139,
                0.13780589703911579422,
                0.08612868564944737139,
                0.05383042853090460712,
                0.08612868564944737139,
                0.05383042853090460712
            };

            double[] xyz_save =
            {
                -0.35309846330877704481, -0.35309846330877704481, 0.544151844011225288800,
                0.00000000000000000000, -0.35309846330877704481, 0.544151844011225288800,
                0.35309846330877704481, -0.35309846330877704481, 0.544151844011225288800,
                -0.35309846330877704481, 0.00000000000000000000, 0.544151844011225288800,
                0.00000000000000000000, 0.00000000000000000000, 0.544151844011225288800,
                0.35309846330877704481, 0.00000000000000000000, 0.544151844011225288800,
                -0.35309846330877704481, 0.35309846330877704481, 0.544151844011225288800,
                0.00000000000000000000, 0.35309846330877704481, 0.544151844011225288800,
                0.35309846330877704481, 0.35309846330877704481, 0.544151844011225288800,
                -0.67969709567986745790, -0.67969709567986745790, 0.12251482265544137787,
                0.00000000000000000000, -0.67969709567986745790, 0.12251482265544137787,
                0.67969709567986745790, -0.67969709567986745790, 0.12251482265544137787,
                -0.67969709567986745790, 0.00000000000000000000, 0.12251482265544137787,
                0.00000000000000000000, 0.00000000000000000000, 0.12251482265544137787,
                0.67969709567986745790, 0.00000000000000000000, 0.12251482265544137787,
                -0.67969709567986745790, 0.67969709567986745790, 0.12251482265544137787,
                0.00000000000000000000, 0.67969709567986745790, 0.12251482265544137787,
                0.67969709567986745790, 0.67969709567986745790, 0.12251482265544137787
            };

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(3 * order, xyz_save, ref xyz);

            return;
        }

        public static void pyramid_unit_o27(ref double[] w, ref double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRAMID__UNIT_O27 returns a 27 point quadrature rule for the unit pyramid.
            //
            //  Discussion:
            //
            //    The integration region is:
            //
            //     - ( 1 - Z ) <= X <= 1 - Z
            //     - ( 1 - Z ) <= Y <= 1 - Z
            //               0 <= Z <= 1.
            //
            //    When Z is zero, the integration region is a square lying in the (X,Y) 
            //    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
            //    radius of the square diminishes, and when Z reaches 1, the square has 
            //    contracted to the single point (0,0,1).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 April 2009
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
            //    Output, double W[27], the weights.
            //
            //    Output, double XYZ[3*27], the abscissas.
            //
        {
            int order = 27;

            double[] w_save =
            {
                0.036374157653908938268,
                0.05819865224625430123,
                0.036374157653908938268,
                0.05819865224625430123,
                0.09311784359400688197,
                0.05819865224625430123,
                0.036374157653908938268,
                0.05819865224625430123,
                0.036374157653908938268,
                0.033853303069413431019,
                0.054165284911061489631,
                0.033853303069413431019,
                0.054165284911061489631,
                0.08666445585769838341,
                0.054165284911061489631,
                0.033853303069413431019,
                0.054165284911061489631,
                0.033853303069413431019,
                0.006933033103838124540,
                0.011092852966140999264,
                0.006933033103838124540,
                0.011092852966140999264,
                0.017748564745825598822,
                0.011092852966140999264,
                0.006933033103838124540,
                0.011092852966140999264,
                0.006933033103838124540
            };

            double[] xyz_save =
            {
                -0.7180557413198889387, -0.7180557413198889387, 0.07299402407314973216,
                0.00000000000000000000, -0.7180557413198889387, 0.07299402407314973216,
                0.7180557413198889387, -0.7180557413198889387, 0.07299402407314973216,
                -0.7180557413198889387, 0.00000000000000000000, 0.07299402407314973216,
                0.00000000000000000000, 0.00000000000000000000, 0.07299402407314973216,
                0.7180557413198889387, 0.00000000000000000000, 0.07299402407314973216,
                -0.7180557413198889387, 0.7180557413198889387, 0.07299402407314973216,
                0.00000000000000000000, 0.7180557413198889387, 0.07299402407314973216,
                0.7180557413198889387, 0.7180557413198889387, 0.07299402407314973216,
                -0.50580870785392503961, -0.50580870785392503961, 0.34700376603835188472,
                0.00000000000000000000, -0.50580870785392503961, 0.34700376603835188472,
                0.50580870785392503961, -0.50580870785392503961, 0.34700376603835188472,
                -0.50580870785392503961, 0.00000000000000000000, 0.34700376603835188472,
                0.00000000000000000000, 0.00000000000000000000, 0.34700376603835188472,
                0.50580870785392503961, 0.00000000000000000000, 0.34700376603835188472,
                -0.50580870785392503961, 0.50580870785392503961, 0.34700376603835188472,
                0.00000000000000000000, 0.50580870785392503961, 0.34700376603835188472,
                0.50580870785392503961, 0.50580870785392503961, 0.34700376603835188472,
                -0.22850430565396735360, -0.22850430565396735360, 0.70500220988849838312,
                0.00000000000000000000, -0.22850430565396735360, 0.70500220988849838312,
                0.22850430565396735360, -0.22850430565396735360, 0.70500220988849838312,
                -0.22850430565396735360, 0.00000000000000000000, 0.70500220988849838312,
                0.00000000000000000000, 0.00000000000000000000, 0.70500220988849838312,
                0.22850430565396735360, 0.00000000000000000000, 0.70500220988849838312,
                -0.22850430565396735360, 0.22850430565396735360, 0.70500220988849838312,
                0.00000000000000000000, 0.22850430565396735360, 0.70500220988849838312,
                0.22850430565396735360, 0.22850430565396735360, 0.70500220988849838312
            };

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(3 * order, xyz_save, ref xyz);

        }

        public static void pyramid_unit_o48(ref double[] w, ref double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRAMID__UNIT_O48 returns a 48 point quadrature rule for the unit pyramid.
            //
            //  Discussion:
            //
            //    The integration region is:
            //
            //     - ( 1 - Z ) <= X <= 1 - Z
            //     - ( 1 - Z ) <= Y <= 1 - Z
            //               0 <= Z <= 1.
            //
            //    When Z is zero, the integration region is a square lying in the (X,Y) 
            //    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
            //    radius of the square diminishes, and when Z reaches 1, the square has 
            //    contracted to the single point (0,0,1).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 April 2009
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
            //    Output, double W[48], the weights.
            //
            //    Output, double XYZ[3*48], the abscissas.
            //
        {
            int order = 48;

            double[] w_save =
            {
                2.01241939442682455E-002,
                2.01241939442682455E-002,
                2.01241939442682455E-002,
                2.01241939442682455E-002,
                2.60351137043010779E-002,
                2.60351137043010779E-002,
                2.60351137043010779E-002,
                2.60351137043010779E-002,
                1.24557795239745531E-002,
                1.24557795239745531E-002,
                1.24557795239745531E-002,
                1.24557795239745531E-002,
                1.87873998794808156E-003,
                1.87873998794808156E-003,
                1.87873998794808156E-003,
                1.87873998794808156E-003,
                4.32957927807745280E-002,
                4.32957927807745280E-002,
                4.32957927807745280E-002,
                4.32957927807745280E-002,
                1.97463249834127288E-002,
                1.97463249834127288E-002,
                1.97463249834127288E-002,
                1.97463249834127288E-002,
                5.60127223523590526E-002,
                5.60127223523590526E-002,
                5.60127223523590526E-002,
                5.60127223523590526E-002,
                2.55462562927473852E-002,
                2.55462562927473852E-002,
                2.55462562927473852E-002,
                2.55462562927473852E-002,
                2.67977366291788643E-002,
                2.67977366291788643E-002,
                2.67977366291788643E-002,
                2.67977366291788643E-002,
                1.22218992265373354E-002,
                1.22218992265373354E-002,
                1.22218992265373354E-002,
                1.22218992265373354E-002,
                4.04197740453215038E-003,
                4.04197740453215038E-003,
                4.04197740453215038E-003,
                4.04197740453215038E-003,
                1.84346316995826843E-003,
                1.84346316995826843E-003,
                1.84346316995826843E-003,
                1.84346316995826843E-003
            };

            double[] xyz_save =
            {
                0.88091731624450909E+00, 0.00000000000000000E+00, 4.85005494469969989E-02,
                -0.88091731624450909E+00, 0.00000000000000000E+00, 4.85005494469969989E-02,
                0.00000000000000000E+00, 0.88091731624450909E+00, 4.85005494469969989E-02,
                0.00000000000000000E+00, -0.88091731624450909E+00, 4.85005494469969989E-02,
                0.70491874112648223E+00, 0.00000000000000000E+00, 0.23860073755186201E+00,
                -0.70491874112648223E+00, 0.00000000000000000E+00, 0.23860073755186201E+00,
                0.00000000000000000E+00, 0.70491874112648223E+00, 0.23860073755186201E+00,
                0.00000000000000000E+00, -0.70491874112648223E+00, 0.23860073755186201E+00,
                0.44712732143189760E+00, 0.00000000000000000E+00, 0.51704729510436798E+00,
                -0.44712732143189760E+00, 0.00000000000000000E+00, 0.51704729510436798E+00,
                0.00000000000000000E+00, 0.44712732143189760E+00, 0.51704729510436798E+00,
                0.00000000000000000E+00, -0.44712732143189760E+00, 0.51704729510436798E+00,
                0.18900486065123448E+00, 0.00000000000000000E+00, 0.79585141789677305E+00,
                -0.18900486065123448E+00, 0.00000000000000000E+00, 0.79585141789677305E+00,
                0.00000000000000000E+00, 0.18900486065123448E+00, 0.79585141789677305E+00,
                0.00000000000000000E+00, -0.18900486065123448E+00, 0.79585141789677305E+00,
                0.36209733410322176E+00, 0.36209733410322176E+00, 4.85005494469969989E-02,
                -0.36209733410322176E+00, 0.36209733410322176E+00, 4.85005494469969989E-02,
                -0.36209733410322176E+00, -0.36209733410322176E+00, 4.85005494469969989E-02,
                0.36209733410322176E+00, -0.36209733410322176E+00, 4.85005494469969989E-02,
                0.76688932060387538E+00, 0.76688932060387538E+00, 4.85005494469969989E-02,
                -0.76688932060387538E+00, 0.76688932060387538E+00, 4.85005494469969989E-02,
                -0.76688932060387538E+00, -0.76688932060387538E+00, 4.85005494469969989E-02,
                0.76688932060387538E+00, -0.76688932060387538E+00, 4.85005494469969989E-02,
                0.28975386476618070E+00, 0.28975386476618070E+00, 0.23860073755186201E+00,
                -0.28975386476618070E+00, 0.28975386476618070E+00, 0.23860073755186201E+00,
                -0.28975386476618070E+00, -0.28975386476618070E+00, 0.23860073755186201E+00,
                0.28975386476618070E+00, -0.28975386476618070E+00, 0.23860073755186201E+00,
                0.61367241226233160E+00, 0.61367241226233160E+00, 0.23860073755186201E+00,
                -0.61367241226233160E+00, 0.61367241226233160E+00, 0.23860073755186201E+00,
                -0.61367241226233160E+00, -0.61367241226233160E+00, 0.23860073755186201E+00,
                0.61367241226233160E+00, -0.61367241226233160E+00, 0.23860073755186201E+00,
                0.18378979287798017E+00, 0.18378979287798017E+00, 0.51704729510436798E+00,
                -0.18378979287798017E+00, 0.18378979287798017E+00, 0.51704729510436798E+00,
                -0.18378979287798017E+00, -0.18378979287798017E+00, 0.51704729510436798E+00,
                0.18378979287798017E+00, -0.18378979287798017E+00, 0.51704729510436798E+00,
                0.38925011625173161E+00, 0.38925011625173161E+00, 0.51704729510436798E+00,
                -0.38925011625173161E+00, 0.38925011625173161E+00, 0.51704729510436798E+00,
                -0.38925011625173161E+00, -0.38925011625173161E+00, 0.51704729510436798E+00,
                0.38925011625173161E+00, -0.38925011625173161E+00, 0.51704729510436798E+00,
                7.76896479525748113E-02, 7.76896479525748113E-02, 0.79585141789677305E+00,
                -7.76896479525748113E-02, 7.76896479525748113E-02, 0.79585141789677305E+00,
                -7.76896479525748113E-02, -7.76896479525748113E-02, 0.79585141789677305E+00,
                7.76896479525748113E-02, -7.76896479525748113E-02, 0.79585141789677305E+00,
                0.16453962988669860E+00, 0.16453962988669860E+00, 0.79585141789677305E+00,
                -0.16453962988669860E+00, 0.16453962988669860E+00, 0.79585141789677305E+00,
                -0.16453962988669860E+00, -0.16453962988669860E+00, 0.79585141789677305E+00,
                0.16453962988669860E+00, -0.16453962988669860E+00, 0.79585141789677305E+00
            };

            typeMethods.r8vec_copy(order, w_save, ref w);
            typeMethods.r8vec_copy(3 * order, xyz_save, ref xyz);

        }

        public static double pyramid_unit_volume()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRAMID__UNIT_VOLUME: volume of a unit pyramid with square base.
            //
            //  Discussion:
            //
            //    The volume of this unit pyramid is 4/3.
            //
            //    The integration region is:
            //
            //      - ( 1 - Z ) <= X <= 1 - Z
            //      - ( 1 - Z ) <= Y <= 1 - Z
            //                0 <= Z <= 1.
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
            //    Output, double PYRAMID__UNIT_VOLUME, the volume of the pyramid.
            //
        {
            double volume;

            volume = 4.0 / 3.0;

            return volume;
        }
    }
}