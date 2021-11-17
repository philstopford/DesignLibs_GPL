using System;
using Burkardt.Types;

namespace Burkardt.Stroud;

public static class Pyramid
{
    public static double pyramid_unit_o01_3d(int setting, Func<int, double, double, double, double> func)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID_UNIT_O01_3D approximates an integral inside the unit pyramid in 3D.
        //
        //  Discussion:
        //
        //    A 1 point degree 1 formula is used.
        //
        //    The (X,Y,Z) integration region can be represented as:
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
        //    22 March 2008
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
        //    user supplied function which is to be integrated.
        //
        //    Output, double RESULT, the approximate integral of the function.
        //
    {
        double quad;
        double result;
        double volume;
        double w;
        double x;
        double y;
        double z;
        //
        //  Quadrature.
        //
        quad = 0.0;

        x = 0.0;
        y = 0.0;
        z = 1.0 / 4.0;
        w = 1.0;

        quad += w * func(setting, x, y, z);
        //
        //  Volume.
        //
        volume = pyramid_unit_volume_3d();
        //
        //  Result.
        //
        result = quad * volume;

        return result;
    }

    public static double pyramid_unit_o05_3d(int setting, Func<int, double, double, double, double> func)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID_UNIT_O05_3D approximates an integral inside the unit pyramid in 3D.
        //
        //  Discussion:
        //
        //    A 5 point formula is used.
        //
        //    The (X,Y,Z) integration region can be represented as:
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
        //    22 March 2008
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
        //    Input, Func< double, double, double, double > func, the name of the 
        //    user supplied function which is to be integrated.
        //
        //    Output, double PYRAMID_UNIT_O05_3D, the approximate integral
        //    of the function.
        //
    {
        int i;
        int order = 5;
        double quad;
        double result;
        double volume;
        double[] w =
        {
            0.21093750000000000000,
            0.21093750000000000000,
            0.21093750000000000000,
            0.21093750000000000000,
            0.15625000000000000000
        };
        double[] x =
        {
            -0.48686449556014765641,
            0.48686449556014765641,
            0.48686449556014765641,
            -0.48686449556014765641,
            0.00000000000000000000
        };
        double[] y =
        {
            -0.48686449556014765641,
            -0.48686449556014765641,
            0.48686449556014765641,
            0.48686449556014765641,
            0.00000000000000000000
        };
        double[] z =
        {
            0.16666666666666666667,
            0.16666666666666666667,
            0.16666666666666666667,
            0.16666666666666666667,
            0.70000000000000000000
        };
        //
        //  Quadrature.
        //
        quad = 0.0;
        for (i = 0; i < order; i++)
        {
            quad += w[i] * func(setting, x[i], y[i], z[i]);
        }

        //
        //  Volume.
        //
        volume = pyramid_unit_volume_3d();
        //
        //  Result.
        //
        result = quad * volume;

        return result;
    }

    public static double pyramid_unit_o06_3d(int setting, Func<int, double, double, double, double> func)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID_UNIT_O06_3D approximates an integral inside the unit pyramid in 3D.
        //
        //  Discussion:
        //
        //    A 6 point formula is used.
        //
        //    The (X,Y,Z) integration region can be represented as:
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
        //    22 March 2008
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
        //    Input, Func< double, double, double, double > func, the name of the 
        //    user supplied function which is to be integrated.
        //
        //    Output, double PYRAMID_UNIT_O06_3D, the approximate integral
        //    of the function.
        //
    {
        int i;
        int order = 6;
        double quad;
        double result;
        double volume;
        double[] w =
        {
            0.21000000000000000000,
            0.21000000000000000000,
            0.21000000000000000000,
            0.21000000000000000000,
            0.06000000000000000000,
            0.10000000000000000000
        };
        double[] x =
        {
            -0.48795003647426658968,
            0.48795003647426658968,
            0.48795003647426658968,
            -0.48795003647426658968,
            0.00000000000000000000,
            0.00000000000000000000
        };
        double[] y =
        {
            -0.48795003647426658968,
            -0.48795003647426658968,
            0.48795003647426658968,
            0.48795003647426658968,
            0.00000000000000000000,
            0.00000000000000000000
        };
        double[] z =
        {
            0.16666666666666666667,
            0.16666666666666666667,
            0.16666666666666666667,
            0.16666666666666666667,
            0.58333333333333333333,
            0.75000000000000000000
        };
        //
        //  Quadrature.
        //
        quad = 0.0;
        for (i = 0; i < order; i++)
        {
            quad += w[i] * func(setting, x[i], y[i], z[i]);
        }

        //
        //  Volume.
        //
        volume = pyramid_unit_volume_3d();
        //
        //  Result.
        //
        result = quad * volume;

        return result;
    }

    public static double pyramid_unit_o08_3d(int setting, Func<int, double, double, double, double> func)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID_UNIT_O08_3D approximates an integral inside the unit pyramid in 3D.
        //
        //  Discussion:
        //
        //    An 8 point formula is used.
        //
        //    The (X,Y,Z) integration region can be represented as:
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
        //    22 March 2008
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
        //    Input, Func< double, double, double, double > func, the name of the 
        //    user supplied function which is to be integrated.
        //
        //    Output, double PYRAMID_UNIT_O08_3D, the approximate integral
        //    of the function.
        //
    {
        int i;
        int order = 8;
        double quad;
        double result;
        double volume;
        double[] w =
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
        double[] x =
        {
            -0.26318405556971359557,
            0.26318405556971359557,
            0.26318405556971359557,
            -0.26318405556971359557,
            -0.50661630334978742377,
            0.50661630334978742377,
            0.50661630334978742377,
            -0.50661630334978742377
        };
        double[] y =
        {
            -0.26318405556971359557,
            -0.26318405556971359557,
            0.26318405556971359557,
            0.26318405556971359557,
            -0.50661630334978742377,
            -0.50661630334978742377,
            0.50661630334978742377,
            0.50661630334978742377
        };
        double[] z =
        {
            0.54415184401122528880,
            0.54415184401122528880,
            0.54415184401122528880,
            0.54415184401122528880,
            0.12251482265544137787,
            0.12251482265544137787,
            0.12251482265544137787,
            0.12251482265544137787
        };
        //
        //  Quadrature.
        //
        quad = 0.0;
        for (i = 0; i < order; i++)
        {
            quad += w[i] * func(setting, x[i], y[i], z[i]);
        }

        //
        //  Volume.
        //
        volume = pyramid_unit_volume_3d();
        //
        //  Result.
        //
        result = quad * volume;

        return result;
    }

    public static double pyramid_unit_o08b_3d(int setting, Func<int, double, double, double, double> func)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID_UNIT_O08B_3D approximates an integral inside the unit pyramid in 3D.
        //
        //  Discussion:
        //
        //    An 8 point formula is used.
        //
        //    The (X,Y,Z) integration region can be represented as:
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
        //    22 March 2008
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
        //    Input, Func< double, double, double, double > func, the name of the 
        //    user supplied function which is to be integrated.
        //
        //    Output, double PYRAMID_UNIT_O08B_3D, the approximate integral
        //    of the function.
        //
    {
        int i;
        int order = 8;
        double quad;
        double result;
        double volume;
        double[] w =
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
        double[] x =
        {
            -0.51197009372656270107,
            0.51197009372656270107,
            0.51197009372656270107,
            -0.51197009372656270107,
            -0.28415447557052037456,
            0.28415447557052037456,
            0.28415447557052037456,
            -0.28415447557052037456
        };
        double[] y =
        {
            -0.51197009372656270107,
            -0.51197009372656270107,
            0.51197009372656270107,
            0.51197009372656270107,
            -0.28415447557052037456,
            -0.28415447557052037456,
            0.28415447557052037456,
            0.28415447557052037456
        };
        double[] z =
        {
            0.11024490204163285720,
            0.11024490204163285720,
            0.11024490204163285720,
            0.11024490204163285720,
            0.518326526529795714229,
            0.518326526529795714229,
            0.518326526529795714229,
            0.518326526529795714229
        };
        //
        //  Quadrature.
        //
        quad = 0.0;
        for (i = 0; i < order; i++)
        {
            quad += w[i] * func(setting, x[i], y[i], z[i]);
        }

        //
        //  Volume.
        //
        volume = pyramid_unit_volume_3d();
        //
        //  Result.
        //
        result = quad * volume;

        return result;
    }

    public static double pyramid_unit_o09_3d(int setting, Func<int, double, double, double, double> func)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID_UNIT_O09_3D approximates an integral inside the unit pyramid in 3D.
        //
        //  Discussion:
        //
        //    A 9 point formula is used.
        //
        //    The (X,Y,Z) integration region can be represented as:
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
        //    22 March 2008
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
        //    Input, Func< double, double, double, double > func, the name of the 
        //    user supplied function which is to be integrated.
        //
        //    Output, double PYRAMID_UNIT_O09_3D, the approximate integral
        //    of the function.
        //
    {
        int i;
        int order = 9;
        double quad;
        double result;
        double volume;
        double[] w =
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
        double[] x =
        {
            -0.52966422253852215131,
            0.52966422253852215131,
            0.52966422253852215131,
            -0.52966422253852215131,
            -0.34819753825720418039,
            0.34819753825720418039,
            0.34819753825720418039,
            -0.34819753825720418039,
            0.00000000000000000000
        };
        double[] y =
        {
            -0.52966422253852215131,
            -0.52966422253852215131,
            0.52966422253852215131,
            0.52966422253852215131,
            -0.34819753825720418039,
            -0.34819753825720418039,
            0.34819753825720418039,
            0.34819753825720418039,
            0.00000000000000000000
        };
        double[] z =
        {
            0.08176876558246862335,
            0.08176876558246862335,
            0.08176876558246862335,
            0.08176876558246862335,
            0.400374091560388519511,
            0.400374091560388519511,
            0.400374091560388519511,
            0.400374091560388519511,
            0.83333333333333333333
        };
        //
        //  Quadrature.
        //
        quad = 0.0;
        for (i = 0; i < order; i++)
        {
            quad += w[i] * func(setting, x[i], y[i], z[i]);
        }

        //
        //  Volume.
        //
        volume = pyramid_unit_volume_3d();
        //
        //  Result.
        //
        result = quad * volume;

        return result;
    }

    public static double pyramid_unit_o13_3d(int setting, Func<int, double, double, double, double> func)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID_UNIT_O13_3D approximates an integral inside the unit pyramid in 3D.
        //
        //  Discussion:
        //
        //    A 13 point formula is used.
        //
        //    The (X,Y,Z) integration region can be represented as:
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
        //    22 March 2008
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
        //    Input, Func< double, double, double, double > func, the name of the 
        //    user supplied function which is to be integrated.
        //
        //    Output, double PYRAMID_UNIT_O13_3D, the approximate integral
        //    of the function.
        //
    {
        int i;
        int order = 13;
        double quad;
        double result;
        double volume;
        double[] w =
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
        double[] x =
        {
            -0.38510399211870384331,
            0.38510399211870384331,
            0.38510399211870384331,
            -0.38510399211870384331,
            -0.40345831960728204766,
            0.40345831960728204766,
            0.00000000000000000000,
            0.00000000000000000000,
            -0.53157877436961973359,
            0.53157877436961973359,
            0.53157877436961973359,
            -0.53157877436961973359,
            0.00000000000000000000
        };
        double[] y =
        {
            -0.38510399211870384331,
            -0.38510399211870384331,
            0.38510399211870384331,
            0.38510399211870384331,
            0.00000000000000000000,
            0.00000000000000000000,
            -0.40345831960728204766,
            0.40345831960728204766,
            -0.53157877436961973359,
            -0.53157877436961973359,
            0.53157877436961973359,
            0.53157877436961973359,
            0.00000000000000000000
        };
        double[] z =
        {
            0.428571428571428571429,
            0.428571428571428571429,
            0.428571428571428571429,
            0.428571428571428571429,
            0.33928571428571428571,
            0.33928571428571428571,
            0.33928571428571428571,
            0.33928571428571428571,
            0.08496732026143790850,
            0.08496732026143790850,
            0.08496732026143790850,
            0.08496732026143790850,
            0.76219701803768503595
        };
        //
        //  Quadrature.
        //
        quad = 0.0;
        for (i = 0; i < order; i++)
        {
            quad += w[i] * func(setting, x[i], y[i], z[i]);
        }

        //
        //  Volume.
        //
        volume = pyramid_unit_volume_3d();
        //
        //  Result.
        //
        result = quad * volume;

        return result;
    }

    public static double pyramid_unit_o18_3d(int setting, Func<int, double, double, double, double> func)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID_UNIT_O18_3D approximates an integral inside the unit pyramid in 3D.
        //
        //  Discussion:
        //
        //    An 18 point formula is used.
        //
        //    The (X,Y,Z) integration region can be represented as:
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
        //    22 March 2008
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
        //    Input, Func< double, double, double, double > func, the name of the 
        //    user supplied function which is to be integrated.
        //
        //    Output, double PYRAMID_UNIT_O18_3D, the approximate integral
        //    of the function.
        //
    {
        int i;
        int order = 18;
        double quad;
        double result;
        double volume;
        double[] w =
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
        double[] x =
        {
            -0.35309846330877704481,
            0.00000000000000000000,
            0.35309846330877704481,
            -0.35309846330877704481,
            0.00000000000000000000,
            0.35309846330877704481,
            -0.35309846330877704481,
            0.00000000000000000000,
            0.35309846330877704481,
            -0.67969709567986745790,
            0.00000000000000000000,
            0.67969709567986745790,
            -0.67969709567986745790,
            0.00000000000000000000,
            0.67969709567986745790,
            -0.67969709567986745790,
            0.00000000000000000000,
            0.67969709567986745790
        };
        double[] y =
        {
            -0.35309846330877704481,
            -0.35309846330877704481,
            -0.35309846330877704481,
            0.00000000000000000000,
            0.00000000000000000000,
            0.00000000000000000000,
            0.35309846330877704481,
            0.35309846330877704481,
            0.35309846330877704481,
            -0.67969709567986745790,
            -0.67969709567986745790,
            -0.67969709567986745790,
            0.00000000000000000000,
            0.00000000000000000000,
            0.00000000000000000000,
            0.67969709567986745790,
            0.67969709567986745790,
            0.67969709567986745790
        };
        double[] z =
        {
            0.544151844011225288800,
            0.544151844011225288800,
            0.544151844011225288800,
            0.544151844011225288800,
            0.544151844011225288800,
            0.544151844011225288800,
            0.544151844011225288800,
            0.544151844011225288800,
            0.544151844011225288800,
            0.12251482265544137787,
            0.12251482265544137787,
            0.12251482265544137787,
            0.12251482265544137787,
            0.12251482265544137787,
            0.12251482265544137787,
            0.12251482265544137787,
            0.12251482265544137787,
            0.12251482265544137787
        };
        //
        //  Quadrature.
        //
        quad = 0.0;
        for (i = 0; i < order; i++)
        {
            quad += w[i] * func(setting, x[i], y[i], z[i]);
        }

        //
        //  Volume.
        //
        volume = pyramid_unit_volume_3d();
        //
        //  Result.
        //
        result = quad * volume;

        return result;
    }

    public static double pyramid_unit_o27_3d(int setting, Func<int, double, double, double, double> func)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID_UNIT_O27_3D approximates an integral inside the unit pyramid in 3D.
        //
        //  Discussion:
        //
        //    A 27 point formula is used.
        //
        //    The (X,Y,Z) integration region can be represented as:
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
        //    22 March 2008
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
        //    Input, Func< double, double, double, double > func, the name of the 
        //    user supplied function which is to be integrated.
        //
        //    Output, double PYRAMID_UNIT_O27_3D, the approximate integral
        //    of the function.
        //
    {
        int i;
        int order = 27;
        double quad;
        double result;
        double volume;
        double[] w =
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
        double[] x =
        {
            -0.7180557413198889387,
            0.00000000000000000000,
            0.7180557413198889387,
            -0.7180557413198889387,
            0.00000000000000000000,
            0.7180557413198889387,
            -0.7180557413198889387,
            0.00000000000000000000,
            0.7180557413198889387,
            -0.50580870785392503961,
            0.00000000000000000000,
            0.50580870785392503961,
            -0.50580870785392503961,
            0.00000000000000000000,
            0.50580870785392503961,
            -0.50580870785392503961,
            0.00000000000000000000,
            0.50580870785392503961,
            -0.22850430565396735360,
            0.00000000000000000000,
            0.22850430565396735360,
            -0.22850430565396735360,
            0.00000000000000000000,
            0.22850430565396735360,
            -0.22850430565396735360,
            0.00000000000000000000,
            0.22850430565396735360
        };
        double[] y =
        {
            -0.7180557413198889387,
            -0.7180557413198889387,
            -0.7180557413198889387,
            0.00000000000000000000,
            0.00000000000000000000,
            0.00000000000000000000,
            0.7180557413198889387,
            0.7180557413198889387,
            0.7180557413198889387,
            -0.50580870785392503961,
            -0.50580870785392503961,
            -0.50580870785392503961,
            0.00000000000000000000,
            0.00000000000000000000,
            0.00000000000000000000,
            0.50580870785392503961,
            0.50580870785392503961,
            0.50580870785392503961,
            -0.22850430565396735360,
            -0.22850430565396735360,
            -0.22850430565396735360,
            0.00000000000000000000,
            0.00000000000000000000,
            0.00000000000000000000,
            0.22850430565396735360,
            0.22850430565396735360,
            0.22850430565396735360
        };
        double[] z =
        {
            0.07299402407314973216,
            0.07299402407314973216,
            0.07299402407314973216,
            0.07299402407314973216,
            0.07299402407314973216,
            0.07299402407314973216,
            0.07299402407314973216,
            0.07299402407314973216,
            0.07299402407314973216,
            0.34700376603835188472,
            0.34700376603835188472,
            0.34700376603835188472,
            0.34700376603835188472,
            0.34700376603835188472,
            0.34700376603835188472,
            0.34700376603835188472,
            0.34700376603835188472,
            0.34700376603835188472,
            0.70500220988849838312,
            0.70500220988849838312,
            0.70500220988849838312,
            0.70500220988849838312,
            0.70500220988849838312,
            0.70500220988849838312,
            0.70500220988849838312,
            0.70500220988849838312,
            0.70500220988849838312
        };
        //
        //  Quadrature.
        //
        quad = 0.0;
        for (i = 0; i < order; i++)
        {
            quad += w[i] * func(setting, x[i], y[i], z[i]);
        }

        //
        //  Volume.
        //
        volume = pyramid_unit_volume_3d();
        //
        //  Result.
        //
        result = quad * volume;

        return result;
    }

    public static double pyramid_unit_o48_3d(int setting, Func<int, double, double, double, double> func)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID_UNIT_O48_3D approximates an integral inside the unit pyramid in 3D.
        //
        //  Discussion:
        //
        //    An 48 point degree 7 formula, Stroud CN:C2:7-1, is used.
        //
        //    The (X,Y,Z) integration region can be represented as:
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
        //    22 March 2008
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
        //    user supplied function which is to be integrated.
        //
        //    Output, double PYRAMID_UNIT_O48_3D, the approximate integral
        //    of the function.
        //
    {
        int i;
        int order = 48;
        double quad;
        double result;
        double volume;
        double[] w =
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
        double[] x =
        {
            0.88091731624450909E+00,
            -0.88091731624450909E+00,
            0.0000000000000000E+00,
            0.0000000000000000E+00,
            0.70491874112648223E+00,
            -0.70491874112648223E+00,
            0.0000000000000000E+00,
            0.0000000000000000E+00,
            0.44712732143189760E+00,
            -0.44712732143189760E+00,
            0.0000000000000000E+00,
            0.0000000000000000E+00,
            0.18900486065123448E+00,
            -0.18900486065123448E+00,
            0.0000000000000000E+00,
            0.0000000000000000E+00,
            0.36209733410322176E+00,
            -0.36209733410322176E+00,
            -0.36209733410322176E+00,
            0.36209733410322176E+00,
            0.76688932060387538E+00,
            -0.76688932060387538E+00,
            -0.76688932060387538E+00,
            0.76688932060387538E+00,
            0.28975386476618070E+00,
            -0.28975386476618070E+00,
            -0.28975386476618070E+00,
            0.28975386476618070E+00,
            0.61367241226233160E+00,
            -0.61367241226233160E+00,
            -0.61367241226233160E+00,
            0.61367241226233160E+00,
            0.18378979287798017E+00,
            -0.18378979287798017E+00,
            -0.18378979287798017E+00,
            0.18378979287798017E+00,
            0.38925011625173161E+00,
            -0.38925011625173161E+00,
            -0.38925011625173161E+00,
            0.38925011625173161E+00,
            7.76896479525748113E-02,
            -7.76896479525748113E-02,
            -7.76896479525748113E-02,
            7.76896479525748113E-02,
            0.16453962988669860E+00,
            -0.16453962988669860E+00,
            -0.16453962988669860E+00,
            0.16453962988669860E+00
        };
        double[] y =
        {
            0.0000000000000000E+00,
            0.0000000000000000E+00,
            0.88091731624450909E+00,
            -0.88091731624450909E+00,
            0.0000000000000000E+00,
            0.0000000000000000E+00,
            0.70491874112648223E+00,
            -0.70491874112648223E+00,
            0.0000000000000000E+00,
            0.0000000000000000E+00,
            0.44712732143189760E+00,
            -0.44712732143189760E+00,
            0.0000000000000000E+00,
            0.0000000000000000E+00,
            0.18900486065123448E+00,
            -0.18900486065123448E+00,
            0.36209733410322176E+00,
            0.36209733410322176E+00,
            -0.36209733410322176E+00,
            -0.36209733410322176E+00,
            0.76688932060387538E+00,
            0.76688932060387538E+00,
            -0.76688932060387538E+00,
            -0.76688932060387538E+00,
            0.28975386476618070E+00,
            0.28975386476618070E+00,
            -0.28975386476618070E+00,
            -0.28975386476618070E+00,
            0.61367241226233160E+00,
            0.61367241226233160E+00,
            -0.61367241226233160E+00,
            -0.61367241226233160E+00,
            0.18378979287798017E+00,
            0.18378979287798017E+00,
            -0.18378979287798017E+00,
            -0.18378979287798017E+00,
            0.38925011625173161E+00,
            0.38925011625173161E+00,
            -0.38925011625173161E+00,
            -0.38925011625173161E+00,
            7.76896479525748113E-02,
            7.76896479525748113E-02,
            -7.76896479525748113E-02,
            -7.76896479525748113E-02,
            0.16453962988669860E+00,
            0.16453962988669860E+00,
            -0.16453962988669860E+00,
            -0.16453962988669860E+00
        };
        double[] z =
        {
            4.85005494469969989E-02,
            4.85005494469969989E-02,
            4.85005494469969989E-02,
            4.85005494469969989E-02,
            0.23860073755186201E+00,
            0.23860073755186201E+00,
            0.23860073755186201E+00,
            0.23860073755186201E+00,
            0.51704729510436798E+00,
            0.51704729510436798E+00,
            0.51704729510436798E+00,
            0.51704729510436798E+00,
            0.79585141789677305E+00,
            0.79585141789677305E+00,
            0.79585141789677305E+00,
            0.79585141789677305E+00,
            4.85005494469969989E-02,
            4.85005494469969989E-02,
            4.85005494469969989E-02,
            4.85005494469969989E-02,
            4.85005494469969989E-02,
            4.85005494469969989E-02,
            4.85005494469969989E-02,
            4.85005494469969989E-02,
            0.23860073755186201E+00,
            0.23860073755186201E+00,
            0.23860073755186201E+00,
            0.23860073755186201E+00,
            0.23860073755186201E+00,
            0.23860073755186201E+00,
            0.23860073755186201E+00,
            0.23860073755186201E+00,
            0.51704729510436798E+00,
            0.51704729510436798E+00,
            0.51704729510436798E+00,
            0.51704729510436798E+00,
            0.51704729510436798E+00,
            0.51704729510436798E+00,
            0.51704729510436798E+00,
            0.51704729510436798E+00,
            0.79585141789677305E+00,
            0.79585141789677305E+00,
            0.79585141789677305E+00,
            0.79585141789677305E+00,
            0.79585141789677305E+00,
            0.79585141789677305E+00,
            0.79585141789677305E+00,
            0.79585141789677305E+00
        };
        //
        //  Quadrature.
        //
        quad = 0.0;
        for (i = 0; i < order; i++)
        {
            quad += w[i] * func(setting, x[i], y[i], z[i]);
        }

        //
        //  Volume.
        //
        volume = pyramid_unit_volume_3d();
        //
        //  Result.
        //
        result = quad * volume;

        return result;
    }

    public static double pyramid_unit_monomial_3d(int alpha, int beta, int gamma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID_UNIT_MONOMIAL_3D: monomial integral in a unit pyramid in 3D.
        //
        //  Discussion:
        //
        //    This function returns the value of the integral of X^ALPHA Y^BETA Z^GAMMA
        //    over the unit pyramid.
        //
        //    The unit pyramid is defined as:
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
        //    Input, int ALPHA, BETA, GAMMA, the exponents of
        //    X, Y and Z in the monomial.
        //
        //    Output, double PYRAMID_UNIT_MONOMIAL_3D, the volume of the pyramid.
        //
    {
        int i;
        int i_hi;
        double value = 0;

        value = 0.0;

        switch (alpha % 2)
        {
            case 0 when beta % 2 == 0:
            {
                i_hi = 2 + alpha + beta;

                for (i = 0; i <= i_hi; i++)
                {
                    value += typeMethods.r8_mop(i) * typeMethods.r8_choose(i_hi, i)
                             / (i + gamma + 1);
                }

                value = value
                    * 2.0 / (alpha + 1)
                    * 2.0 / (beta + 1);
                break;
            }
        }

        return value;
    }

    public static double pyramid_unit_volume_3d()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID_UNIT_VOLUME_3D: volume of a unit pyramid with square base in 3D.
        //
        //  Integration region:
        //
        //    - ( 1 - Z ) <= X <= 1 - Z
        //    - ( 1 - Z ) <= Y <= 1 - Z
        //              0 <= Z <= 1.
        //
        //  Discussion:
        //
        //    The volume of this unit pyramid is 4/3.
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
        //    Output, double PYRAMID_UNIT_VOLUME_3D, the volume of the pyramid.
        //
    {
        double volume;

        volume = 4.0 / 3.0;

        return volume;
    }

    public static double pyramid_volume_3d(double r, double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID_VOLUME_3D returns the volume of a pyramid with square base in 3D.
        //
        //  Integration region:
        //
        //    - ( H - Z ) * R <= X <= ( H - Z ) * R
        //    - ( H - Z ) * R <= Y <= ( H - Z ) * R
        //                  0 <= Z <= H.
        //
        //  Discussion:
        //
        //    A pyramid with square base can be regarded as the upper half of a
        //    3D octahedron.
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
        //    Input, double R, the "radius" of the pyramid, that is, half the
        //    length of one of the sides of the square base.
        //
        //    Input, double H, the height of the pyramid.
        //
        //    Output, double PYRAMID_VOLUME_3D, the volume of the pyramid.
        //
    {
        double value = 0;

        value = 4.0 / 3.0 * h * r * r;

        return value;
    }

}