using System;
using Burkardt.Stroud;

namespace StroudTest;

public static class tests
{
    public static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests BALL_F1_ND, BALL_F3_ND.
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
    {
        double[] center;
        int i;
        int n;
        const int N_MAX = 3;
        string name = "";
        int num;
        double result1 = 0;
        double result2 = 0;
        double r;

        r = 2.0;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  For integrals in a ball in ND:");
        Console.WriteLine("  BALL_F1_ND approximates the integral;");
        Console.WriteLine("  BALL_F3_ND approximates the integral.");

        for (n = 2; n <= N_MAX; n++)
        {
            center = new double[n];
            for (i = 0; i < n; i++)
            {
                center[i] = i + 1;
            }

            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension N = " + n + "");
            Console.WriteLine("  Ball center:");
            for (i = 0; i < n; i++)
            {
            }

            Console.WriteLine("");
            Console.WriteLine("  Ball radius = " + r + "");
            Console.WriteLine("  Ball volume = " + Ball.ball_volume_nd(n, r) + "");
            Console.WriteLine("");
            Console.WriteLine("    Rule:      F1          F3");
            Console.WriteLine("    F(X)");
            Console.WriteLine("");

            num = functions.function_nd_num();

            for (i = 1; i <= num; i++)
            {
                int function_nd_index = i;
                functions.function_nd_name(function_nd_index, ref name);

                result1 = Ball.ball_f1_nd(function_nd_index, functions.function_nd, n, center, r);
                result2 = Ball.ball_f3_nd(function_nd_index, functions.function_nd, n, center, r);
                Console.WriteLine("  " + name
                                       + "  " + result1.ToString().PadLeft(14)
                                       + "  " + result2.ToString().PadLeft(14) + "");
            }

        }
    }

    public static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests BALL_MONOMIAL_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] center = { 0.0, 0.0, 0.0 };
        int dim_num = 3;
        int[] p = new int[3];
        double result1 = 0;
        double result2 = 0;
        double r = 2.0;
        string string_ = "";
        int test;
        int test_num = 4;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  For the integral of a monomial in a ball in ND:");
        Console.WriteLine("  BALL_MONOMIAL_ND approximates the integral.");
        Console.WriteLine("  BALL_F1_ND, which can handle general integrands,");
        Console.WriteLine("    will be used for comparison.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension N = " + dim_num + "");
        Console.WriteLine("  Ball radius = " + r + "");
        Console.WriteLine("  Ball volume = " + Ball.ball_volume_nd(dim_num, r) + "");
        Console.WriteLine("");
        Console.WriteLine("    Rule:     MONOMIAL    F1");
        Console.WriteLine("    F(X)");
        Console.WriteLine("");

        for (test = 1; test <= test_num; test++)
        {
            switch (test)
            {
                case 1:
                    string_ = "         1";
                    p[0] = 0;
                    p[1] = 0;
                    p[2] = 0;
                    result2 = Ball.ball_f1_nd(0, mono.mono_000_3d, dim_num, center, r);
                    break;
                case 2:
                    string_ = "       xyz";
                    p[0] = 1;
                    p[1] = 1;
                    p[2] = 1;
                    result2 = Ball.ball_f1_nd(0, mono.mono_111_3d, dim_num, center, r);
                    break;
                case 3:
                    string_ = "   x^2 z^2";
                    p[0] = 2;
                    p[1] = 0;
                    p[2] = 2;
                    result2 = Ball.ball_f1_nd(0, mono.mono_202_3d, dim_num, center, r);
                    break;
                case 4:
                    string_ = " x^4y^2z^2";
                    p[0] = 4;
                    p[1] = 2;
                    p[2] = 2;
                    result2 = Ball.ball_f1_nd(0, mono.mono_422_3d, dim_num, center, r);
                    break;
            }

            result1 = Ball.ball_monomial_nd(dim_num, p, r);

            Console.WriteLine("  " + string_
                                   + "  " + result1.ToString().PadLeft(14)
                                   + "  " + result2.ToString().PadLeft(14) + "");
        }

    }

    public static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests BALL_UNIT_**_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        string name = "";
        int num;
        double result1 = 0;
        double result2 = 0;
        double result3;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  For integrals in the unit ball in 3D:");
        Console.WriteLine("  BALL_UNIT_07_3D uses a formula of degree 7;");
        Console.WriteLine("  BALL_UNIT_14_3D uses a formula of degree 14;");
        Console.WriteLine("  BALL_UNIT_15_3D uses a formula of degree 15.");
        Console.WriteLine("");
        Console.WriteLine("  Unit ball volume = " + Ball.ball_unit_volume_nd(3) + "");
        Console.WriteLine("");
        Console.WriteLine("    Rule:      #7             #14           #15");
        Console.WriteLine("    F(X)");
        Console.WriteLine("");

        num = functions.function_3d_num();

        for (i = 1; i <= num; i++)
        {
            int function_3d_index = i;
            functions.function_3d_name(function_3d_index, ref name);

            result1 = Ball.ball_unit_07_3d(function_3d_index, functions.function_3d);
            result2 = Ball.ball_unit_14_3d(function_3d_index, functions.function_3d);
            result3 = Ball.ball_unit_15_3d(function_3d_index, functions.function_3d);

            Console.WriteLine("  " + name
                                   + "  " + result1.ToString().PadLeft(14)
                                   + "  " + result2.ToString().PadLeft(14)
                                   + "  " + result3.ToString().PadLeft(14) + "");
        }
    }

    public static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests BALL_UNIT_F1_ND, BALL_UNIT_F3_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int n;
        const int N_MAX = 3;
        string name = "";
        int num;
        double result1 = 0;
        double result2 = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  For integrals inside the unit ball in ND:");
        Console.WriteLine("  BALL_UNIT_F1_ND approximates the integral;");
        Console.WriteLine("  BALL_UNIT_F3_ND approximates the integral.");
        Console.WriteLine("");

        for (n = 2; n <= N_MAX; n++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension N = " + n + "");
            Console.WriteLine("  Unit ball volume = " + Ball.ball_unit_volume_nd(n) + "");
            Console.WriteLine("");
            Console.WriteLine("");
            Console.WriteLine("    Rule:      F1          F3");
            Console.WriteLine("    F(X)");
            Console.WriteLine("");

            num = functions.function_nd_num();

            for (i = 1; i <= num; i++)
            {
                int function_nd_index = i;
                functions.function_nd_name(function_nd_index, ref name);

                result1 = Ball.ball_unit_f1_nd(function_nd_index, functions.function_nd, n);
                result2 = Ball.ball_unit_f3_nd(function_nd_index, functions.function_nd, n);

                Console.WriteLine("  " + name
                                       + result1.ToString().PadLeft(14)
                                       + result2.ToString().PadLeft(14) + "");
            }
        }
    }

    public static void test045()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST045 tests BALL_UNIT_VOLUME_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim_num = 3;

        Console.WriteLine("");
        Console.WriteLine("TEST045");
        Console.WriteLine("  In 3 dimensions:");
        Console.WriteLine("  BALL_UNIT_VOLUME_3D gets the volume of the unit ball.");
        Console.WriteLine("  BALL_UNIT_VOLUME_ND will be called for comparison.");
        Console.WriteLine("");
        Console.WriteLine("    N    Volume    Method");
        Console.WriteLine("");

        Console.WriteLine("  " + dim_num.ToString().PadLeft(3)
                               + "  " + Ball.ball_unit_volume_3d().ToString().PadLeft(14)
                               + "  BALL_UNIT_VOLUME_3D");

        Console.WriteLine("  " + dim_num.ToString().PadLeft(3)
                               + "  " + Ball.ball_unit_volume_nd(dim_num).ToString().PadLeft(14)
                               + "  BALL_UNIT_VOLUME_ND");

    }

    public static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests BALL_UNIT_VOLUME_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim_num;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  BALL_UNIT_VOLUME_ND computes the volume");
        Console.WriteLine("    of the unit ball in ND.");
        Console.WriteLine("");
        Console.WriteLine("    N    Volume");
        Console.WriteLine("");

        for (dim_num = 2; dim_num <= 10; dim_num++)
        {
            Console.WriteLine("  " + dim_num.ToString().PadLeft(3)
                                   + "  " + Ball.ball_unit_volume_nd(dim_num) + "");
        }
    }

    public static void test052()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST052 tests BALL_VOLUME_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int n = 3;
        double r;

        Console.WriteLine("");
        Console.WriteLine("TEST052");
        Console.WriteLine("  In 3 dimensions:");
        Console.WriteLine("  BALL_VOLUME_3D computes the volume of a unit ball.");
        Console.WriteLine("  BALL_VOLUME_ND will be called for comparison.");
        Console.WriteLine("");
        Console.WriteLine("    N    R      Volume    Method");
        Console.WriteLine("");

        r = 1.0;

        for (i = 1; i <= 3; i++)
        {
            Console.WriteLine("  " + n.ToString().PadLeft(3)
                                   + "  " + r.ToString().PadLeft(14)
                                   + "  " + Ball.ball_volume_3d(r)
                                   + "  BALL_VOLUME_3D");

            Console.WriteLine("  " + n.ToString().PadLeft(3)
                                   + "  " + r.ToString().PadLeft(14)
                                   + "  " + Ball.ball_volume_nd(n, r)
                                   + "  " + "BALL_VOLUME_ND");

            r *= 2.0;
        }
    }

    public static void test054()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST054 tests BALL_VOLUME_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int n;
        double r;

        Console.WriteLine("");
        Console.WriteLine("TEST054");
        Console.WriteLine("  BALL_UNIT_VOLUME_ND computes the volume of");
        Console.WriteLine("    the unit ball in N dimensions.");
        Console.WriteLine("");
        Console.WriteLine("    N        R      Volume");
        Console.WriteLine("");

        for (n = 2; n <= 10; n++)
        {
            r = 0.5;
            for (i = 1; i <= 3; i++)
            {
                Console.WriteLine("  " + n.ToString().PadLeft(3)
                                       + "  " + r.ToString().PadLeft(14)
                                       + "  " + Ball.ball_volume_nd(n, r) + "");
                r *= 2.0;
            }
        }
    }

    public static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests CIRCLE_ANNULUS.
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
    {
        double area;
        double[] center = new double[2];
        double[] center_test =
            {
                0.0, 0.0,
                0.0, 0.0
            }
            ;
        int dim;
        int dim_num = 2;
        int i;
        int j;
        string name = "";
        int num;
        int nr;
        double radius1;
        double[] radius1_test =
            {
                0.0, 1.0
            }
            ;
        double radius2;
        double[] radius2_test =
            {
                1.0, 2.0
            }
            ;
        double result;
        int test_num = 2;

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  CIRCLE_ANNULUS estimates integrals");
        Console.WriteLine("    in a circular annulus.");
        Console.WriteLine("");
        Console.WriteLine("        F       CENTER         Radius1   Radius2   NR  Result");
        Console.WriteLine("");

        for (i = 0; i < test_num; i++)
        {
            for (dim = 0; dim < dim_num; dim++)
            {
                center[dim] = center_test[dim + i * dim_num];
            }

            radius1 = radius1_test[i];
            radius2 = radius2_test[i];

            area = Circle.circle_annulus_area_2d(radius1, radius2);

            Console.WriteLine("");
            Console.WriteLine("  " + "   Area"
                                   + center[0].ToString().PadLeft(10)
                                   + center[1].ToString().PadLeft(10)
                                   + radius1.ToString().PadLeft(10)
                                   + radius2.ToString().PadLeft(10)
                                   + "  " + area.ToString().PadLeft(10) + "");

            num = functions.function_2d_num();

            for (j = 1; j <= num; j++)
            {
                int function_2d_index = j;

                for (nr = 1; nr <= 4; nr++)
                {
                    result = Circle.circle_annulus(function_2d_index, functions.function_2d, center, radius1,
                        radius2, nr);

                    functions.function_2d_name(function_2d_index, ref name);

                    Console.WriteLine("  " + name
                                           + center[0].ToString().PadLeft(10)
                                           + center[1].ToString().PadLeft(10)
                                           + radius1.ToString().PadLeft(10)
                                           + radius2.ToString().PadLeft(10)
                                           + nr.ToString().PadLeft(2)
                                           + result.ToString().PadLeft(10) + "");
                }
            }
        }
    }

    public static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests CIRCLE_ANNULUS, CIRCLE_RT_SET, CIRCLE_RT_SUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double area;
        double[] center = new double[2];
        double[] center_test =
            {
                0.0, 0.0,
                0.0, 0.0,
                0.0, 0.0
            }
            ;
        int dim = 0;
        int dim_num = 2;
        int i;
        int j;
        string name = "";
        int nc = 0;
        int num = 0;
        int nr = 0;
        int nr2 = 0;
        int nt = 0;
        double[] ra = new double[5];
        double radius1;
        double[] radius1_test =
            {
                0.0, 1.0, 1.0
            }
            ;
        double radius2;
        double[] radius2_test =
            {
                1.0, 2.0, 3.0
            }
            ;
        double result1 = 0;
        double result2 = 0;
        double result3;
        double[] rw = new double[5];
        int rule = 0;
        double[] ta = new double[20];
        int test_num = 3;
        double[] tw = new double[20];
        double zw = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  CIRCLE_ANNULUS estimates integrals in a");
        Console.WriteLine("    circular annulus.");
        Console.WriteLine("  CIRCLE_RT_SET sets up a rule for a circle;");
        Console.WriteLine("  CIRCLE_RT_SUM applies the rule.");
        Console.WriteLine("");
        Console.WriteLine("  RESULT1 = CIRCLE_ANNULUS result.");
        Console.WriteLine("  RESULT2 = Difference of two CIRCLE_RT_SUM results.");
        Console.WriteLine("");
        Console.WriteLine("        F      CENTER       Radius1   Radius2   Result1 Result2");
        Console.WriteLine("");

        for (i = 0; i < test_num; i++)
        {
            for (dim = 0; dim < dim_num; dim++)
            {
                center[dim] = center_test[dim + i * dim_num];
            }

            radius1 = radius1_test[i];
            radius2 = radius2_test[i];

            area = Circle.circle_annulus_area_2d(radius1, radius2);

            name = "   Area";
            Console.WriteLine("");
            Console.WriteLine("  " + name
                                   + center[0].ToString().PadLeft(11)
                                   + center[1].ToString().PadLeft(11)
                                   + radius1.ToString().PadLeft(11)
                                   + radius2.ToString().PadLeft(11)
                                   + area.ToString().PadLeft(11) + "");

            rule = 9;
            Circle.circle_rt_size(rule, ref nr2, ref nt, ref nc);
            Circle.circle_rt_set(rule, nr2, nt, nc, ref ra, ref rw, ref ta, ref tw, ref zw);

            num = functions.function_2d_num();

            for (j = 1; j <= num; j++)
            {
                int function_2d_index = j;
                functions.function_2d_name(function_2d_index, ref name);

                nr = 5;
                result1 = Circle.circle_annulus(function_2d_index, functions.function_2d, center, radius1, radius2,
                    nr);

                result2 = Circle.circle_rt_sum(function_2d_index, functions.function_2d, center, radius1, nr2, ra,
                    rw, nt,
                    ta, tw, zw);

                result3 = Circle.circle_rt_sum(function_2d_index, functions.function_2d, center, radius2, nr2, ra,
                    rw, nt,
                    ta, tw, zw);

                Console.WriteLine("  " + name
                                       + center[0].ToString().PadLeft(11)
                                       + center[1].ToString().PadLeft(11)
                                       + radius1.ToString().PadLeft(11)
                                       + radius2.ToString().PadLeft(11)
                                       + result1.ToString().PadLeft(11)
                                       + (result3 - result2).ToString().PadLeft(11) + "");
            }
        }
    }

    public static void test085()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST085 tests CIRCLE_ANNULUS_AREA_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double area;
        double[] center = new double[2];
        double[] center_test =
        {
            0.0, 0.0,
            1.0, 0.0,
            3.0, 4.0
        };
        int dim;
        int dim_num = 2;
        int i;
        int ntest = 3;
        double radius1;
        double[] radius1_test = { 0.0, 1.0, 1.0 };
        double radius2;
        double[] radius2_test = { 1.0, 2.0, 3.0 };

        Console.WriteLine("");
        Console.WriteLine("TEST085");
        Console.WriteLine("  CIRCLE_ANNULUS_AREA_2D computes the area of a");
        Console.WriteLine("    circular annulus.");
        Console.WriteLine("");
        Console.WriteLine("      CENTER       Radius1   Radius2   Area");
        Console.WriteLine("");

        for (i = 0; i < ntest; i++)
        {
            for (dim = 0; dim < dim_num; dim++)
            {
                center[dim] = center_test[dim + i * dim_num];
            }

            radius1 = radius1_test[i];
            radius2 = radius2_test[i];

            area = Circle.circle_annulus_area_2d(radius1, radius2);

            Console.WriteLine("");
            Console.WriteLine("  " + center[0].ToString().PadLeft(9)
                                   + "  " + center[1].ToString().PadLeft(9)
                                   + "  " + radius1.ToString().PadLeft(9)
                                   + "  " + radius2.ToString().PadLeft(9)
                                   + "  " + area.ToString().PadLeft(9) + "");
        }

    }

    public static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests CIRCLE_ANNULUS_SECTOR, CIRCLE_RT_SET, CIRCLE_RT_SUM.
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
    {
        double as1 = 0;
        double as2 = 0;
        double as3 = 0;
        double as4 = 0;
        double[] center = new double[2];
        int j;
        int nc = 0;
        string name = "";
        int num = 0;
        int nr = 0;
        int nr2 = 0;
        int nt = 0;
        double pi = 3.141592653589793;
        double[] ra = new double[5];
        double radius = 0;
        double radius1a = 0;
        double radius2a = 0;
        double radius1b = 0;
        double radius2b = 0;
        double radius1c = 0;
        double radius2c = 0;
        double radius1d = 0;
        double radius2d = 0;
        double result1 = 0;
        double result2 = 0;
        int rule = 0;
        double[] rw = new double[5];
        double[] ta = new double[20];
        double theta1a = 0;
        double theta2a = 0;
        double theta1b = 0;
        double theta2b = 0;
        double theta1c = 0;
        double theta2c = 0;
        double theta1d = 0;
        double theta2d = 0;
        double[] tw = new double[20];
        double zw = 0;

        nr = 5;

        rule = 9;
        Circle.circle_rt_size(rule, ref nr2, ref nt, ref nc);
        Circle.circle_rt_set(rule, nr2, nt, nc, ref ra, ref rw, ref ta, ref tw, ref zw);

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  CIRCLE_ANNULUS_SECTOR estimates an integral in a");
        Console.WriteLine("    circular annulus sector.");
        Console.WriteLine("  CIRCLE_RT_SET sets an integration rule in a circle.");
        Console.WriteLine("  CIRCLE_RT_SUM uses an integration rule in a circle.");
        Console.WriteLine("");
        Console.WriteLine("  To test CIRCLE_ANNULUS_SECTOR, we estimate an integral");
        Console.WriteLine("  over 4 annular sectors that make up the unit circle,");
        Console.WriteLine("  and add to get RESULT1.");
        Console.WriteLine("");
        Console.WriteLine("  We will also estimate the integral over the unit circle");
        Console.WriteLine("  using CIRCLE_RT_SET and CIRCLE_RT_SUM to get RESULT2.");
        Console.WriteLine("");
        Console.WriteLine("  We will then compare RESULT1 and RESULT2.");
        Console.WriteLine("");
        Console.WriteLine("  CIRCLE_ANNULUS_SECTOR computations will use NR = " + nr + "");
        Console.WriteLine("  CIRCLE_RT_SET/CIRCLE_RT_SUM will use rule " + rule + "");
        Console.WriteLine("");
        Console.WriteLine("  RESULT1 is the sum of Annulus Sector calculations.");
        Console.WriteLine("  RESULT2 is for CIRCLE_RT_SET/CIRCLE_RT_SUM.");
        Console.WriteLine("");

        center[0] = 0.0;
        center[1] = 0.0;
        radius = 1.0;

        radius1a = 0.0;
        radius2a = 0.25;
        theta1a = 0.0;
        theta2a = 0.5 * pi;

        radius1b = 0.0;
        radius2b = 0.25;
        theta1b = 0.5 * pi;
        theta2b = 2.0 * pi;

        radius1c = 0.25;
        radius2c = 1.0;
        theta1c = 0.0;
        theta2c = 0.25 * pi;

        radius1d = 0.25;
        radius2d = 1.0;
        theta1d = 0.25 * pi;
        theta2d = 2.0 * pi;

        Console.WriteLine("");
        Console.WriteLine("       F  Result1  Result2");
        Console.WriteLine("");

        num = functions.function_2d_num();

        for (j = 1; j <= num; j++)
        {
            int function_2d_index = j;

            as1 = Circle.circle_annulus_sector(function_2d_index, functions.function_2d, center, radius1a, radius2a,
                theta1a,
                theta2a, nr);

            as2 = Circle.circle_annulus_sector(function_2d_index, functions.function_2d, center, radius1b, radius2b,
                theta1b,
                theta2b, nr);

            as3 = Circle.circle_annulus_sector(function_2d_index, functions.function_2d, center, radius1c, radius2c,
                theta1c,
                theta2c, nr);

            as4 = Circle.circle_annulus_sector(function_2d_index, functions.function_2d, center, radius1d, radius2d,
                theta1d,
                theta2d, nr);

            result1 = as1 + as2 + as3 + as4;

            result2 = Circle.circle_rt_sum(function_2d_index, functions.function_2d, center, radius, nr2, ra, rw,
                nt,
                ta, tw, zw);

            functions.function_2d_name(function_2d_index, ref name);

            Console.WriteLine("  " + name
                                   + result1.ToString().PadLeft(14)
                                   + result2.ToString().PadLeft(14) + "");
        }

    }
}