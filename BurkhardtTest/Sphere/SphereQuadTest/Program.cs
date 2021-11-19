using System;
using Burkardt.SphereNS;

namespace SphereQuadTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SPHERE_QUAD_TEST.
        //
        //  Discussion:
        //
        //    SPHERE_QUAD_TEST tests the SPHERE_QUAD library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SPHERE_QUAD_TEST");
        Console.WriteLine("  Test the SPHERE_QUAD library.");

        int[] e_save = new int[3];

        test01(ref e_save);
        test02(ref e_save);
        test03(ref e_save);
        test04(ref e_save);
        test05(ref e_save);

        Console.WriteLine("");
        Console.WriteLine("SPHERE_QUAD_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(ref int[] e_save)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests SPHERE01_QUAD_LL*.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e = new int[3];
        double exact;
        double h = 0;
        int h_test;
        int i;
        int n_llc = 0;
        int n_llm = 0;
        int n_llv = 0;
        int n_mc = 0;
        double result_llc;
        double result_llm;
        double result_llv;
        double result_mc;
        int seed;

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Approximate the integral of a function on the unit sphere.");
        Console.WriteLine("");
        Console.WriteLine("  SPHERE01_QUAD_MC uses a Monte Carlo method.");
        Console.WriteLine("  SPHERE01_QUAD_LLC uses centroids of spherical triangles.");
        Console.WriteLine("  SPHERE01_QUAD_LLM uses midsides of spherical triangles.");
        Console.WriteLine("  SPHERE01_QUAD_LLV uses vertices of spherical triangles.");
        Console.WriteLine("");
        Console.WriteLine("  H              QUAD_MC       QUAD_LLC      QUAD_LLM      QUAD_LLV         EXACT");

        for (i = 0; i <= 17; i++)
        {
            switch (i)
            {
                case 0:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 0;
                    break;
                case 1:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 0;
                    break;
                case 2:
                    e[0] = 1;
                    e[1] = 0;
                    e[2] = 0;
                    break;
                case 3:
                    e[0] = 0;
                    e[1] = 1;
                    e[2] = 0;
                    break;
                case 4:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 1;
                    break;
                case 5:
                    e[0] = 2;
                    e[1] = 0;
                    e[2] = 0;
                    break;
                case 6:
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 2;
                    break;
                case 7:
                    e[0] = 2;
                    e[1] = 2;
                    e[2] = 2;
                    break;
                case 8:
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 4;
                    break;
                case 9:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 6;
                    break;
                case 10:
                    e[0] = 1;
                    e[1] = 2;
                    e[2] = 4;
                    break;
                case 11:
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 2;
                    break;
                case 12:
                    e[0] = 6;
                    e[1] = 2;
                    e[2] = 0;
                    break;
                case 13:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 8;
                    break;
                case 14:
                    e[0] = 6;
                    e[1] = 0;
                    e[2] = 4;
                    break;
                case 15:
                    e[0] = 4;
                    e[1] = 6;
                    e[2] = 2;
                    break;
                case 16:
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 8;
                    break;
                case 17:
                    e[0] = 16;
                    e[1] = 0;
                    e[2] = 0;
                    break;
            }

            polyterm_exponent("SET", ref e_save, ref e);

            switch (i)
            {
                case 0:
                    Console.WriteLine("");
                    Console.WriteLine("Point counts per method:");
                    break;
                default:
                    polyterm_exponent("PRINT", ref e_save, ref e);
                    break;
            }

            for (h_test = 1; h_test <= 3; h_test++)
            {
                h = h_test switch
                {
                    1 => 1.0,
                    2 => 0.1,
                    3 => 0.01,
                    _ => h
                };

                n_mc = Quad.sphere01_quad_mc_size(h);

                result_mc = Quad.sphere01_quad_mc(polyterm_value_3d, h, ref seed, n_mc);

                result_llc = Quad.sphere01_quad_llc(polyterm_value_3d, h, ref n_llc);

                result_llm = Quad.sphere01_quad_llm(polyterm_value_3d, h, ref n_llm);

                result_llv = Quad.sphere01_quad_llv(polyterm_value_3d, h, ref n_llv);

                exact = Integrals.sphere01_monomial_integral(e);

                switch (i)
                {
                    case 0:
                        Console.WriteLine("  " + h.ToString().PadLeft(12)
                                               + "  " + n_mc.ToString().PadLeft(12)
                                               + "  " + n_llc.ToString().PadLeft(12)
                                               + "  " + n_llm.ToString().PadLeft(12)
                                               + "  " + n_llv.ToString().PadLeft(12) + "");
                        break;
                    default:
                        Console.WriteLine("  " + h.ToString().PadLeft(12)
                                               + "  " + result_mc.ToString().PadLeft(12)
                                               + "  " + result_llc.ToString().PadLeft(12)
                                               + "  " + result_llm.ToString().PadLeft(12)
                                               + "  " + result_llv.ToString().PadLeft(12)
                                               + "  " + exact.ToString().PadLeft(12) + "");
                        break;
                }
            }
        }

    }

    private static void test02(ref int[] e_save)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests SPHERE01_QUAD_ICOS1C.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e = new int[3];
        double error;
        double exact;
        int factor;
        int factor_log;
        int i;
        int n = 0;
        double result;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Approximate the integral of a function on the unit sphere.");
        Console.WriteLine("  SPHERE01_QUAD_ICOS1C uses centroids of spherical triangles.");
        Console.WriteLine("");
        Console.WriteLine("FACTOR         N        QUAD          EXACT         ERROR");

        for (i = 1; i <= 17; i++)
        {
            switch (i)
            {
                case 1:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 0;
                    break;
                case 2:
                    e[0] = 1;
                    e[1] = 0;
                    e[2] = 0;
                    break;
                case 3:
                    e[0] = 0;
                    e[1] = 1;
                    e[2] = 0;
                    break;
                case 4:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 1;
                    break;
                case 5:
                    e[0] = 2;
                    e[1] = 0;
                    e[2] = 0;
                    break;
                case 6:
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 2;
                    break;
                case 7:
                    e[0] = 2;
                    e[1] = 2;
                    e[2] = 2;
                    break;
                case 8:
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 4;
                    break;
                case 9:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 6;
                    break;
                case 10:
                    e[0] = 1;
                    e[1] = 2;
                    e[2] = 4;
                    break;
                case 11:
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 2;
                    break;
                case 12:
                    e[0] = 6;
                    e[1] = 2;
                    e[2] = 0;
                    break;
                case 13:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 8;
                    break;
                case 14:
                    e[0] = 6;
                    e[1] = 0;
                    e[2] = 4;
                    break;
                case 15:
                    e[0] = 4;
                    e[1] = 6;
                    e[2] = 2;
                    break;
                case 16:
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 8;
                    break;
                case 17:
                    e[0] = 16;
                    e[1] = 0;
                    e[2] = 0;
                    break;
            }

            polyterm_exponent("SET", ref e_save, ref e);

            polyterm_exponent("PRINT", ref e_save, ref e);

            factor = 1;
            for (factor_log = 0; factor_log <= 5; factor_log++)
            {
                result = Quad.sphere01_quad_icos1c(factor, polyterm_value_3d, ref n);

                exact = Integrals.sphere01_monomial_integral(e);

                error = Math.Abs(exact - result);

                Console.WriteLine("  " + factor.ToString().PadLeft(4)
                                       + "  " + n.ToString().PadLeft(8)
                                       + "  " + result.ToString().PadLeft(14)
                                       + "  " + exact.ToString().PadLeft(14)
                                       + "  " + error.ToString().PadLeft(14) + "");

                factor *= 2;
            }
        }
    }

    private static void test03(ref int[] e_save)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests SPHERE01_QUAD_ICOS1M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e = new int[3];
        double error;
        double exact;
        int factor;
        int factor_log;
        int i;
        int n = 0;
        double result;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Approximate the integral of a function on the unit sphere.");
        Console.WriteLine("  SPHERE01_QUAD_ICOS1M uses midpoints of sides of spherical triangles.");
        Console.WriteLine("");
        Console.WriteLine("FACTOR         N        QUAD          EXACT         ERROR");

        for (i = 1; i <= 17; i++)
        {
            switch (i)
            {
                case 1:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 0;
                    break;
                case 2:
                    e[0] = 1;
                    e[1] = 0;
                    e[2] = 0;
                    break;
                case 3:
                    e[0] = 0;
                    e[1] = 1;
                    e[2] = 0;
                    break;
                case 4:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 1;
                    break;
                case 5:
                    e[0] = 2;
                    e[1] = 0;
                    e[2] = 0;
                    break;
                case 6:
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 2;
                    break;
                case 7:
                    e[0] = 2;
                    e[1] = 2;
                    e[2] = 2;
                    break;
                case 8:
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 4;
                    break;
                case 9:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 6;
                    break;
                case 10:
                    e[0] = 1;
                    e[1] = 2;
                    e[2] = 4;
                    break;
                case 11:
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 2;
                    break;
                case 12:
                    e[0] = 6;
                    e[1] = 2;
                    e[2] = 0;
                    break;
                case 13:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 8;
                    break;
                case 14:
                    e[0] = 6;
                    e[1] = 0;
                    e[2] = 4;
                    break;
                case 15:
                    e[0] = 4;
                    e[1] = 6;
                    e[2] = 2;
                    break;
                case 16:
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 8;
                    break;
                case 17:
                    e[0] = 16;
                    e[1] = 0;
                    e[2] = 0;
                    break;
            }

            polyterm_exponent("SET", ref e_save, ref e);

            polyterm_exponent("PRINT", ref e_save, ref e);

            factor = 1;
            for (factor_log = 0; factor_log <= 5; factor_log++)
            {
                result = Quad.sphere01_quad_icos1m(factor, polyterm_value_3d, ref n);

                exact = Integrals.sphere01_monomial_integral(e);

                error = Math.Abs(exact - result);

                Console.WriteLine("  " + factor.ToString().PadLeft(4)
                                       + "  " + n.ToString().PadLeft(8)
                                       + "  " + result.ToString().PadLeft(14)
                                       + "  " + exact.ToString().PadLeft(14)
                                       + "  " + error.ToString().PadLeft(14) + "");

                factor *= 2;
            }
        }
    }

    private static void test04(ref int[] e_save)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests SPHERE01_QUAD_ICOS1V.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e = new int[3];
        double error;
        double exact;
        int factor;
        int factor_log;
        int i;
        int n = 0;
        double result;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  Approximate the integral of a function on the unit sphere.");
        Console.WriteLine("  SPHERE01_QUAD_ICOS1V uses vertices of spherical triangles.");
        Console.WriteLine("");
        Console.WriteLine("FACTOR         N        QUAD          EXACT         ERROR");

        for (i = 1; i <= 17; i++)
        {
            switch (i)
            {
                case 1:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 0;
                    break;
                case 2:
                    e[0] = 1;
                    e[1] = 0;
                    e[2] = 0;
                    break;
                case 3:
                    e[0] = 0;
                    e[1] = 1;
                    e[2] = 0;
                    break;
                case 4:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 1;
                    break;
                case 5:
                    e[0] = 2;
                    e[1] = 0;
                    e[2] = 0;
                    break;
                case 6:
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 2;
                    break;
                case 7:
                    e[0] = 2;
                    e[1] = 2;
                    e[2] = 2;
                    break;
                case 8:
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 4;
                    break;
                case 9:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 6;
                    break;
                case 10:
                    e[0] = 1;
                    e[1] = 2;
                    e[2] = 4;
                    break;
                case 11:
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 2;
                    break;
                case 12:
                    e[0] = 6;
                    e[1] = 2;
                    e[2] = 0;
                    break;
                case 13:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 8;
                    break;
                case 14:
                    e[0] = 6;
                    e[1] = 0;
                    e[2] = 4;
                    break;
                case 15:
                    e[0] = 4;
                    e[1] = 6;
                    e[2] = 2;
                    break;
                case 16:
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 8;
                    break;
                case 17:
                    e[0] = 16;
                    e[1] = 0;
                    e[2] = 0;
                    break;
            }

            polyterm_exponent("SET", ref e_save, ref e);

            polyterm_exponent("PRINT", ref e_save, ref e);

            factor = 1;
            for (factor_log = 0; factor_log <= 5; factor_log++)
            {
                result = Quad.sphere01_quad_icos1v(factor, polyterm_value_3d, ref n);

                exact = Integrals.sphere01_monomial_integral(e);

                error = Math.Abs(exact - result);

                Console.WriteLine("  " + factor.ToString().PadLeft(4)
                                       + "  " + n.ToString().PadLeft(8)
                                       + "  " + result.ToString().PadLeft(14)
                                       + "  " + exact.ToString().PadLeft(14)
                                       + "  " + error.ToString().PadLeft(14) + "");

                factor *= 2;
            }
        }
    }

    private static void test05(ref int[] e_save)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests SPHERE01_QUAD_ICOS2V.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e = new int[3];
        double error;
        double exact;
        int factor;
        int factor_log;
        int i;
        int n = 0;
        double result;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  Approximate the integral of a function on the unit sphere.");
        Console.WriteLine("  SPHERE01_QUAD_ICOS2V uses vertices of spherical triangles.");
        Console.WriteLine("");
        Console.WriteLine("FACTOR         N        QUAD          EXACT         ERROR");

        for (i = 1; i <= 17; i++)
        {
            switch (i)
            {
                case 1:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 0;
                    break;
                case 2:
                    e[0] = 1;
                    e[1] = 0;
                    e[2] = 0;
                    break;
                case 3:
                    e[0] = 0;
                    e[1] = 1;
                    e[2] = 0;
                    break;
                case 4:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 1;
                    break;
                case 5:
                    e[0] = 2;
                    e[1] = 0;
                    e[2] = 0;
                    break;
                case 6:
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 2;
                    break;
                case 7:
                    e[0] = 2;
                    e[1] = 2;
                    e[2] = 2;
                    break;
                case 8:
                    e[0] = 0;
                    e[1] = 2;
                    e[2] = 4;
                    break;
                case 9:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 6;
                    break;
                case 10:
                    e[0] = 1;
                    e[1] = 2;
                    e[2] = 4;
                    break;
                case 11:
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 2;
                    break;
                case 12:
                    e[0] = 6;
                    e[1] = 2;
                    e[2] = 0;
                    break;
                case 13:
                    e[0] = 0;
                    e[1] = 0;
                    e[2] = 8;
                    break;
                case 14:
                    e[0] = 6;
                    e[1] = 0;
                    e[2] = 4;
                    break;
                case 15:
                    e[0] = 4;
                    e[1] = 6;
                    e[2] = 2;
                    break;
                case 16:
                    e[0] = 2;
                    e[1] = 4;
                    e[2] = 8;
                    break;
                case 17:
                    e[0] = 16;
                    e[1] = 0;
                    e[2] = 0;
                    break;
            }

            polyterm_exponent("SET", ref e_save, ref e);

            polyterm_exponent("PRINT", ref e_save, ref e);

            factor = 1;
            for (factor_log = 0; factor_log <= 5; factor_log++)
            {
                result = Quad.sphere01_quad_icos2v(factor, polyterm_value_3d, ref n);

                exact = Integrals.sphere01_monomial_integral(e);

                error = Math.Abs(exact - result);

                Console.WriteLine("  " + factor.ToString().PadLeft(4)
                                       + "  " + n.ToString().PadLeft(8)
                                       + "  " + result.ToString().PadLeft(14)
                                       + "  " + exact.ToString().PadLeft(14)
                                       + "  " + error.ToString().PadLeft(14) + "");

                factor *= 2;
            }
        }
    }

    private static void polyterm_exponent(string action, ref int[] e_save, ref int[] e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYTERM_EXPONENT gets or sets the exponents for the polynomial term.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string ACTION.
        //    'GET' asks the routine to return the current values in E.
        //    'SET' asks the routine to set the current values to E.
        //
        //    Input/output, int E[3], storage used to set or get values.
        //
    {
        int i;

        switch (action[0])
        {
            case 'G':
            {
                for (i = 0; i < 3; i++)
                {
                    e[i] = e_save[i];
                }

                break;
            }
            case 'P':
            {
                Console.WriteLine("");

                switch (e_save[0])
                {
                    case 0 when e_save[1] == 0 && e_save[2] == 0:
                        Console.WriteLine("P(X,Y,Z) = 1");
                        break;
                    default:
                    {
                        string cout = "P(X,Y,Z) = ";

                        switch (e_save[0])
                        {
                            case 0:
                                break;
                            case 1:
                                cout += " X";
                                break;
                            default:
                                cout += " X^" + e_save[0];
                                break;
                        }

                        switch (e_save[1])
                        {
                            case 0:
                                break;
                            case 1:
                                cout += " Y";
                                break;
                            default:
                                cout += " Y^" + e_save[1];
                                break;
                        }

                        switch (e_save[2])
                        {
                            case 0:
                                break;
                            case 1:
                                cout += " Z";
                                break;
                            default:
                                cout += " Z^" + e_save[2];
                                break;
                        }

                        Console.WriteLine(cout);
                        break;
                    }
                }

                break;
            }
            case 'S':
            {
                for (i = 0; i < 3; i++)
                {
                    e_save[i] = e[i];
                }

                break;
            }
        }
    }

    private static double[] polyterm_value_3d(int n, double[] x, double[] f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYTERM_VALUE_3D evaluates a single polynomial term in 3D.
        //
        //  Discussion:
        //
        //    The polynomial term has the form:
        //
        //      F(X) = X(1)^E(1) * X(2)^E(2) * X(3)^E(3)
        //
        //    The exponents E(1:3) are set by calling POLYTERM_EXPONENT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[3*N], the points where the polynomial term 
        //    is to be evaluated.
        //
        //    Output, double F[N], the value of the polynomial term.
        //
    {
        int[] e = new int[3];
        int[] e_save = new int[3];
        int i;
        int j;

        polyterm_exponent("GET", ref e_save, ref e);

        for (j = 0; j < n; j++)
        {
            f[j] = 1.0;
        }

        for (i = 0; i < 3; i++)
        {
            if (e[i] != 0)
            {
                for (j = 0; j < n; j++)
                {
                    f[j] *= Math.Pow(x[i + j * 3], e[i]);
                }
            }
        }

        return f;
    }
}