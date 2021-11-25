using System;
using Burkardt.SolveNS;
using Burkardt.Types;

namespace SubsetTestNS;

public static class CongruenceTest
{
    public static void congruence_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONGRUENCE_TEST tests CONGRUENCE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int TEST_NUM = 20;

        int[] a_test =  {
                1027, 1027, 1027, 1027, -1027,
                -1027, -1027, -1027, 6, 0,
                0, 0, 1, 1, 1,
                1024, 0, 0, 5, 2
            }
            ;
        int[] b_test =  {
                712, 712, -712, -712, 712,
                712, -712, -712, 8, 0,
                1, 1, 0, 0, 1,
                -15625, 0, 3, 0, 4
            }
            ;
        int[] c_test =  {
                7, -7, 7, -7, 7,
                -7, 7, -7, 50, 0,
                0, 1, 0, 1, 0,
                11529, 1, 11, 19, 7
            }
            ;
        bool error = false;
        int test_i;

        Console.WriteLine("");
        Console.WriteLine("CONGRUENCE_TEST");
        Console.WriteLine("  CONGRUENCE solves a congruence equation:");
        Console.WriteLine("    A * X = C mod ( B )");
        Console.WriteLine("");
        Console.WriteLine("   I        A         B         C         X     Mod ( A*X-C,B)");
        Console.WriteLine("");

        for (test_i = 1; test_i < TEST_NUM; test_i++)
        {
            int a = a_test[test_i];
            int b = b_test[test_i];
            int c = c_test[test_i];

            int x = Congruence.congruence(a, b, c, ref error);

            switch (error)
            {
                case true:
                    Console.WriteLine("  "
                                      + test_i.ToString().PadLeft(2) + "  "
                                      + a.ToString().PadLeft(10) + "  "
                                      + b.ToString().PadLeft(10) + "  "
                                      + c.ToString().PadLeft(10) + "  "
                                      + "(An error occurred)");
                    break;
                default:
                {
                    int result = b != 0 ? typeMethods.i4_modp(a * x - c, b) : 0;

                    Console.WriteLine("  "
                                      + test_i.ToString().PadLeft(2) + "  "
                                      + a.ToString().PadLeft(10) + "  "
                                      + b.ToString().PadLeft(10) + "  "
                                      + c.ToString().PadLeft(10) + "  "
                                      + x.ToString().PadLeft(10) + "  "
                                      + result.ToString().PadLeft(10) + "");
                    break;
                }
            }

        }
    }

}