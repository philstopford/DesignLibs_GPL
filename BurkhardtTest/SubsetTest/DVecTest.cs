using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SubsetTest
{
    public static class DVecTest
    {
        public static void dvec_add_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DVEC_ADD_TEST tests DVEC_ADD;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 10;

            int[] dvec1 = new int[N];
            int[] dvec2 = new int[N];
            int[] dvec3 = new int[N];
            int i;
            int j;
            int k;
            int l;
            int seed = 123456789;
            int test;
            int test_num = 10;

            Console.WriteLine("");
            Console.WriteLine("DVEC_ADD_TEST");
            Console.WriteLine("  DVEC_ADD adds decimal vectors representing integers;");
            Console.WriteLine("");
            Console.WriteLine("        I        J        I + J    DVEC_ADD");
            Console.WriteLine("");

            for (test = 1; test <= test_num; test++)
            {
                i = UniformRNG.i4_uniform_ab(-100, 100, ref seed);
                j = UniformRNG.i4_uniform_ab(-100, 100, ref seed);

                k = i + j;

                typeMethods.i4_to_dvec(i, N, ref dvec1);
                typeMethods.i4_to_dvec(j, N, ref dvec2);
                typeMethods.dvec_add(N, dvec1, dvec2, ref dvec3);
                l = typeMethods.dvec_to_i4(N, ref dvec3);

                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                  + "  " + j.ToString().PadLeft(8)
                                  + "  " + k.ToString().PadLeft(8)
                                  + "  " + l.ToString().PadLeft(8) + "");
            }

        }

        public static void dvec_complementx_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DVEC_COMPLEMENTX_TEST tests DVEC_COMPLEMENTX;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 10;

            int[] dvec1 = new int[N];
            int[] dvec2 = new int[N];
            int i;
            int j;
            int seed = 123456789;
            int test;
            int test_num = 5;

            Console.WriteLine("");
            Console.WriteLine("DVEC_COMPLEMENTX_TEST");
            Console.WriteLine("  DVEC_COMPLEMENTX returns the ten's complement");
            Console.WriteLine("  of a (signed) decimal vector;");
            Console.WriteLine("");

            for (test = 1; test <= test_num; test++)
            {
                i = UniformRNG.i4_uniform_ab(-100, 100, ref seed);

                typeMethods.i4_to_dvec(i, N, ref dvec1);

                typeMethods.dvec_complementx(N, dvec1, ref dvec2);

                j = typeMethods.dvec_to_i4(N, ref dvec2);

                Console.WriteLine("");
                Console.WriteLine("  I = " + "  " + i + "");
                Console.WriteLine("  J = " + "  " + j + "");
                typeMethods.dvec_print(N, dvec1, "");
                typeMethods.dvec_print(N, dvec2, "");

            }
        }

        public static void dvec_mul_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DVEC_MUL_TEST tests DVEC_MUL;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 10;

            int[] dvec1 = new int[N];
            int[] dvec2 = new int[N];
            int[] dvec3 = new int[N];
            int i;
            int j;
            int k;
            int l;
            int n2 = 0;
            int seed = 123456789;
            int test;
            int test_num = 10;
            int test2;
            int test2_num = 2;

            Console.WriteLine("");
            Console.WriteLine("DVEC_MUL_TEST");
            Console.WriteLine("  DVEC_MUL multiplies decimal vectors");
            Console.WriteLine("  representing integers;");

            for (test2 = 1; test2 <= test2_num; test2++)
            {
                if (test2 == 1)
                {
                    n2 = N;
                }
                else if (test2 == 2)
                {
                    n2 = 6;

                    Console.WriteLine("");
                    Console.WriteLine("  NOW REPEAT THE TEST...");
                    Console.WriteLine("");
                    Console.WriteLine("  but use too few digits to represent big products.");
                    Console.WriteLine("  This corresponds to an \"overflow\".");
                    Console.WriteLine("  The result here should get the final decimal");
                    Console.WriteLine("  digits correctly, though.");
                }

                Console.WriteLine("");
                Console.WriteLine("        I        J        I * J  DVEC_MUL");
                Console.WriteLine("");

                for (test = 1; test <= test_num; test++)
                {
                    i = UniformRNG.i4_uniform_ab(-1000, 1000, ref seed);
                    j = UniformRNG.i4_uniform_ab(-1000, 1000, ref seed);

                    k = i * j;

                    typeMethods.i4_to_dvec(i, n2, ref dvec1);
                    typeMethods.i4_to_dvec(j, n2, ref dvec2);
                    typeMethods.dvec_mul(n2, dvec1, dvec2, ref dvec3);
                    l = typeMethods.dvec_to_i4(n2, ref dvec3);

                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                      + "  " + j.ToString().PadLeft(8)
                                      + "  " + k.ToString().PadLeft(8)
                                      + "  " + l.ToString().PadLeft(8) + "");
                }
            }
        }

        public static void dvec_print_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DVEC_PRINT_TEST tests DVEC_PRINT;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] dvec =  {
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                3, 4, 1, 7, 7, 5, 5, 0, 0, 9
            }
            ;
            int n = 20;

            Console.WriteLine("");
            Console.WriteLine("DVEC_PRINT_TEST");
            Console.WriteLine("  DVEC_PRINT prints a (signed) decimal vector;");

            typeMethods.dvec_print(n, dvec, "  The DVEC:");
        }

        public static void dvec_sub_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DVEC_SUB_TEST tests DVEC_SUB;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 10;

            int[] dvec1 = new int[N];
            int[] dvec2 = new int[N];
            int[] dvec4 = new int[N];
            int i;
            int j;
            int k;
            int l;
            int seed = 123456789;
            int test;
            int test_num = 10;

            Console.WriteLine("");
            Console.WriteLine("DVEC_SUB_TEST");
            Console.WriteLine("  DVEC_SUB subtracts decimal vectors representing integers;");
            Console.WriteLine("");
            Console.WriteLine("        I        J        I - J    DVEC_SUB");
            Console.WriteLine("");

            for (test = 1; test <= test_num; test++)
            {
                i = UniformRNG.i4_uniform_ab(-100, 100, ref seed);
                j = UniformRNG.i4_uniform_ab(-100, 100, ref seed);

                k = i - j;

                typeMethods.i4_to_dvec(i, N, ref dvec1);
                typeMethods.i4_to_dvec(j, N, ref dvec2);
                typeMethods.dvec_sub(N, dvec1, dvec2, ref dvec4);
                l = typeMethods.dvec_to_i4(N, ref dvec4);

                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                  + "  " + j.ToString().PadLeft(8)
                                  + "  " + k.ToString().PadLeft(8)
                                  + "  " + l.ToString().PadLeft(8) + "");
            }
            
        }

        public static void dvec_to_i4_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DVEC_TO_I4_TEST tests DVEC_TO_I4;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] dvec = new int[6];
            int i;
            int i1;
            int i2;
            int n;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("DVEC_TO_I4_TEST");
            Console.WriteLine("  DVEC_TO_I4 converts a DVEC to an I4;");
            Console.WriteLine("");
            Console.WriteLine("         I4 => DVEC => I4");
            Console.WriteLine("");

            seed = 123456789;
            i1 = UniformRNG.i4_uniform_ab(-10000, 10000, ref seed);

            n = 6;
            typeMethods.i4_to_dvec(i1, n, ref dvec);

            i2 = typeMethods.dvec_to_i4(n, ref dvec);

            string cout = "  " + i1.ToString().PadLeft(6) + "  ";
            for (i = n - 1; 0 <= i; i--)
            {
                cout += dvec[i].ToString().PadLeft(2);
            }

            Console.WriteLine(cout + "  " + i2.ToString().PadLeft(6) + "");
        }

    }
}