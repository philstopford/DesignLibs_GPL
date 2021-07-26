using System;
using Burkardt.NiederreiterNS;

namespace Niederreiter2Test
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for NIEDERREITER2_TEST.
            //
            //  Discussion:
            //
            //    NIEDERREITER2_TEST tests the NIEDERREITER2 library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 March 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("NIEDERREITER2_TEST");
            Console.WriteLine("  Test the NIEDERREITER2 library.");

            test01();
            test02();

            Console.WriteLine("");
            Console.WriteLine("NIEDERREITER2_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests NIEDERREITER2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 March 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int dim_num;
            int i;
            int j;
            double[] r;
            int seed;
            int seed_in;
            int seed_out;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  NIEDERREITER2 computes the next element of ");
            Console.WriteLine("  a Niederreiter quasirandom sequence using base 2.");
            Console.WriteLine("");
            Console.WriteLine("  In this test, we call NIEDERREITER2 repeatedly.");

            for (dim_num = 2; dim_num <= 4; dim_num++)
            {
                r = new double[dim_num];

                seed = 0;

                Console.WriteLine("");
                Console.WriteLine("  Using dimension DIM_NUM =   " + dim_num + "");
                Console.WriteLine("");
                Console.WriteLine("  Seed  Seed   Niederreiter2");
                Console.WriteLine("  In    Out");
                Console.WriteLine("");

                for (i = 0; i <= 110; i++)
                {
                    seed_in = seed;
                    Niederreiter2.niederreiter2(dim_num, ref seed, ref r);
                    seed_out = seed;
                    if (i <= 11 || 95 <= i)
                    {
                        string cout = seed_in.ToString().PadLeft(6) + "  ";
                        cout += seed_out.ToString().PadLeft(6) + " ";
                        for (j = 0; j < dim_num; j++)
                        {
                            cout += r[j].ToString().PadLeft(10) + "  ";
                        }

                        Console.WriteLine(cout);
                    }
                    else if (i == 12)
                    {
                        Console.WriteLine("......................");
                    }

                }
            }
        }

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests NIEDERREITER2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 March 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;

            int i;
            int j;
            double[] r = new double[DIM_NUM];
            int seed;
            int seed_in;
            int seed_out;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  NIEDERREITER2 computes the next element of ");
            Console.WriteLine("  a Niederreiter quasirandom sequence using base 2.");
            Console.WriteLine("");
            Console.WriteLine("  In this test, we demonstrate how the SEED can be");
            Console.WriteLine("  manipulated to skip ahead in the sequence, or");
            Console.WriteLine("  to come back to any part of the sequence.");

            Console.WriteLine("");
            Console.WriteLine("  Using dimension DIM_NUM =   " + DIM_NUM + "");

            seed = 0;

            Console.WriteLine("");
            Console.WriteLine("  Seed  Seed   Niederreiter2");
            Console.WriteLine("  In    Out");
            Console.WriteLine("");

            for (i = 0; i <= 10; i++)
            {
                seed_in = seed;
                Niederreiter2.niederreiter2(DIM_NUM, ref seed, ref r);
                seed_out = seed;
                string cout = seed_in.ToString().PadLeft(6) + "  ";
                cout += seed_out.ToString().PadLeft(6) + " ";
                for (j = 0; j < DIM_NUM; j++)
                {
                    cout += r[j].ToString().PadLeft(10) + "  ";
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  Jump ahead by increasing SEED:");

            seed = 100;

            Console.WriteLine("");
            Console.WriteLine("  Seed  Seed   Niederreiter2");
            Console.WriteLine("  In    Out");
            Console.WriteLine("");

            for (i = 0; i <= 5; i++)
            {
                seed_in = seed;
                Niederreiter2.niederreiter2(DIM_NUM, ref seed, ref r);
                seed_out = seed;
                seed_out = seed;
                string cout = seed_in.ToString().PadLeft(6) + "  ";
                cout += seed_out.ToString().PadLeft(6) + " ";
                for (j = 0; j < DIM_NUM; j++)
                {
                    cout += r[j].ToString().PadLeft(10) + "  ";
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  Jump back by decreasing SEED:");

            seed = 3;

            Console.WriteLine("");
            Console.WriteLine("  Seed  Seed   Niederreiter2");
            Console.WriteLine("  In    Out");
            Console.WriteLine("");

            for (i = 0; i <= 10; i++)
            {
                seed_in = seed;
                Niederreiter2.niederreiter2(DIM_NUM, ref seed, ref r);
                seed_out = seed;
                seed_out = seed;
                string cout = seed_in.ToString().PadLeft(6) + "  ";
                cout += seed_out.ToString().PadLeft(6) + " ";
                for (j = 0; j < DIM_NUM; j++)
                {
                    cout += r[j].ToString().PadLeft(10) + "  ";
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  Jump ahead by increasing SEED:");

            seed = 98;

            Console.WriteLine("");
            Console.WriteLine("  Seed  Seed   Niederreiter2");
            Console.WriteLine("  In    Out");
            Console.WriteLine("");

            for (i = 1; i <= 5; i++)
            {
                seed_in = seed;
                Niederreiter2.niederreiter2(DIM_NUM, ref seed, ref r);
                seed_out = seed;
                seed_out = seed;
                string cout = seed_in.ToString().PadLeft(6) + "  ";
                cout += seed_out.ToString().PadLeft(6) + " ";
                for (j = 0; j < DIM_NUM; j++)
                {
                    cout += r[j].ToString().PadLeft(10) + "  ";
                }

                Console.WriteLine(cout);
            }
        }
    }
}