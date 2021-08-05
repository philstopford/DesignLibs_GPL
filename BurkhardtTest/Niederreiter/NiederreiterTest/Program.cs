using System;
using Burkardt.NiederreiterNS;

namespace NiederreiterTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for NIEDERREITER_TEST.
            //
            //  Discussion:
            //
            //    NIEDERREITER_TEST tests the NIEDERREITER library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 June 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int base_;
            int dim_num;

            Console.WriteLine("");
            Console.WriteLine("NIEDERREITER_TEST");
            Console.WriteLine("  Test the NIEDERREITER routines.");

            base_ = 2;
            test01(base_);

            base_ = 3;
            test01(base_);

            base_ = 13;
            test01(base_);

            base_ = 2;
            test02(base_);

            base_ = 3;
            test02(base_);

            base_ = 2;
            dim_num = 20;
            test03(base_, dim_num);

            base_ = 2;
            dim_num = 29;
            test03(base_, dim_num);
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("NIEDERREITER_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01(int base_)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests NIEDERREITER.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 September 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int BASE, the base_ to use in the computation.
            //    BASE should be a prime, or a power of a prime.
            //
        {
            const int dim_max = 4;

            int dim;
            int dim_num;
            int i;
            double[] r = new double[dim_max];
            int seed;
            int seed_in;
            int seed_out;
            Niederreiter.NiederReiterData data = new Niederreiter.NiederReiterData();

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  NIEDERREITER computes the next element of ");
            Console.WriteLine("  a Niederreiter quasirandom sequence using base_ BASE.");
            Console.WriteLine("");
            Console.WriteLine("  In this test, we call NIEDERREITER repeatedly.");
            Console.WriteLine("");
            Console.WriteLine("  Using base_ BASE =      " + base_ + "");

            for (dim_num = 2; dim_num <= dim_max; dim_num++)
            {
                seed = 0;

                Console.WriteLine("");
                Console.WriteLine("  Using dimension DIM_NUM =   " + dim_num + "");
                Console.WriteLine("");
                Console.WriteLine("    Seed    Seed     Niederreiter");
                Console.WriteLine("      In     Out");
                Console.WriteLine("");
                for (i = 0; i <= 110; i++)
                {
                    seed_in = seed;
                    Niederreiter.niederreiter(ref data, dim_num, base_, ref seed, ref r);
                    seed_out = seed;
                    if (i <= 11 || 95 <= i)
                    {
                        string cout = "  " + seed_in.ToString().PadLeft(8)
                                          + "  " + seed_out.ToString().PadLeft(8);
                        for (dim = 0; dim < dim_num; dim++)
                        {
                            cout += r[dim].ToString().PadLeft(10);
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

        static void test02(int base_)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests NIEDERREITER.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 September 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int BASE, the base_ to use in the computation.
            //    BASE should be a prime, or a power of a prime.
            //
        {
            const int dim_num = 3;

            int dim;
            int i;
            double[] r = new double[dim_num];
            int seed;
            int seed_in;
            int seed_out;
            string cout = "";
            Niederreiter.NiederReiterData data = new Niederreiter.NiederReiterData();

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  NIEDERREITER computes the next element of");
            Console.WriteLine("  a Niederreiter quasirandom sequence using base_ BASE.");
            Console.WriteLine("");
            Console.WriteLine("  In this test, we demonstrate how the SEED can be");
            Console.WriteLine("  manipulated to skip ahead in the sequence, or");
            Console.WriteLine("  to come back to any part of the sequence.");

            Console.WriteLine("");
            Console.WriteLine("  Using base_ BASE =           " + base_ + "");
            Console.WriteLine("  Using dimension DIM_NUM =   " + dim_num + "");

            seed = 0;

            Console.WriteLine("");
            Console.WriteLine("    Seed    Seed     Niederreiter");
            Console.WriteLine("      In     Out");
            Console.WriteLine("");
            for (i = 0; i <= 10; i++)
            {
                seed_in = seed;
                Niederreiter.niederreiter(ref data, dim_num, base_, ref seed, ref r);
                seed_out = seed;
                cout = "  " + seed_in.ToString().PadLeft(8)
                                  + "  " + seed_out.ToString().PadLeft(8);
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += r[dim].ToString().PadLeft(10);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  Jump ahead by increasing SEED:");
            Console.WriteLine("");

            seed = 100;

            Console.WriteLine("");
            Console.WriteLine("    Seed    Seed     Niederreiter");
            Console.WriteLine("      In     Out");
            Console.WriteLine("");
            for (i = 1; i <= 5; i++)
            {
                seed_in = seed;
                Niederreiter.niederreiter(ref data, dim_num, base_, ref seed, ref r);
                seed_out = seed;
                cout = "  " + seed_in.ToString().PadLeft(8)
                                  + "  " + seed_out.ToString().PadLeft(8);
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += r[dim].ToString().PadLeft(10);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  Jump back by decreasing SEED:");
            Console.WriteLine("");

            seed = 3;

            Console.WriteLine("");
            Console.WriteLine("    Seed    Seed     Niederreiter");
            Console.WriteLine("      In     Out");
            Console.WriteLine("");
            for (i = 0; i <= 10; i++)
            {
                seed_in = seed;
                Niederreiter.niederreiter(ref data, dim_num, base_, ref seed, ref r);
                seed_out = seed;
                cout = "  " + seed_in.ToString().PadLeft(8)
                                  + "  " + seed_out.ToString().PadLeft(8);
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += r[dim].ToString().PadLeft(10);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  Jump ahead by increasing SEED:");
            Console.WriteLine("");

            seed = 98;

            Console.WriteLine("");
            Console.WriteLine("    Seed    Seed     Niederreiter");
            Console.WriteLine("      In     Out");
            Console.WriteLine("");
            cout = "";
            for (i = 1; i <= 5; i++)
            {
                seed_in = seed;
                Niederreiter.niederreiter(ref data, dim_num, base_, ref seed, ref r);
                seed_out = seed;
                cout = "  " + seed_in.ToString().PadLeft(8)
                                  + "  " + seed_out.ToString().PadLeft(8);
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += r[dim].ToString().PadLeft(10);
                }

                Console.WriteLine(cout);
            }

        }

        static void test03(int base_, int dim_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 tests NIEDERREITER.
            //
            //  Discussion:
            //
            //    Simply verify that a few terms of a sequence of given dimension
            //    can be computed.  Most recently, the NIEDERREITER code was set
            //    up to handle up to dimension 50...we think.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 June 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int BASE, the base_ to use in the computation.
            //    BASE should be a prime, or a power of a prime.
            //
            //    Input, int DIM, the spatial dimension.
            //
        {
            int dim;
            int i;
            double[] r;
            int seed;
            int seed_in;
            int seed_out;
            Niederreiter.NiederReiterData data = new Niederreiter.NiederReiterData();

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  NIEDERREITER computes the next element of");
            Console.WriteLine("  a Niederreiter quasirandom sequence using base_ BASE.");
            Console.WriteLine("");
            Console.WriteLine("  In this test, we simply generate ten elements in a given base_");
            Console.WriteLine("  and dimension.");
            Console.WriteLine("  manipulated to skip ahead in the sequence, or");
            Console.WriteLine("  to come back to any part of the sequence.");

            Console.WriteLine("");
            Console.WriteLine("  Using base_ BASE =           " + base_ + "");
            Console.WriteLine("  Using dimension DIM_NUM =   " + dim_num + "");

            seed = 0;
            r = new double[dim_num];

            Console.WriteLine("");
            Console.WriteLine("    Seed    Seed     Niederreiter");
            Console.WriteLine("      In     Out");
            Console.WriteLine("");
            string cout = "";
            for (i = 0; i <= 10; i++)
            {
                seed_in = seed;
                Niederreiter.niederreiter(ref data, dim_num, base_, ref seed, ref r);
                seed_out = seed;
                cout = "  " + seed_in.ToString().PadLeft(8)
                                  + "  " + seed_out.ToString().PadLeft(8);
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += r[dim].ToString().PadLeft(10);
                    if (((dim + 1) % 5 == 0) && (dim != dim_num))
                    {
                        Console.WriteLine(cout);
                        cout = "                    ";
                    }
                }

                Console.WriteLine(cout);
            }
        }
    }
}