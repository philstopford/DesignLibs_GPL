using System;
using Burkardt;
using Burkardt.PolynomialNS;
using Burkardt.RandomNS;

namespace FaureTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for FAURE_TEST.
            //
            //  Discussion:
            //
            //    FAURE_TEST tests the FAURE library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("FAURE_TEST");
            Console.WriteLine("  Test the FAURE library.");

            test005();
            test006();
            test01();
            test02();
            test03();

            Console.WriteLine("");
            Console.WriteLine("FAURE_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test005()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST005 tests BINOMIAL_TABLE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] coef;
            int i;
            int j;
            int m = 10;
            int n = 7;
            int qs = 7;

            Console.WriteLine("");
            Console.WriteLine("TEST005");
            Console.WriteLine("  BINOMIAL_TABLE computes a table of binomial.");
            Console.WriteLine("  coefficients mod QS.");
            Console.WriteLine("");
            Console.WriteLine("  Here, QS = " + qs + "");

            coef = Coefficients.binomial_table(qs, m, n);

            Console.WriteLine("");
            string cout = "   I/J";
            for (j = 0; j <= n; j++)
            {
                cout += j.ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
            Console.WriteLine("");

            for (i = 0; i <= m; i++)
            {
                cout = "  " + i.ToString().PadLeft(2) + "  ";
                for (j = 0; j <= n; j++)
                {
                    cout += coef[i + j * (m + 1)].ToString().PadLeft(8);
                }

                Console.WriteLine(cout);
            }
        }

        static void test006()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST006 tests I4_LOG_I4.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i4;
            int j4;

            Console.WriteLine("");
            Console.WriteLine("TEST006");
            Console.WriteLine("  I4_LOG_R8: whole part of log base B,");
            Console.WriteLine("");
            Console.WriteLine("        I4        J4 I4_LOG_J4");
            Console.WriteLine("");

            for (j4 = 2; j4 <= 5; j4++)
            {
                for (i4 = 0; i4 <= 10; i4++)
                {
                    Console.WriteLine("  " + i4.ToString().PadLeft(8)
                        + "  " + j4.ToString().PadLeft(8)
                        + "  " + Math.Log(i4, j4).ToString().PadLeft(8) + "");
                }

                Console.WriteLine("");
            }
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests FAURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_MAX = 4;

            int dim;
            int dim_num;
            int i;
            int qs;
            double[] r = new double[DIM_MAX];
            int seed;
            int seed_in;
            int seed_out;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  FAURE computes the next element of a Faure sequence.");
            Console.WriteLine("");
            Console.WriteLine("  In this test, we call FAURE repeatedly.");

            FaureData data = new FaureData();
            
            for (dim_num = 2; dim_num <= DIM_MAX; dim_num++)
            {

                seed = -1;
                qs = Prime.prime_ge(dim_num);

                Console.WriteLine("");
                Console.WriteLine("  Using dimension DIM_NUM =   " + dim_num + "");
                Console.WriteLine("  The underlying base is QS = " + qs + "");
                Console.WriteLine("");
                Console.WriteLine("  Seed  Seed   Faure");
                Console.WriteLine("  In    Out");
                Console.WriteLine("");

                for (i = 1; i <= 10; i++)
                {
                    seed_in = seed;
                    Faure.faure(ref data, dim_num, ref seed, ref r);
                    seed_out = seed;
                    string cout = seed_in.ToString().PadLeft(6) + "  "
                        + seed_out.ToString().PadLeft(6) + "  ";
                    for (dim = 0; dim < dim_num; dim++)
                    {
                        cout += r[dim].ToString().PadLeft(10) + "  ";
                    }

                    Console.WriteLine(cout);
                }

            }
        }

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests FAURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;

            int dim;
            int i;
            int qs;
            double[] r = new double[DIM_NUM];
            int seed;
            int seed_in;
            int seed_out;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  FAURE computes the next element of a Faure sequence.");
            Console.WriteLine("");
            Console.WriteLine("  In this test, we demonstrate how the SEED can be");
            Console.WriteLine("  manipulated to skip ahead in the sequence, or");
            Console.WriteLine("  to come back to any part of the sequence.");

            qs = Prime.prime_ge(DIM_NUM);

            Console.WriteLine("");
            Console.WriteLine("  Using dimension DIM_NUM =   " + DIM_NUM + "");
            Console.WriteLine("  The underlying base is QS = " + qs + "");

            Console.WriteLine("");
            Console.WriteLine("  Note that on the first call to FAURE, if");
            Console.WriteLine("  SEED is negative, it is reset to a value that");
            Console.WriteLine("  is the recommended starting point:");

            seed = -1;
            Console.WriteLine("");
            Console.WriteLine("  Seed  Seed   Faure");
            Console.WriteLine("  In    Out");
            Console.WriteLine("");

            FaureData data = new FaureData();
            
            for (i = 1; i <= 5; i++)
            {
                seed_in = seed;
                Faure.faure(ref data, DIM_NUM, ref seed, ref r);
                seed_out = seed;
                string cout = seed_in.ToString().PadLeft(6) + "  "
                    + seed_out.ToString().PadLeft(6) + "  ";
                for (dim = 0; dim < DIM_NUM; dim++)
                {
                    cout += r[dim].ToString().PadLeft(10) + "  ";
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  However, if the input value of SEED is 0,");
            Console.WriteLine("  then no initial skipping is done.");

            seed = 0;

            Console.WriteLine("");
            Console.WriteLine("  Seed  Seed   Faure");
            Console.WriteLine("  In    Out");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                seed_in = seed;
                Faure.faure(ref data, DIM_NUM, ref seed, ref r);
                seed_out = seed;
                string cout = seed_in.ToString().PadLeft(6) + "  "
                    + seed_out.ToString().PadLeft(6) + "  ";
                for (dim = 0; dim < DIM_NUM; dim++)
                {
                    cout += r[dim].ToString().PadLeft(10) + "  ";
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  Jump ahead by increasing SEED:");
            Console.WriteLine("");

            seed = 100;

            Console.WriteLine("");
            Console.WriteLine("  Seed  Seed   Faure");
            Console.WriteLine("  In    Out");
            Console.WriteLine("");

            for (i = 1; i <= 5; i++)
            {
                seed_in = seed;
                Faure.faure(ref data, DIM_NUM, ref seed, ref r);
                seed_out = seed;
                string cout= seed_in.ToString().PadLeft(6) + "  "
                    + seed_out.ToString().PadLeft(6) + "  ";
                for (dim = 0; dim < DIM_NUM; dim++)
                {
                    cout += r[dim].ToString().PadLeft(10) + "  ";
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  Jump back by decreasing SEED:");
            Console.WriteLine("");

            seed = 3;

            Console.WriteLine("");
            Console.WriteLine("  Seed  Seed   Faure");
            Console.WriteLine("  In    Out");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                seed_in = seed;
                Faure.faure(ref data, DIM_NUM, ref seed, ref r);
                seed_out = seed;
                string cout = seed_in.ToString().PadLeft(6) + "  "
                    + seed_out.ToString().PadLeft(6) + "  ";
                for (dim = 0; dim < DIM_NUM; dim++)
                {
                    cout += r[dim].ToString().PadLeft(10) + "  ";
                }

                Console.WriteLine(cout);
            }
        }

        static void test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 tests FAURE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int dim;
            int dim_base = 10;
            int dim_num;
            int i;
            int qs;
            double[] r;
            int seed;
            int seed_in;
            int seed_out;

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  FAURE computes the next element of a Faure sequence.");
            Console.WriteLine("");
            Console.WriteLine("  In this test, we try some large dimensions.");

            FaureData data = new FaureData();
            
            for (dim_num = dim_base; dim_num <= 6 * dim_base; dim_num = dim_num + dim_base)
            {
                r = new double[dim_num];

                seed = -1;
                qs = Prime.prime_ge(dim_num);

                Console.WriteLine("");
                Console.WriteLine("  Using dimension DIM_NUM =   " + dim_num + "");
                Console.WriteLine("  The underlying base is QS = " + qs + "");
                Console.WriteLine("");
                Console.WriteLine("  Seed  Seed   Faure");
                Console.WriteLine("  In    Out");
                Console.WriteLine("");

                for (i = 1; i <= 2; i++)
                {
                    seed_in = seed;
                    Faure.faure(ref data, dim_num, ref seed, ref r);
                    seed_out = seed;
                    Console.WriteLine("  " + seed_in.ToString().PadLeft(8)
                        + "  " + seed_out.ToString().PadLeft(8) + "");
                    string cout = "                    ";
                    for (dim = 0; dim < dim_num; dim++)
                    {
                        cout += "  " + r[dim].ToString().PadLeft(10);
                        if ((dim + 1) % 5 == 0 || dim + 1 == dim_num)
                        {
                            Console.WriteLine(cout);
                        }

                        if ((dim + 1) % 5 == 0 && dim + 1 < dim_num)
                        {
                            cout += "                    ";
                        }
                    }
                }
            }
        }
    }
}