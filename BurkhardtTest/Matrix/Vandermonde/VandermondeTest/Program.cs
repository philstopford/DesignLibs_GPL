using System;
using Burkardt;
using Burkardt.Types;
using Burkardt.Uniform;

namespace VandermondeTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            /*
            Purpose:
            
            MAIN is the main program for VANDERMONDE_TEST.
            
            Discussion:
            
            VANDERMONDE_TEST tests the VANDERMONDE library.
            
            Licensing:
            
            This code is distributed under the GNU LGPL license.
            
            Modified:
            
            01 May 2014
            
            Author:
            
            John Burkardt
            */
        {
            Console.WriteLine("");
            Console.WriteLine("VANDERMONDE_TEST");
            Console.WriteLine("  Test the VANDERMONDE library.");

            bivand1_test();
            bivand2_test();
            dvand_test();
            dvandprg_test();
            pvand_test();
            pvandprg_test();

            Console.WriteLine("");
            Console.WriteLine("VANDERMONDE_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void bivand1_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BIVAND1_TEST tests BIVAND1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            double[] a;
            double[] alpha = {1.0, 2.0, 3.0};
            double[] beta = {10.0, 20.0, 30.0};
            int n = N;
            int n2;

            Console.WriteLine("");
            Console.WriteLine("BIVAND1_TEST:");
            Console.WriteLine("  Compute a bidimensional Vandermonde matrix");
            Console.WriteLine("  associated with polynomials of");
            Console.WriteLine("  total degree less than N.");

            typeMethods.r8vec_print(n, alpha, "  Vandermonde vector ALPHA:");
            typeMethods.r8vec_print(n, beta, "  Vandermonde vector BETA:");

            a = VandermondeMatrix.bivand1(n, alpha, beta);

            n2 = (n * (n + 1)) / 2;
            typeMethods.r8mat_print(n2, n2, a, "  Bidimensional Vandermonde matrix:");

        }

        static void bivand2_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BIVAND2_TEST tests BIVAND2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 3;

            double[] a;
            double[] alpha = {1.0, 2.0, 3.0};
            double[] beta = {10.0, 20.0, 30.0};
            int n = N;
            int n2;

            Console.WriteLine("");
            Console.WriteLine("BIVAND2_TEST:");
            Console.WriteLine("  Compute a bidimensional Vandermonde matrix");
            Console.WriteLine("  associated with polynomials of maximum degree less than N.");

            typeMethods.r8vec_print(n, alpha, "  Vandermonde vector ALPHA:");
            typeMethods.r8vec_print(n, beta, "  Vandermonde vector BETA:");

            a = VandermondeMatrix.bivand2(n, alpha, beta);

            n2 = n * n;
            typeMethods.r8mat_print(n2, n2, a, "  Bidimensional Vandermonde matrix:");

        }

        static void dvand_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DVAND_TEST tests DVAND.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 5;

            double[] a;
            double[] alpha = null;
            double[] alpha1 = {0.0, 1.0, 2.0, 3.0, 4.0};
            double[] b;
            int n = N;
            int seed = 12345;
            int test;
            double[] x;
            double[] x1 = {5.0, 3.0, 4.0, 1.0, 2.0};

            Console.WriteLine("");
            Console.WriteLine("DVAND_TEST:");
            Console.WriteLine("  Solve a Vandermonde linear system A'*x=b");

            for (test = 1; test <= 2; test++)
            {
                if (test == 1)
                {
                    alpha = typeMethods.r8vec_copy_new(n, alpha1);
                }
                else if (test == 2)
                {
                    alpha = UniformRNG.r8vec_uniform_01_new(n, ref seed);
                }

                typeMethods.r8vec_print(n, alpha, "  Vandermonde vector ALPHA:");

                a = VandermondeMatrix.vand1(n, alpha);

                x = typeMethods.r8vec_copy_new(n, x1);
                b = typeMethods.r8mat_mtv_new(n, n, a, x);
                typeMethods.r8vec_print(n, b, "  Right hand side B:");

                x = VandermondeMatrix.dvand(n, alpha, b);
                typeMethods.r8vec_print(n, x, "  Solution X:");
            }
        }

        static void dvandprg_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DVANDPRG_TEST tests DVANDPRG.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 5;

            double[] a;
            double[] alpha = null;
            double[] alpha1 = {0.0, 1.0, 2.0, 3.0, 4.0};
            double[] b;
            double[] c;
            double[] m;
            int n = N;
            int nsub;
            int seed;
            int test;
            double[] x;
            double[] x1 = {5.0, 3.0, 4.0, 1.0, 2.0};

            Console.WriteLine("");
            Console.WriteLine("DVANDPRG_TEST:");
            Console.WriteLine("  Solve a Vandermonde linear system A'*x=b");
            Console.WriteLine("  progressively.");
            Console.WriteLine("  First we use ALPHA = 0, 1, 2, 3, 4.");
            Console.WriteLine("  Then we choose ALPHA at random.");

            for (test = 1; test <= 2; test++)
            {
                if (test == 1)
                {
                    alpha = typeMethods.r8vec_copy_new(n, alpha1);
                }
                else if (test == 2)
                {
                    seed = 123456789;
                    alpha = UniformRNG.r8vec_uniform_01_new(n, ref seed);
                }

                typeMethods.r8vec_print(n, alpha, "  Vandermonde vector ALPHA:");

                a = VandermondeMatrix.vand1(n, alpha);

                x = typeMethods.r8vec_copy_new(n, x1);
                b = typeMethods.r8mat_mtv_new(n, n, a, x);
                typeMethods.r8vec_print(n, b, "  Right hand side B:");

                x = new double[n];
                c = new double[n];
                m = new double[n];

                for (nsub = 1; nsub <= n; nsub++)
                {
                    VandermondeMatrix.dvandprg(nsub, alpha, b, ref x, ref c, ref m);
                    typeMethods.r8vec_print(nsub, x, "  Solution X:");
                }
            }
        }

        static void pvand_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PVAND_TEST tests PVAND.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 5;

            double[] a;
            double[] alpha = null;
            double[] alpha1 = {0.0, 1.0, 2.0, 3.0, 4.0};
            double[] b;
            int n = N;
            int seed;
            int test;
            double[] x;
            double[] x1 = {5.0, 3.0, 4.0, 1.0, 2.0};

            Console.WriteLine("");
            Console.WriteLine("PVAND_TEST:");
            Console.WriteLine("  Solve a Vandermonde linear system A*x=b");

            for (test = 1; test <= 2; test++)
            {
                if (test == 1)
                {
                    alpha = typeMethods.r8vec_copy_new(n, alpha1);
                }
                else if (test == 2)
                {
                    seed = 123456789;
                    alpha = UniformRNG.r8vec_uniform_01_new(n, ref seed);
                }

                typeMethods.r8vec_print(n, alpha, "  Vandermonde vector ALPHA:");

                a = VandermondeMatrix.vand1(n, alpha);

                x = typeMethods.r8vec_copy_new(n, x1);
                b = typeMethods.r8mat_mv_new(n, n, a, x);
                typeMethods.r8vec_print(n, b, "  Right hand side B:");

                x = VandermondeMatrix.pvand(n, alpha, b);
                typeMethods.r8vec_print(n, x, "  Solution X:");

            }
        }

        static void pvandprg_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PVANDPRG_TEST tests PVANDPRG.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 5;

            double[] a;
            double[] alpha = null;
            double[] alpha1 = {0.0, 1.0, 2.0, 3.0, 4.0};
            double[] b;
            double[] c;
            double[] m;
            int n = N;
            int nsub;
            int seed;
            int test;
            double[] x;
            double[] x1 = {5.0, 3.0, 4.0, 1.0, 2.0};

            Console.WriteLine("");
            Console.WriteLine("PVANDPRG_TEST:");
            Console.WriteLine("  Solve a Vandermonde linear system A*x=b");

            for (test = 1; test <= 2; test++)
            {
                if (test == 1)
                {
                    alpha = typeMethods.r8vec_copy_new(n, alpha1);
                }
                else if (test == 2)
                {
                    seed = 123456789;
                    alpha = UniformRNG.r8vec_uniform_01_new(n, ref seed);
                }

                typeMethods.r8vec_print(n, alpha, "  Vandermonde vector ALPHA:");

                a = VandermondeMatrix.vand1(n, alpha);

                x = typeMethods.r8vec_copy_new(n, x1);
                b = typeMethods.r8mat_mv_new(n, n, a, x);
                typeMethods.r8vec_print(n, b, "  Right hand side B:");

                x = new double[n];
                c = new double[n];
                m = new double[n];

                for (nsub = 1; nsub <= n; nsub++)
                {
                    VandermondeMatrix.pvandprg(nsub, alpha, b, ref x, ref c, ref m);
                    typeMethods.r8vec_print(nsub, x, "  Solution X:");
                }
            }

        }
    }
}