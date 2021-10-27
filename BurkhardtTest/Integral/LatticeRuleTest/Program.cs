using System;
using Burkardt;
using Burkardt.Function;
using Burkardt.IntegralNS;
using Burkardt.Types;

namespace LatticeRuleTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for LATTICE_RULE_TEST.
            //
            //  Discussion:
            //
            //    LATTICE_RULE_TEST tests the LATTICE_RULE library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 August 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("LATTICE_RULE_TEST");
            Console.WriteLine("  Test the LATTICE_RULE library.");

            test01();
            test02();
            test03();
            test04();
            test05();
            test06();
            test07();
            test08();
            test085();
            test09();

            test10();
            test11();
            test12();
            test13();
            test14();

            Console.WriteLine("");
            Console.WriteLine("LATTICE_RULE_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");

        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests FIBONACCI_LATTICE_Q.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 November 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] a;
            double[] b;
            int dim;
            int dim_num = 2;
            double error;
            double exact;
            int k;
            int m;
            double quad;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  FIBONACCI_LATTICE_Q applies a Fibonacci lattice rule");
            Console.WriteLine("  to integrate a function over the unit square.");
            Console.WriteLine("  These Fibonacci rules are only available in 2D.");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM = " + dim_num + "");

            a = new double[dim_num];
            b = new double[dim_num];

            for (dim = 0; dim < dim_num; dim++)
            {
                a[dim] = 0.0;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                b[dim] = 1.0;
            }

            exact = FibonacciLattice.e_01_2d(dim_num, a, b);

            Console.WriteLine("");
            Console.WriteLine("         K         M      EXACT     ESTIMATE  ERROR");
            Console.WriteLine("");

            for (k = 3; k <= 18; k++)
            {
                m = Helpers.fibonacci(k);

                quad = FibonacciLattice.fibonacci_lattice_q(k, FibonacciLattice.f_01_2d);

                error = Math.Abs(exact - quad);

                Console.WriteLine("  " + k.ToString().PadLeft(8)
                                       + "  " + m.ToString().PadLeft(8)
                                       + "  " + exact.ToString().PadLeft(10)
                                       + "  " + quad.ToString().PadLeft(10)
                                       + "  " + error.ToString().PadLeft(10) + "");
            }

        }

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests FIBONACCI_LATTICE_T.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 November 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] a;
            double[] b;
            int dim;
            int dim_num = 2;
            double error;
            double exact;
            int k;
            int m;
            double quad;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  FIBONACCI_LATTICE_T applies a symmetric Fibonacci lattice rule");
            Console.WriteLine("  to integrate a function over the unit square.");
            Console.WriteLine("  These Fibonacci rules are only available in 2D.");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM = " + dim_num + "");

            a = new double[dim_num];
            b = new double[dim_num];

            for (dim = 0; dim < dim_num; dim++)
            {
                a[dim] = 0.0;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                b[dim] = 1.0;
            }

            exact = FibonacciLattice.e_01_2d(dim_num, a, b);

            Console.WriteLine("");
            Console.WriteLine("         K         M      EXACT     ESTIMATE  ERROR");
            Console.WriteLine("");

            for (k = 3; k <= 18; k++)
            {
                m = Helpers.fibonacci(k);

                quad = FibonacciLattice.fibonacci_lattice_t(k, FibonacciLattice.f_01_2d);

                error = Math.Abs(exact - quad);

                Console.WriteLine("  " + k.ToString().PadLeft(8)
                                       + "  " + m.ToString().PadLeft(8)
                                       + "  " + exact.ToString().PadLeft(10)
                                       + "  " + quad.ToString().PadLeft(10)
                                       + "  " + error.ToString().PadLeft(10) + "");
            }

        }

        static void test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 tests FIBONACCI_LATTICE_B.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 November 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] a;
            double[] b;
            int dim;
            int dim_num = 2;
            double error;
            double exact;
            int k;
            int m;
            double quad;

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  FIBONACCI_LATTICE_B applies an optimal Fibonacci lattice rule");
            Console.WriteLine("  to integrate a function over the unit square.");
            Console.WriteLine("  These Fibonacci rules are only available in 2D.");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM = " + dim_num + "");

            a = new double[dim_num];
            b = new double[dim_num];

            for (dim = 0; dim < dim_num; dim++)
            {
                a[dim] = 0.0;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                b[dim] = 1.0;
            }

            exact = FibonacciLattice.e_01_2d(dim_num, a, b);

            Console.WriteLine("");
            Console.WriteLine("         K         M      EXACT     ESTIMATE  ERROR");
            Console.WriteLine("");

            for (k = 3; k <= 18; k++)
            {
                m = Helpers.fibonacci(k);

                quad = FibonacciLattice.fibonacci_lattice_b(k, FibonacciLattice.f_01_2d);

                error = Math.Abs(exact - quad);

                Console.WriteLine("  " + k.ToString().PadLeft(8)
                                       + "  " + m.ToString().PadLeft(8)
                                       + "  " + exact.ToString().PadLeft(10)
                                       + "  " + quad.ToString().PadLeft(10)
                                       + "  " + error.ToString().PadLeft(10) + "");
            }

        }

        static void test04()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04 tests FIBONACCI_LATTICE_Q1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 November 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] a;
            double[] b;
            int dim;
            int dim_num = 2;
            double error;
            double exact;
            int k;
            int m;
            double quad;

            Console.WriteLine("");
            Console.WriteLine("TEST04");
            Console.WriteLine("  FIBONACCI_LATTICE_Q1 applies a Fibonacci lattice rule");
            Console.WriteLine("  to integrate a function over the unit square.");
            Console.WriteLine("  A nonlinear coordinate transformation is applied.");
            Console.WriteLine("  These Fibonacci rules are only available in 2D.");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM = " + dim_num + "");

            a = new double[dim_num];
            b = new double[dim_num];

            for (dim = 0; dim < dim_num; dim++)
            {
                a[dim] = 0.0;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                b[dim] = 1.0;
            }

            exact = FibonacciLattice.e_01_2d(dim_num, a, b);

            Console.WriteLine("");
            Console.WriteLine("         K         M      EXACT     ESTIMATE  ERROR");
            Console.WriteLine("");

            for (k = 3; k <= 18; k++)
            {
                m = Helpers.fibonacci(k);

                quad = FibonacciLattice.fibonacci_lattice_q1(k, FibonacciLattice.f_01_2d);

                error = Math.Abs(exact - quad);

                Console.WriteLine("  " + k.ToString().PadLeft(8)
                                       + "  " + m.ToString().PadLeft(8)
                                       + "  " + exact.ToString().PadLeft(10)
                                       + "  " + quad.ToString().PadLeft(10)
                                       + "  " + error.ToString().PadLeft(10) + "");
            }

        }

        static void test05()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST05 tests FIBONACCI_LATTICE_Q2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 November 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] a;
            double[] b;
            int dim;
            int dim_num = 2;
            double error;
            double exact;
            int k;
            int m;
            double quad;

            Console.WriteLine("");
            Console.WriteLine("TEST05");
            Console.WriteLine("  FIBONACCI_LATTICE_Q2 applies a Fibonacci lattice rule");
            Console.WriteLine("  to integrate a function over the unit square.");
            Console.WriteLine("  A nonlinear coordinate transformation is applied.");
            Console.WriteLine("  These Fibonacci rules are only available in 2D.");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM = " + dim_num + "");

            a = new double[dim_num];
            b = new double[dim_num];

            for (dim = 0; dim < dim_num; dim++)
            {
                a[dim] = 0.0;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                b[dim] = 1.0;
            }

            exact = FibonacciLattice.e_01_2d(dim_num, a, b);

            Console.WriteLine("");
            Console.WriteLine("         K         M      EXACT     ESTIMATE  ERROR");
            Console.WriteLine("");

            for (k = 3; k <= 18; k++)
            {
                m = Helpers.fibonacci(k);

                quad = FibonacciLattice.fibonacci_lattice_q2(k, FibonacciLattice.f_01_2d);

                error = Math.Abs(exact - quad);

                Console.WriteLine("  " + k.ToString().PadLeft(8)
                                       + "  " + m.ToString().PadLeft(8)
                                       + "  " + exact.ToString().PadLeft(10)
                                       + "  " + quad.ToString().PadLeft(10)
                                       + "  " + error.ToString().PadLeft(10) + "");
            }

        }

        static void test06()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST06 tests FIBONACCI_LATTICE_Q3.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 November 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] a;
            double[] b;
            int dim;
            int dim_num = 2;
            double error;
            double exact;
            int k;
            int m;
            double quad;

            Console.WriteLine("");
            Console.WriteLine("TEST06");
            Console.WriteLine("  FIBONACCI_LATTICE_Q3 applies a Fibonacci lattice rule");
            Console.WriteLine("  to integrate a function over the unit square.");
            Console.WriteLine("  A nonlinear coordinate transformation is applied.");
            Console.WriteLine("  These Fibonacci rules are only available in 2D.");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM = " + dim_num + "");

            a = new double[dim_num];
            b = new double[dim_num];

            for (dim = 0; dim < dim_num; dim++)
            {
                a[dim] = 0.0;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                b[dim] = 1.0;
            }

            exact = FibonacciLattice.e_01_2d(dim_num, a, b);

            Console.WriteLine("");
            Console.WriteLine("         K         M      EXACT     ESTIMATE  ERROR");
            Console.WriteLine("");

            for (k = 3; k <= 18; k++)
            {
                m = Helpers.fibonacci(k);

                quad = FibonacciLattice.fibonacci_lattice_q3(k, FibonacciLattice.f_01_2d);

                error = Math.Abs(exact - quad);

                Console.WriteLine("  " + k.ToString().PadLeft(8)
                                       + "  " + m.ToString().PadLeft(8)
                                       + "  " + exact.ToString().PadLeft(10)
                                       + "  " + quad.ToString().PadLeft(10)
                                       + "  " + error.ToString().PadLeft(10) + "");
            }

        }

        static void test07()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST07 tests LATTICE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 November 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Ian Sloan, Stephen Joe,
            //    Lattice Methods for Multiple Integration,
            //    Oxford, 1994, page 18.
            //
        {
            double[] a;
            double[] b;
            int dim;
            int dim_num = 2;
            double error;
            double exact;
            int i;
            int m;
            double quad;
            int[] z;
            ;

            Console.WriteLine("");
            Console.WriteLine("TEST07");
            Console.WriteLine("  LATTICE applies a lattice rule to integrate");
            Console.WriteLine("  a function over the unit hypercube.");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM = " + dim_num + "");
            Console.WriteLine("  The lattice rule order M will vary.");

            z = new int[dim_num];

            z[0] = 1;
            z[1] = 2;

            a = new double[dim_num];
            b = new double[dim_num];

            for (dim = 0; dim < dim_num; dim++)
            {
                a[dim] = 0.0;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                b[dim] = 1.0;
            }

            typeMethods.i4vec_print(dim_num, z, "  The lattice generator vector:");

            Console.WriteLine("");
            Console.WriteLine("         I         M      EXACT     ESTIMATE  ERROR");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                m = Prime.prime(3 * i);

                quad = Lattice.lattice(dim_num, m, z, FibonacciLattice.f_01_2d);

                exact = FibonacciLattice.e_01_2d(dim_num, a, b);

                error = Math.Abs(exact - quad);

                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                       + "  " + m.ToString().PadLeft(8)
                                       + "  " + exact.ToString().PadLeft(10)
                                       + "  " + quad.ToString().PadLeft(10)
                                       + "  " + error.ToString().PadLeft(10) + "");
            }

        }

        static void test08()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST08 tests LATTICE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 November 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Ian Sloan, Stephen Joe,
            //    Lattice Methods for Multiple Integration,
            //    Oxford, 1994, page 18.
            //
        {
            double[] a;
            double[] b;
            int dim;
            int dim_num = 2;
            double error;
            double exact;
            int i;
            int m = 53;
            double quad;
            int[] z;
            ;

            Console.WriteLine("");
            Console.WriteLine("TEST08");
            Console.WriteLine("  LATTICE applies a lattice rule to integrate");
            Console.WriteLine("  a function over the unit hypercube.");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM = " + dim_num + "");
            Console.WriteLine("  The lattice rule order M will vary.");
            Console.WriteLine("  The lattice generator vector Z will vary.");

            z = new int[dim_num];

            z[0] = 1;

            a = new double[dim_num];
            b = new double[dim_num];

            for (dim = 0; dim < dim_num; dim++)
            {
                a[dim] = 0.0;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                b[dim] = 1.0;
            }

            typeMethods.i4vec_print(dim_num, z, "  The lattice generator vector:");

            Console.WriteLine("");
            Console.WriteLine("         M      Z[0]      Z[1]      EXACT     ESTIMATE  ERROR");
            Console.WriteLine("");

            for (i = 1; i <= m - 1; i++)
            {
                z[1] = i;

                quad = Lattice.lattice(dim_num, m, z, FibonacciLattice.f_01_2d);

                exact = FibonacciLattice.e_01_2d(dim_num, a, b);

                error = Math.Abs(exact - quad);

                Console.WriteLine("  " + m.ToString().PadLeft(8)
                                       + "  " + z[0].ToString().PadLeft(8)
                                       + "  " + z[1].ToString().PadLeft(8)
                                       + "  " + exact.ToString().PadLeft(10)
                                       + "  " + quad.ToString().PadLeft(10)
                                       + "  " + error.ToString().PadLeft(10) + "");
            }

        }

        static void test085()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST085 tests LATTICE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 November 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Ian Sloan, Stephen Joe,
            //    Lattice Methods for Multiple Integration,
            //    Oxford, 1994, page 18.
            //
        {
            double[] a;
            double[] b;
            int dim;
            int dim_num = 2;
            double error;
            double exact;
            int k;
            int m;
            double quad;
            int[] z;
            ;

            Console.WriteLine("");
            Console.WriteLine("TEST085");
            Console.WriteLine("  LATTICE is a lattice rule for periodic functions.");
            Console.WriteLine("  However, we apply it to a nonperiodic function");
            Console.WriteLine("  just to see how it does.");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM = " + dim_num + "");

            z = new int[dim_num];

            z[0] = 1;
            z[1] = 2;

            a = new double[dim_num];
            b = new double[dim_num];

            for (dim = 0; dim < dim_num; dim++)
            {
                a[dim] = 0.0;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                b[dim] = 1.0;
            }

            typeMethods.i4vec_print(dim_num, z, "  The lattice generator vector:");

            Console.WriteLine("");
            Console.WriteLine("         I         M      EXACT     ESTIMATE  ERROR");
            Console.WriteLine("");

            for (k = 3; k <= 18; k++)
            {
                m = Helpers.fibonacci(k);

                quad = Lattice.lattice(dim_num, m, z, FibonacciLattice.f_01_2d);

                exact = FibonacciLattice.e_01_2d(dim_num, a, b);

                error = Math.Abs(exact - quad);

                Console.WriteLine("  " + k.ToString().PadLeft(8)
                                       + "  " + m.ToString().PadLeft(8)
                                       + "  " + exact.ToString().PadLeft(10)
                                       + "  " + quad.ToString().PadLeft(10)
                                       + "  " + error.ToString().PadLeft(10) + "");
            }

        }

        static void test09()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST09 tests LATTICE_NP0.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 November 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Ian Sloan, Stephen Joe,
            //    Lattice Methods for Multiple Integration,
            //    Oxford, 1994, page 18.
            //
        {
            double[] a;
            double[] b;
            int dim;
            int dim_num = 2;
            double error;
            double exact;
            int k;
            int m;
            double quad;
            int[] z;
            ;

            Console.WriteLine("");
            Console.WriteLine("TEST09");
            Console.WriteLine("  LATTICE_NP0 applies a lattice rule to a");
            Console.WriteLine("  nonperiodic function by reflecting the function");
            Console.WriteLine("  about the midpoint and averaging.");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM = " + dim_num + "");

            z = new int[dim_num];

            z[0] = 1;
            z[1] = 2;

            a = new double[dim_num];
            b = new double[dim_num];

            for (dim = 0; dim < dim_num; dim++)
            {
                a[dim] = 0.0;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                b[dim] = 1.0;
            }

            typeMethods.i4vec_print(dim_num, z, "  The lattice generator vector:");

            Console.WriteLine("");
            Console.WriteLine("         I         M      EXACT     ESTIMATE  ERROR");
            Console.WriteLine("");

            for (k = 3; k <= 18; k++)
            {
                m = Helpers.fibonacci(k);

                quad = Lattice.lattice_np0(dim_num, m, z, FibonacciLattice.f_01_2d);

                exact = FibonacciLattice.e_01_2d(dim_num, a, b);

                error = Math.Abs(exact - quad);

                Console.WriteLine("  " + k.ToString().PadLeft(8)
                                       + "  " + m.ToString().PadLeft(8)
                                       + "  " + exact.ToString().PadLeft(10)
                                       + "  " + quad.ToString().PadLeft(10)
                                       + "  " + error.ToString().PadLeft(10) + "");
            }

        }

        static void test10()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST10 tests LATTICE_NP1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 November 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Ian Sloan, Stephen Joe,
            //    Lattice Methods for Multiple Integration,
            //    Oxford, 1994, page 18.
            //
        {
            double[] a;
            double[] b;
            int dim;
            int dim_num = 2;
            double error;
            double exact;
            int k;
            int m;
            double quad;
            int[] z;
            ;

            Console.WriteLine("");
            Console.WriteLine("TEST10");
            Console.WriteLine("  LATTICE_NP1 applies a lattice rule to a");
            Console.WriteLine("  nonperiodic function using a nonlinear transformation");
            Console.WriteLine("  to integrate a function over the unit square.");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM = " + dim_num + "");

            z = new int[dim_num];

            z[0] = 1;
            z[1] = 2;

            a = new double[dim_num];
            b = new double[dim_num];

            for (dim = 0; dim < dim_num; dim++)
            {
                a[dim] = 0.0;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                b[dim] = 1.0;
            }

            typeMethods.i4vec_print(dim_num, z, "  The lattice generator vector:");

            Console.WriteLine("");
            Console.WriteLine("         I         M      EXACT     ESTIMATE  ERROR");
            Console.WriteLine("");

            for (k = 3; k <= 18; k++)
            {
                m = Helpers.fibonacci(k);

                quad = Lattice.lattice_np1(dim_num, m, z, FibonacciLattice.f_01_2d);

                exact = FibonacciLattice.e_01_2d(dim_num, a, b);

                error = Math.Abs(exact - quad);

                Console.WriteLine("  " + k.ToString().PadLeft(8)
                                       + "  " + m.ToString().PadLeft(8)
                                       + "  " + exact.ToString().PadLeft(10)
                                       + "  " + quad.ToString().PadLeft(10)
                                       + "  " + error.ToString().PadLeft(10) + "");
            }

        }

        static void test11()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST11 tests MONTE_CARLO.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 November 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] a;
            double[] b;
            int dim;
            int dim_num = 2;
            double error;
            double exact;
            int k;
            int m;
            double quad;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("TEST11");
            Console.WriteLine("  MONTE_CARLO applies a Monte Carlo scheme");
            Console.WriteLine("  to estimate the integral of a function");
            Console.WriteLine("  over the unit hypercube.");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM = " + dim_num + "");

            a = new double[dim_num];
            b = new double[dim_num];

            for (dim = 0; dim < dim_num; dim++)
            {
                a[dim] = 0.0;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                b[dim] = 1.0;
            }

            seed = 123456789;

            exact = FibonacciLattice.e_01_2d(dim_num, a, b);

            Console.WriteLine("");
            Console.WriteLine("         K         M      EXACT     ESTIMATE  ERROR");
            Console.WriteLine("");

            for (k = 2; k <= 5; k++)
            {
                m = (int) Math.Pow(10, k);

                quad = MonteCarlo.monte_carlo(dim_num, m, FibonacciLattice.f_01_2d, ref seed);

                error = Math.Abs(exact - quad);

                Console.WriteLine("  " + k.ToString().PadLeft(8)
                                       + "  " + m.ToString().PadLeft(8)
                                       + "  " + exact.ToString().PadLeft(10)
                                       + "  " + quad.ToString().PadLeft(10)
                                       + "  " + error.ToString().PadLeft(10) + "");
            }

        }

        static void test12()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST12 tests Lattice.lattice_print.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 November 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Ian Sloan, Stephen Joe,
            //    Lattice Methods for Multiple Integration,
            //    Oxford, 1994, page 18.
            //
        {
            int dim_num = 2;
            int m = 8;
            int[] z;

            z = new int[dim_num];

            z[0] = 1;
            z[1] = 3;

            Console.WriteLine("");
            Console.WriteLine("TEST12");
            Console.WriteLine("  Lattice.lattice_print prints out the lattice generated");
            Console.WriteLine("  by a single generator vector.");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM = " + dim_num + "");

            typeMethods.i4vec_print(dim_num, z, "  The generator vector:");

            Lattice.lattice_print(dim_num, m, z, "  The total lattice:");


        }

        static void test13()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST13 tests FibonacciLattice.find_z20.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 November 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Ian Sloan, Stephen Joe,
            //    Lattice Methods for Multiple Integration,
            //    Oxford, 1994, page 18.
            //
        {
            int dim;
            int dim_num = 2;
            int i;
            int m;
            int[] z;

            Console.WriteLine("");
            Console.WriteLine("TEST13");
            Console.WriteLine("  FibonacciLattice.find_z20 finds the optimal lattice generator Z");
            Console.WriteLine("  with Fourier coefficient smoothness ALPHA = 2,");
            Console.WriteLine("'  and copy exponent 0,");
            Console.WriteLine("  for a rank 1 \"method of good lattice points\" rule.");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM = " + dim_num + "");
            Console.WriteLine("");
            Console.WriteLine("     M      Z(1)  Z(2)");
            Console.WriteLine("");
            Console.WriteLine("  (M = Fibonacci)");
            Console.WriteLine("");

            for (i = 3; i <= 10; i++)
            {
                m = Helpers.fibonacci(i);

                z = FibonacciLattice.find_z20(dim_num, m);

                string cout = "  " + m.ToString().PadLeft(8);
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += "  " + z[dim].ToString().PadLeft(8);
                }

                Console.WriteLine(cout);

            }

            Console.WriteLine("");
            Console.WriteLine("  (M = 2**K)");
            Console.WriteLine("");

            for (i = 2; i <= 10; i++)
            {
                m = (int) Math.Pow(2, i);

                z = FibonacciLattice.find_z20(dim_num, m);

                string cout = "  " + m.ToString().PadLeft(8);
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += "  " + z[dim].ToString().PadLeft(8);
                }

                Console.WriteLine(cout);

            }

            Console.WriteLine("");
            Console.WriteLine("  (M = 3*2**K)");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                m = 3 * (int) Math.Pow(2, i);

                z = FibonacciLattice.find_z20(dim_num, m);

                string cout = "  " + m.ToString().PadLeft(8);
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += "  " + z[dim].ToString().PadLeft(8);
                }

                Console.WriteLine(cout);

            }

            Console.WriteLine("");
            Console.WriteLine("  (M = Prime)");
            Console.WriteLine("");

            for (i = 3; i <= 10; i++)
            {
                m = Prime.prime(10 * i);

                z = FibonacciLattice.find_z20(dim_num, m);

                string cout = "  " + m.ToString().PadLeft(8);
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += "  " + z[dim].ToString().PadLeft(8);
                }

                Console.WriteLine(cout);

            }

        }

        static void test14()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST14 tests FIBONACCI_LATTICE_Q_NODES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    26 January 2005.
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Ian Sloan, Stephen Joe,
            //    Lattice Methods for Multiple Integration,
            //    Oxford, 1994, page 18.
            //
        {
            int dim_num = 2;
            int k;
            int m;
            double[] x;

            k = 12;
            m = Helpers.fibonacci(k);

            Console.WriteLine("");
            Console.WriteLine("TEST14");
            Console.WriteLine("  FIBONACCI_LATTICE_Q_NODES...");
            Console.WriteLine("");
            Console.WriteLine("  The spatial dimension DIM_NUM = " + dim_num + "");
            Console.WriteLine("  The Fibonacci index K =   " + k + "");
            Console.WriteLine("  The Fibonacci value M =   " + m + "");

            x = FibonacciLattice.fibonacci_lattice_q_nodes(k);

            typeMethods.r8mat_transpose_print(dim_num, m, x, "  The Fibonacci lattice nodes:");

        }
    }
}