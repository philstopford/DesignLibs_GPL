using System;
using Burkardt;
using Burkardt.PointsNS;
using Burkardt.Types;

namespace LebesgueTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for LEBESGUE_TEST.
            //
            //  Discussion:
            //
            //    LEBESGUE_TEST tests the LEBESGUE library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 March 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("LEBESGUE_TEST");
            Console.WriteLine("  Test the LEBESGUE library.");

            test01();
            test02();
            test03();
            test04();
            test05();
            test06();
            test07();
            test08();
            test09();

            Console.WriteLine("");
            Console.WriteLine("LEBESGUE_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        public static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEBESGUE_TEST01 looks at Chebyshev1 points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string filename = "chebyshev1";
            double[] l;
            string label = "Chebyshev1 points for N = 11";
            int n;
            int n_max = 11;
            int nfun = 501;
            double[] x;
            double[] xfun;

            Console.WriteLine("");
            Console.WriteLine("LEBESGUE_TEST01:");
            Console.WriteLine("  Analyze Chebyshev1 points.");

            xfun = typeMethods.r8vec_linspace_new(nfun, -1.0, +1.0);

            l = new double[nfun];

            for (n = 1; n <= n_max; n++)
            {
                x = Points.chebyshev1(n);
                l[n - 1] = Lebesgue.lebesgue_constant(n, x, nfun, xfun);

            }

            typeMethods.r8vec_print(n_max, l,
                "  Chebyshev1 Lebesgue constants for N = 1 to 11:");
            //
            //  Examine one case more closely.
            //
            n = 11;
            x = Points.chebyshev1(n);
            typeMethods.r8vec_print(n, x, "  Chebyshev1 points for N = 11");

            Lebesgue.lebesgue_plot(n, x, nfun, xfun, label, filename);
        }

        public static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEBESGUE_TEST02 looks at Chebyshev2 points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string filename = "chebyshev2";
            double[] l;
            string label = "Chebyshev2 points for N = 11";
            int n;
            int n_max = 11;
            int nfun = 501;
            double[] x;
            double[] xfun;

            Console.WriteLine("");
            Console.WriteLine("LEBESGUE_TEST02:");
            Console.WriteLine("  Analyze Chebyshev2 points.");

            xfun = typeMethods.r8vec_linspace_new(nfun, -1.0, +1.0);

            l = new double[nfun];

            for (n = 1; n <= n_max; n++)
            {
                x = Points.chebyshev2(n);
                l[n - 1] = Lebesgue.lebesgue_constant(n, x, nfun, xfun);

            }

            typeMethods.r8vec_print(n_max, l,
                "  Chebyshev2 Lebesgue constants for N = 1 to 11:");
            //
            //  Examine one case more closely.
            //
            n = 11;
            x = Points.chebyshev2(n);
            typeMethods.r8vec_print(n, x, "  Chebyshev2 points for N = 11");

            Lebesgue.lebesgue_plot(n, x, nfun, xfun, label, filename);
        }

        public static void test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEBESGUE_TEST03 looks at Chebyshev3 points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 March 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string filename = "chebyshev3";
            double[] l;
            string label = "Chebyshev3 points for N = 11";
            int n;
            int n_max = 11;
            int nfun = 501;
            double[] x;
            double[] xfun;

            Console.WriteLine("");
            Console.WriteLine("LEBESGUE_TEST03:");
            Console.WriteLine("  Analyze Chebyshev3 points.");

            xfun = typeMethods.r8vec_linspace_new(nfun, -1.0, +1.0);

            l = new double[nfun];

            for (n = 1; n <= n_max; n++)
            {
                x = Points.chebyshev3(n);
                l[n - 1] = Lebesgue.lebesgue_constant(n, x, nfun, xfun);

            }

            typeMethods.r8vec_print(n_max, l,
                "  Chebyshev3 Lebesgue constants for N = 1 to 11:");
            //
            //  Examine one case more closely.
            //
            n = 11;
            x = Points.chebyshev3(n);
            typeMethods.r8vec_print(n, x, "  Chebyshev3 points for N = 11");

            Lebesgue.lebesgue_plot(n, x, nfun, xfun, label, filename);
        }

        public static void test04()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEBESGUE_TEST04 looks at Chebyshev4 points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string filename = "chebyshev4";
            double[] l;
            string label = "Chebyshev4 points for N = 11";
            int n;
            int n_max = 11;
            int nfun = 501;
            double[] x;
            double[] xfun;

            Console.WriteLine("");
            Console.WriteLine("LEBESGUE_TEST04:");
            Console.WriteLine("  Analyze Chebyshev4 points.");

            xfun = typeMethods.r8vec_linspace_new(nfun, -1.0, +1.0);

            l = new double[nfun];

            for (n = 1; n <= n_max; n++)
            {
                x = Points.chebyshev4(n);
                l[n - 1] = Lebesgue.lebesgue_constant(n, x, nfun, xfun);
            }

            typeMethods.r8vec_print(n_max, l,
                "  Chebyshev4 Lebesgue constants for N = 1 to 11:");
            //
            //  Examine one case more closely.
            //
            n = 11;
            x = Points.chebyshev4(n);
            typeMethods.r8vec_print(n, x, "  Chebyshev4 points for N = 11");

            Lebesgue.lebesgue_plot(n, x, nfun, xfun, label, filename);
        }

        public static void test05()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEBESGUE_TEST05 looks at Equidistant1 points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string filename = "equidistant1";
            double[] l;
            string label = "Equidistant1 points for N = 11";
            int n;
            int n_max = 11;
            int nfun = 501;
            double[] x;
            double[] xfun;

            Console.WriteLine("");
            Console.WriteLine("LEBESGUE_TEST05:");
            Console.WriteLine("  Analyze Equidistant1 points.");

            xfun = typeMethods.r8vec_linspace_new(nfun, -1.0, +1.0);

            l = new double[nfun];

            for (n = 1; n <= n_max; n++)
            {
                x = Points.equidistant1(n);
                l[n - 1] = Lebesgue.lebesgue_constant(n, x, nfun, xfun);
            }

            typeMethods.r8vec_print(n_max, l,
                "  Equidistant1 Lebesgue constants for N = 1 to 11:");
            //
            //  Examine one case more closely.
            //
            n = 11;
            x = Points.equidistant1(n);
            typeMethods.r8vec_print(n, x, "  Equidistant1 points for N = 11");

            Lebesgue.lebesgue_plot(n, x, nfun, xfun, label, filename);
        }

        public static void test06()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEBESGUE_TEST06 looks at Equidistant2 points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string filename = "equidistant2";
            double[] l;
            string label = "Equidistant2 points for N = 11";
            int n;
            int n_max = 11;
            int nfun = 501;
            double[] x;
            double[] xfun;

            Console.WriteLine("");
            Console.WriteLine("LEBESGUE_TEST06:");
            Console.WriteLine("  Analyze Equidistant2 points.");

            xfun = typeMethods.r8vec_linspace_new(nfun, -1.0, +1.0);

            l = new double[nfun];

            for (n = 1; n <= n_max; n++)
            {
                x = Points.equidistant2(n);
                l[n - 1] = Lebesgue.lebesgue_constant(n, x, nfun, xfun);
            }

            typeMethods.r8vec_print(n_max, l,
                "  Equidistant2 Lebesgue constants for N = 1 to 11:");
            //
            //  Examine one case more closely.
            //
            n = 11;
            x = Points.equidistant2(n);
            typeMethods.r8vec_print(n, x, "  Equidistant2 points for N = 11");

            Lebesgue.lebesgue_plot(n, x, nfun, xfun, label, filename);
        }

        public static void test07()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEBESGUE_TEST07 looks at Equidistant3 points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string filename = "equidistant3";
            double[] l;
            string label = "Equidistant3 points for N = 11";
            int n;
            int n_max = 11;
            int nfun = 501;
            double[] x;
            double[] xfun;

            Console.WriteLine("");
            Console.WriteLine("LEBESGUE_TEST07:");
            Console.WriteLine("  Analyze Equidistant3 points.");

            xfun = typeMethods.r8vec_linspace_new(nfun, -1.0, +1.0);

            l = new double[nfun];

            for (n = 1; n <= n_max; n++)
            {
                x = Points.equidistant3(n);
                l[n - 1] = Lebesgue.lebesgue_constant(n, x, nfun, xfun);
            }

            typeMethods.r8vec_print(n_max, l,
                "  Equidistant3 Lebesgue constants for N = 1 to 11:");
            //
            //  Examine one case more closely.
            //
            n = 11;
            x = Points.equidistant3(n);
            typeMethods.r8vec_print(n, x, "  Equidistant3 points for N = 11");

            Lebesgue.lebesgue_plot(n, x, nfun, xfun, label, filename);
        }

        public static void test08()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEBESGUE_TEST08 looks at Fejer 1 points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string filename = "fejer1";
            double[] l;
            string label = "Fejer1 points for N = 11";
            int n;
            int n_max = 11;
            int nfun = 501;
            double[] x;
            double[] xfun;

            Console.WriteLine("");
            Console.WriteLine("LEBESGUE_TEST08:");
            Console.WriteLine("  Analyze Fejer1 points.");

            xfun = typeMethods.r8vec_linspace_new(nfun, -1.0, +1.0);

            l = new double[nfun];

            for (n = 1; n <= n_max; n++)
            {
                x = Points.fejer1(n);
                l[n - 1] = Lebesgue.lebesgue_constant(n, x, nfun, xfun);
            }

            typeMethods.r8vec_print(n_max, l,
                "  Fejer1 Lebesgue constants for N = 1 to 11:");
            //
            //  Examine one case more closely.
            //
            n = 11;
            x = Points.fejer1(n);
            typeMethods.r8vec_print(n, x, "  Fejer1 points for N = 11");

            Lebesgue.lebesgue_plot(n, x, nfun, xfun, label, filename);
        }

        public static void test09()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEBESGUE_TEST09 looks at Fejer2 points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 March 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string filename = "fejer2";
            double[] l;
            string label = "Fejer2 points for N = 11";
            int n;
            int n_max = 11;
            int nfun = 501;
            double[] x;
            double[] xfun;

            Console.WriteLine("");
            Console.WriteLine("LEBESGUE_TEST09:");
            Console.WriteLine("  Analyze Fejer2 points.");

            xfun = typeMethods.r8vec_linspace_new(nfun, -1.0, +1.0);

            l = new double[nfun];

            for (n = 1; n <= n_max; n++)
            {
                x = Points.fejer2(n);
                l[n - 1] = Lebesgue.lebesgue_constant(n, x, nfun, xfun);
            }

            typeMethods.r8vec_print(n_max, l,
                "  Fejer2 Lebesgue constants for N = 1 to 11:");
            //
            //  Examine one case more closely.
            //
            n = 11;
            x = Points.fejer2(n);
            typeMethods.r8vec_print(n, x, "  Fejer2 points for N = 11");

            Lebesgue.lebesgue_plot(n, x, nfun, xfun, label, filename);

        }
    }
}