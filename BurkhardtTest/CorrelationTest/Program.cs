﻿using System;
using Burkardt;
using Burkardt.CorrelationNS;
using Burkardt.Types;

namespace CorrelationTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for CORRELATION_TEST.
            //
            //  Discussion:
            //
            //    CORRELATION_TEST tests the CORRELATION library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 November 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("CORRELATION_TEST");
            Console.WriteLine("  Test the CORRELATION library.");

            correlation_test01();
            correlation_test02();
            correlation_test03();
            correlation_test04();
            correlation_test05();
            correlation_test06();

            Console.WriteLine("");
            Console.WriteLine("CORRELATION_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }


        static void correlation_test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CORRELATION_TEST01 plots the correlation functions.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 November 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] c;
            int n;
            double[] rho;
            double rho0;

            Console.WriteLine("");
            Console.WriteLine("CORRELATION_TEST01");
            Console.WriteLine("  Make plots of correlation functions.");
            Console.WriteLine("");

            n = 101;
            //
            //  besselj
            //
            rho0 = 1.0;
            rho = typeMethods.r8vec_linspace_new(n, -8.0, 8.0);
            c = Correlation.correlation_besselj(n, rho, rho0);
            Plot.correlation_plot(n, rho, c, "besselj", "Bessel J correlation");
            //
            //  besselk
            //
            rho0 = 1.0;
            rho = typeMethods.r8vec_linspace_new(n, -4.0, 4.0);
            c = Correlation.correlation_besselk(n, rho, rho0);
            Plot.correlation_plot(n, rho, c, "besselk", "Bessel K correlation");
            //
            //  circular
            //
            rho0 = 1.0;
            rho = typeMethods.r8vec_linspace_new(n, -2.0, 2.0);
            c = Correlation.correlation_circular(n, rho, rho0);
            Plot.correlation_plot(n, rho, c, "circular", "Circular correlation");
            //
            //  constant
            //
            rho0 = 1.0;
            rho = typeMethods.r8vec_linspace_new(n, -2.0, 2.0);
            c = Correlation.correlation_constant(n, rho, rho0);
            Plot.correlation_plot(n, rho, c, "constant", "Constant correlation");
            //
            //  cubic
            //
            rho0 = 1.0;
            rho = typeMethods.r8vec_linspace_new(n, -2.0, 2.0);
            c = Correlation.correlation_cubic(n, rho, rho0);
            Plot.correlation_plot(n, rho, c, "cubic", "Cubic correlation");
            //
            //  damped_cosine
            //
            rho0 = 1.0;
            rho = typeMethods.r8vec_linspace_new(n, -6.0, 6.0);
            c = Correlation.correlation_damped_cosine(n, rho, rho0);
            Plot.correlation_plot(n, rho, c, "damped_cosine", "Damped cosine correlation");
            //
            //  damped_sine
            //
            rho0 = 1.0;
            rho = typeMethods.r8vec_linspace_new(n, -12.0, 12.0);
            c = Correlation.correlation_damped_sine(n, rho, rho0);
            Plot.correlation_plot(n, rho, c, "damped_sine", "Damped sine correlation");
            //
            //  exponential
            //
            rho0 = 1.0;
            rho = typeMethods.r8vec_linspace_new(n, -2.0, 2.0);
            c = Correlation.correlation_exponential(n, rho, rho0);
            Plot.correlation_plot(n, rho, c, "exponential", "Exponential correlation");
            //
            //  gaussian
            //
            rho0 = 1.0;
            rho = typeMethods.r8vec_linspace_new(n, -2.0, 2.0);
            c = Correlation.correlation_gaussian(n, rho, rho0);
            Plot.correlation_plot(n, rho, c, "gaussian", "Gaussian correlation");
            //
            //  hole
            //
            rho0 = 1.0;
            rho = typeMethods.r8vec_linspace_new(n, -6.0, 6.0);
            c = Correlation.correlation_hole(n, rho, rho0);
            Plot.correlation_plot(n, rho, c, "hole", "Hole correlation");
            //
            //  linear
            //
            rho0 = 1.0;
            rho = typeMethods.r8vec_linspace_new(n, -2.0, 2.0);
            c = Correlation.correlation_linear(n, rho, rho0);
            Plot.correlation_plot(n, rho, c, "linear", "Linear correlation");
            //
            //  matern, nu = 2.5
            //
            rho0 = 1.0;
            rho = typeMethods.r8vec_linspace_new(n, -2.0, 2.0);
            c = Correlation.correlation_matern(n, rho, rho0);
            Plot.correlation_plot(n, rho, c, "matern", "Matern correlation (NU = 2.5)");
            //
            //  pentaspherical
            //
            rho0 = 1.0;
            rho = typeMethods.r8vec_linspace_new(n, -2.0, 2.0);
            c = Correlation.correlation_pentaspherical(n, rho, rho0);
            Plot.correlation_plot(n, rho, c, "pentaspherical", "Pentaspherical correlation");
            //
            //  power, e = 2.0
            //
            rho0 = 1.0;
            rho = typeMethods.r8vec_linspace_new(n, -2.0, 2.0);
            c = Correlation.correlation_power(n, rho, rho0);
            Plot.correlation_plot(n, rho, c, "power", "Power correlation");
            //
            //  rational_quadratic
            //
            rho0 = 1.0;
            rho = typeMethods.r8vec_linspace_new(n, -4.0, 4.0);
            c = Correlation.correlation_rational_quadratic(n, rho, rho0);
            Plot.correlation_plot(n, rho, c, "rational_quadratic", "Rational quadratic correlation");
            //
            //  spherical
            //
            rho0 = 1.0;
            rho = typeMethods.r8vec_linspace_new(n, -2.0, 2.0);
            c = Correlation.correlation_spherical(n, rho, rho0);
            Plot.correlation_plot(n, rho, c, "spherical", "Spherical correlation");
            //
            //  white_noise
            //
            rho0 = 1.0;
            rho = typeMethods.r8vec_linspace_new(n, -2.0, 2.0);
            c = Correlation.correlation_white_noise(n, rho, rho0);
            Plot.correlation_plot(n, rho, c, "white_noise", "White noise correlation");

            return;
        }


        static void correlation_test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CORRELATION_TEST02 plots sample paths with SAMPLE_PATHS_CHOLESKY.
            //
            //  Discussion:
            //
            //    Most paths will be blue, but make the LAST one red so that there will
            //    always be one distinguished path that is easy to follow.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 November 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n;
            int n2;
            double[] rho;
            double rho0 = 1.0;
            double rhomax = 10.0;
            double rhomin = 0.0;
            int seed;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("CORRELATION_TEST02");
            Console.WriteLine("  SAMPLE_PATHS_CHOLESKY generates sample paths from the");
            Console.WriteLine("  correlation matrix, factored using the Cholesky factor.");
            Console.WriteLine("  It requires that the correlation matrix is nonnegative definite.");
            Console.WriteLine("");
            Console.WriteLine("  SAMPLE_PATHS_EIGEN generates sample paths from the");
            Console.WriteLine("  correlation matrix, factored using the eigen factorization.");
            Console.WriteLine("  If the correlation matrix is not nonnegative definite,");
            Console.WriteLine("  we simply suppress negative eigenvalues.");
            Console.WriteLine("");

            n = 101;
            n2 = 3;
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            //
            //  besselj
            //  Use EIGEN, because CHOLESKY fails.
            //
            seed = 123456789;
            x = SamplePaths.sample_paths_eigen(n, n2, rhomax, rho0, Correlation.correlation_besselj, ref seed);
            Paths.paths_plot(n, n2, rho, x, "besselj", "Bessel J correlation");
            //
            //  besselk
            //
            seed = 123456789;
            x = SamplePaths.sample_paths_cholesky(n, n2, rhomax, rho0, Correlation.correlation_besselk, ref seed);
            Paths.paths_plot(n, n2, rho, x, "besselk", "Bessel K correlation");
            //
            //  circular
            //
            seed = 123456789;
            x = SamplePaths.sample_paths_cholesky(n, n2, rhomax, rho0, Correlation.correlation_circular, ref seed);
            Paths.paths_plot(n, n2, rho, x, "circular", "Circular correlation");
            //
            //  constant
            //
            seed = 123456789;
            x = SamplePaths.sample_paths_cholesky(n, n2, rhomax, rho0, Correlation.correlation_constant, ref seed);
            Paths.paths_plot(n, n2, rho, x, "constant", "Constant correlation");
            //
            //  cubic
            //
            seed = 123456789;
            x = SamplePaths.sample_paths_cholesky(n, n2, rhomax, rho0, Correlation.correlation_cubic, ref seed);
            Paths.paths_plot(n, n2, rho, x, "cubic", "Cubic correlation");
            //
            //  damped_cosine
            //
            seed = 123456789;
            x = SamplePaths.sample_paths_cholesky(n, n2, rhomax, rho0, Correlation.correlation_damped_cosine, ref seed);
            Paths.paths_plot(n, n2, rho, x, "damped_cosine", "Damped cosine correlation");
            //
            //  damped_sine
            //  Use EIGEN, because CHOLESKY fails.
            //
            seed = 123456789;
            x = SamplePaths.sample_paths_eigen(n, n2, rhomax, rho0, Correlation.correlation_damped_sine, ref seed);
            Paths.paths_plot(n, n2, rho, x, "damped_sine", "Damped sine correlation");
            //
            //  exponential
            //
            seed = 123456789;
            x = SamplePaths.sample_paths_cholesky(n, n2, rhomax, rho0, Correlation.correlation_exponential, ref seed);
            Paths.paths_plot(n, n2, rho, x, "exponential", "Exponential correlation");
            //
            //  gaussian
            //  Use EIGEN, because CHOLESKY fails.
            //
            seed = 123456789;
            x = SamplePaths.sample_paths_eigen(n, n2, rhomax, rho0, Correlation.correlation_gaussian, ref seed);
            Paths.paths_plot(n, n2, rho, x, "gaussian", "Gaussian correlation");
            //
            //  hole
            //
            seed = 123456789;
            x = SamplePaths.sample_paths_cholesky(n, n2, rhomax, rho0, Correlation.correlation_hole, ref seed);
            Paths.paths_plot(n, n2, rho, x, "hole", "Hole correlation");
            //
            //  linear
            //
            seed = 123456789;
            x = SamplePaths.sample_paths_cholesky(n, n2, rhomax, rho0, Correlation.correlation_linear, ref seed);
            Paths.paths_plot(n, n2, rho, x, "linear", "Linear correlation");
            //
            //  matern ( nu = 2.5 )
            //
            seed = 123456789;
            x = SamplePaths.sample_paths_cholesky(n, n2, rhomax, rho0, Correlation.correlation_matern, ref seed);
            Paths.paths_plot(n, n2, rho, x, "matern", "Matern correlation (nu=2.5)");
            //
            //  pentaspherical
            //
            seed = 123456789;
            x = SamplePaths.sample_paths_cholesky(n, n2, rhomax, rho0, Correlation.correlation_pentaspherical, ref seed);
            Paths.paths_plot(n, n2, rho, x, "pentaspherical", "Pentaspherical correlation");
            //
            //  power ( e = 2.0 )
            //
            seed = 123456789;
            x = SamplePaths.sample_paths_cholesky(n, n2, rhomax, rho0, Correlation.correlation_power, ref seed);
            Paths.paths_plot(n, n2, rho, x, "power", "Power correlation (e=2.0)");
            //
            //  rational_quadratic
            //  Use EIGEN, because CHOLESKY fails.
            //
            seed = 123456789;
            x = SamplePaths.sample_paths_eigen(n, n2, rhomax, rho0, Correlation.correlation_rational_quadratic, ref seed);
            Paths.paths_plot(n, n2, rho, x, "rational_quadratic", "Rational quadratic correlation");
            //
            //  spherical
            //
            seed = 123456789;
            x = SamplePaths.sample_paths_cholesky(n, n2, rhomax, rho0, Correlation.correlation_spherical, ref seed);
            Paths.paths_plot(n, n2, rho, x, "spherical", "Spherical correlation");
            //
            //  white_noise
            //
            seed = 123456789;
            x = SamplePaths.sample_paths_cholesky(n, n2, rhomax, rho0, Correlation.correlation_white_noise, ref seed);
            Paths.paths_plot(n, n2, rho, x, "white_noise", "White noise correlation");


            return;
        }


        static void correlation_test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CORRELATION_TEST03 plots a correlation function for several values of RH00.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 November 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] c;
            double[] cvec;
            int i;
            int j;
            int n;
            int n2;
            double[] rho;
            double[] rho0;
            double rhomax;
            double rhomin;

            Console.WriteLine("");
            Console.WriteLine("CORRELATION_TEST03");
            Console.WriteLine("  Make plots of correlation functions for");
            Console.WriteLine("  a range of correlation lengths.");
            Console.WriteLine("");
            //
            //  besselj
            //
            n = 101;
            n2 = 5;
            rho0 = new double[n2];
            rho0[0] = 1.0;
            rho0[1] = 1.5;
            rho0[2] = 2.0;
            rho0[3] = 4.0;
            rho0[4] = 8.0;
            rhomin = -8.0;
            rhomax = +8.0;
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            c = new double[n * n2];
            for (j = 0; j < n2; j++)
            {
                cvec = Correlation.correlation_besselj(n, rho, rho0[j]);
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = cvec[i];
                }

            }

            Plot.correlation_plots(n, n2, rho, rho0, c, "besselj", "Bessel J correlation");
            //
            //  besselk
            //
            n = 101;
            n2 = 5;
            rho0 = new double[n2];
            rho0[0] = 1.0;
            rho0[1] = 1.5;
            rho0[2] = 2.0;
            rho0[3] = 4.0;
            rho0[4] = 8.0;
            rhomin = -4.0;
            rhomax = +4.0;
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            c = new double[n * n2];
            for (j = 0; j < n2; j++)
            {
                cvec = Correlation.correlation_besselk(n, rho, rho0[j]);
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = cvec[i];
                }

            }

            Plot.correlation_plots(n, n2, rho, rho0, c, "besselk", "Bessel K correlation");
            //
            //  circular
            //
            n = 101;
            n2 = 6;
            rho0 = new double[n2];
            rho0[0] = 0.5;
            rho0[1] = 1.0;
            rho0[2] = 1.5;
            rho0[3] = 2.0;
            rho0[4] = 4.0;
            rho0[5] = 8.0;
            rhomin = -4.0;
            rhomax = +4.0;
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            c = new double[n * n2];
            for (j = 0; j < n2; j++)
            {
                cvec = Correlation.correlation_circular(n, rho, rho0[j]);
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = cvec[i];
                }
            }

            Plot.correlation_plots(n, n2, rho, rho0, c, "circular", "Circular correlation");
            //
            //  constant
            //  1 plot is enough
            //
            n = 101;
            n2 = 1;
            rho0 = new double[n2];
            rho0[0] = 1.0;
            rhomin = -2.0;
            rhomax = +2.0;
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            c = new double[n * n2];
            for (j = 0; j < n2; j++)
            {
                cvec = Correlation.correlation_constant(n, rho, rho0[j]);
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = cvec[i];
                }

            }

            Plot.correlation_plots(n, n2, rho, rho0, c, "constant", "Constant correlation");
            //
            //  cubic
            //
            n = 101;
            n2 = 6;
            rho0 = new double[n2];
            rho0[0] = 0.5;
            rho0[1] = 1.0;
            rho0[2] = 1.5;
            rho0[3] = 2.0;
            rho0[4] = 4.0;
            rho0[5] = 8.0;
            rhomin = -8.0;
            rhomax = +8.0;
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            c = new double[n * n2];
            for (j = 0; j < n2; j++)
            {
                cvec = Correlation.correlation_cubic(n, rho, rho0[j]);
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = cvec[i];
                }
            }

            Plot.correlation_plots(n, n2, rho, rho0, c, "cubic", "Cubic correlation");
            //
            //  damped_cosine
            //
            n = 101;
            n2 = 6;
            rho0 = new double[n2];
            rho0[0] = 0.5;
            rho0[1] = 1.0;
            rho0[2] = 1.5;
            rho0[3] = 2.0;
            rho0[4] = 4.0;
            rho0[5] = 8.0;
            rhomin = -6.0;
            rhomax = +6.0;
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            c = new double[n * n2];
            for (j = 0; j < n2; j++)
            {
                cvec = Correlation.correlation_damped_cosine(n, rho, rho0[j]);
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = cvec[i];
                }

            }

            Plot.correlation_plots(n, n2, rho, rho0, c, "damped_cosine", "Damped cosine correlation");
            //
            //  damped_sine
            //
            n = 101;
            n2 = 6;
            rho0 = new double[n2];
            rho0[0] = 0.5;
            rho0[1] = 1.0;
            rho0[2] = 1.5;
            rho0[3] = 2.0;
            rho0[4] = 4.0;
            rho0[5] = 8.0;
            rhomin = -8.0;
            rhomax = +8.0;
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            c = new double[n * n2];
            for (j = 0; j < n2; j++)
            {
                cvec = Correlation.correlation_damped_sine(n, rho, rho0[j]);
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = cvec[i];
                }
            }

            Plot.correlation_plots(n, n2, rho, rho0, c, "damped_sine", "Damped sine correlation");
            //
            //  exponential
            //
            n = 101;
            n2 = 7;
            rho0 = new double[n2];
            rho0[0] = 0.25;
            rho0[1] = 0.5;
            rho0[2] = 1.0;
            rho0[3] = 1.5;
            rho0[4] = 2.0;
            rho0[5] = 4.0;
            rho0[6] = 8.0;
            rhomin = -2.0;
            rhomax = +2.0;
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            c = new double[n * n2];
            for (j = 0; j < n2; j++)
            {
                cvec = Correlation.correlation_exponential(n, rho, rho0[j]);
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = cvec[i];
                }

            }

            Plot.correlation_plots(n, n2, rho, rho0, c, "exponential", "Exponential correlation");
            //
            //  gaussian
            //
            n = 101;
            n2 = 7;
            rho0 = new double[n2];
            rho0[0] = 0.25;
            rho0[1] = 0.5;
            rho0[2] = 1.0;
            rho0[3] = 1.5;
            rho0[4] = 2.0;
            rho0[5] = 4.0;
            rho0[6] = 8.0;
            rhomin = -2.0;
            rhomax = +2.0;
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            c = new double[n * n2];
            for (j = 0; j < n2; j++)
            {
                cvec = Correlation.correlation_gaussian(n, rho, rho0[j]);
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = cvec[i];
                }
            }

            Plot.correlation_plots(n, n2, rho, rho0, c, "gaussian", "Gaussian correlation");
            //
            //  hole
            //
            n = 101;
            n2 = 7;
            rho0 = new double[n2];
            rho0[0] = 0.25;
            rho0[1] = 0.5;
            rho0[2] = 1.0;
            rho0[3] = 1.5;
            rho0[4] = 2.0;
            rho0[5] = 4.0;
            rho0[6] = 8.0;
            rhomin = -2.0;
            rhomax = +2.0;
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            c = new double[n * n2];
            for (j = 0; j < n2; j++)
            {
                cvec = Correlation.correlation_hole(n, rho, rho0[j]);
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = cvec[i];
                }
            }

            Plot.correlation_plots(n, n2, rho, rho0, c, "hole", "Hole correlation");
            //
            //  linear
            //
            n = 101;
            n2 = 6;
            rho0 = new double[n2];
            rho0[0] = 0.5;
            rho0[1] = 1.0;
            rho0[2] = 1.5;
            rho0[3] = 2.0;
            rho0[4] = 4.0;
            rho0[5] = 8.0;
            rhomin = -2.0;
            rhomax = +2.0;
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            c = new double[n * n2];
            for (j = 0; j < n2; j++)
            {
                cvec = Correlation.correlation_linear(n, rho, rho0[j]);
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = cvec[i];
                }
            }

            Plot.correlation_plots(n, n2, rho, rho0, c, "linear", "Linear correlation");
            //
            //  matern, nu = 2.5
            //
            n = 101;
            n2 = 6;
            rho0 = new double[n2];
            rho0[0] = 0.5;
            rho0[1] = 1.0;
            rho0[2] = 1.5;
            rho0[3] = 2.0;
            rho0[4] = 4.0;
            rho0[5] = 8.0;
            rhomin = -2.0;
            rhomax = +2.0;
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            c = new double[n * n2];
            for (j = 0; j < n2; j++)
            {
                cvec = Correlation.correlation_matern(n, rho, rho0[j]);
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = cvec[i];
                }
            }

            Plot.correlation_plots(n, n2, rho, rho0, c, "matern", "Matern correlation (NU = 2.5)");
            //
            //  pentaspherical
            //
            n = 101;
            n2 = 6;
            rho0 = new double[n2];
            rho0[0] = 0.5;
            rho0[1] = 1.0;
            rho0[2] = 1.5;
            rho0[3] = 2.0;
            rho0[4] = 4.0;
            rho0[5] = 8.0;
            rhomin = -2.0;
            rhomax = +2.0;
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            c = new double[n * n2];
            for (j = 0; j < n2; j++)
            {
                cvec = Correlation.correlation_pentaspherical(n, rho, rho0[j]);
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = cvec[i];
                }
            }

            Plot.correlation_plots(n, n2, rho, rho0, c, "pentaspherical", "Pentaspherical correlation");
            //
            //  power, e = 2.0
            //
            n = 101;
            n2 = 6;
            rho0 = new double[n2];
            rho0[0] = 0.5;
            rho0[1] = 1.0;
            rho0[2] = 1.5;
            rho0[3] = 2.0;
            rho0[4] = 4.0;
            rho0[5] = 8.0;
            rhomin = -2.0;
            rhomax = +2.0;
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            c = new double[n * n2];
            for (j = 0; j < n2; j++)
            {
                cvec = Correlation.correlation_power(n, rho, rho0[j]);
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = cvec[i];
                }
            }

            Plot.correlation_plots(n, n2, rho, rho0, c, "power", "Power correlation (E = 2.0)");
            //
            //  rational_quadratic
            //
            n = 101;
            n2 = 6;
            rho0 = new double[n2];
            rho0[0] = 0.5;
            rho0[1] = 1.0;
            rho0[2] = 1.5;
            rho0[3] = 2.0;
            rho0[4] = 4.0;
            rho0[5] = 8.0;
            rhomin = -4.0;
            rhomax = +4.0;
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            c = new double[n * n2];
            for (j = 0; j < n2; j++)
            {
                cvec = Correlation.correlation_rational_quadratic(n, rho, rho0[j]);
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = cvec[i];
                }
            }

            Plot.correlation_plots(n, n2, rho, rho0, c, "rational_quadratic", "Rational quadratic correlation");
            //
            //  spherical
            //
            n = 101;
            n2 = 6;
            rho0 = new double[n2];
            rho0[0] = 0.5;
            rho0[1] = 1.0;
            rho0[2] = 1.5;
            rho0[3] = 2.0;
            rho0[4] = 4.0;
            rho0[5] = 8.0;
            rhomin = -8.0;
            rhomax = +8.0;
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            c = new double[n * n2];
            for (j = 0; j < n2; j++)
            {
                cvec = Correlation.correlation_spherical(n, rho, rho0[j]);
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = cvec[i];
                }
            }

            Plot.correlation_plots(n, n2, rho, rho0, c, "spherical", "Spherical correlation");
            //
            //  white_noise
            //  1 plot is enough
            //
            n = 101;
            n2 = 1;
            rho0 = new double[n2];
            rho0[0] = 1.0;
            rhomin = -2.0;
            rhomax = +2.0;
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            c = new double[n * n2];
            for (j = 0; j < n2; j++)
            {
                cvec = Correlation.correlation_white_noise(n, rho, rho0[j]);
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = cvec[i];
                }
            }

            Plot.correlation_plots(n, n2, rho, rho0, c, "white_noise", "White noise correlation");
        }

        static void correlation_test04()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CORRELATION_TEST04 converts between covariance and correlation matrices.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 November 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] c;
            double[] k;
            double[] k2;
            int n;
            double[] sigma;

            Console.WriteLine("");
            Console.WriteLine("CORRELATION_TEST04");
            Console.WriteLine("  Convert between a correlation and a covariance matrix.");

            n = 5;
            k = Matrix.minij(n, n);

            typeMethods.r8mat_print(n, n, k, "  Covariance matrix K:");

            c = new double[n * n];
            sigma = new double[n];

            Correlation.covariance_to_correlation(n, k, ref c, ref sigma);

            typeMethods.r8mat_print(n, n, c, "  Correlation matrix C:");

            typeMethods.r8vec_print(n, sigma, "  Variances:");

            k2 = Correlation.correlation_to_covariance(n, c, sigma);

            typeMethods.r8mat_print(n, n, k2, "  Recovered covariance matrix K2:");
        }

        static void correlation_test05()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CORRELATION_TEST05 calls CORRELATION_BROWNIAN_DISPLAY.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 November 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("CORRELATION_TEST05");
            Console.WriteLine("  CORRELATION_BROWNIAN_DISPLAY displays 4 slices of");
            Console.WriteLine("  the Brownian correlation function.");

            Correlation.correlation_brownian_display();
        }

        static void correlation_test06()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CORRELATION_TEST06 plots sample paths with SAMPLE_PATHS2_CHOLESKY/EIGEN/FFT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 November 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n;
            int n2;
            double[] rho;
            double rho0;
            double rhomax;
            double rhomin;
            int seed;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("CORRELATION_TEST06");
            Console.WriteLine("  For non-stationary correlation functions:");
            Console.WriteLine("");
            Console.WriteLine("  SAMPLE_PATHS2_CHOLESKY generates sample paths from the");
            Console.WriteLine("  correlation matrix, factored using the Cholesky factor.");
            Console.WriteLine("  It requires that the correlation matrix is nonnegative definite.");
            Console.WriteLine("");
            Console.WriteLine("  SAMPLE_PATHS2_EIGEN generates sample paths from the");
            Console.WriteLine("  correlation matrix, factored using the eigen factorization.");
            Console.WriteLine("  If the correlation matrix is not nonnegative definite,");
            Console.WriteLine("  we simply suppress negative eigenvalues.");
            Console.WriteLine("");
            Console.WriteLine("  SAMPLE_PATHS2_FFT generates sample paths from the");
            Console.WriteLine("  correlation matrix, factored using the FFT factorization");
            Console.WriteLine("  of the correlation matrix after embedding in a circulant.");
            Console.WriteLine("");
            /*
            brownian
            */
            n = 101;
            n2 = 3;
            rhomin = 0.0;
            rhomax = 10.0;
            rho0 = 1.0;
            seed = 123456789;
            x = SamplePaths.sample_paths2_cholesky(n, n2, rhomin, rhomax, rho0, Correlation.correlation_brownian, ref seed);
            rho = typeMethods.r8vec_linspace_new(n, rhomin, rhomax);
            Paths.paths_plot(n, n2, rho, x, "brownian", "Brownian correlation");
        }
    }
}