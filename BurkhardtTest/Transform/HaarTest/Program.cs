using System;
using Burkardt.Transform;
using Burkardt.Types;
using Burkardt.Uniform;

namespace HaarTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for HAAR_TEST.
            //
            //  Discussion:
            //
            //    HAAR_TEST tests the HAAR library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 March 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("HAAR_TEST");
            Console.WriteLine("  Test the HAAR library.");

            test01();
            test02();

            Console.WriteLine("");
            Console.WriteLine("HAAR_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests HAAR_1D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 March 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double err;
            int i;
            int n;
            int seed;
            double[] u;
            double[] v;
            double[] w;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  HAAR_1D computes the Haar transform of a vector.");
            //
            //  Random data.
            //
            n = 16;
            seed = 123456789;
            u = UniformRNG.r8vec_uniform_01_new(n, ref seed);
            v = typeMethods.r8vec_copy_new(n, u);

            Haar.haar_1d(n, ref v);

            w = typeMethods.r8vec_copy_new(n, v);
            Haar.haar_1d_inverse(n, ref w);

            Console.WriteLine("");
            Console.WriteLine("   i      U(i)        H(U)(i)  Hinv(H(U))(i)");
            Console.WriteLine("");
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                    + "  " + u[i].ToString().PadLeft(10)
                    + "  " + v[i].ToString().PadLeft(10)
                    + "  " + w[i].ToString().PadLeft(10) + "");
            }

            //
            //  Constant signal.
            //
            n = 8;
            u = typeMethods.r8vec_ones_new(n);
            v = typeMethods.r8vec_copy_new(n, u);

            Haar.haar_1d(n, ref v);

            w = typeMethods.r8vec_copy_new(n, v);
            Haar.haar_1d_inverse(n, ref w);

            Console.WriteLine("");
            Console.WriteLine("   i      U(i)        H(U)(i)  Hinv(H(U))(i)");
            Console.WriteLine("");
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                    + "  " + u[i].ToString().PadLeft(10)
                    + "  " + v[i].ToString().PadLeft(10)
                    + "  " + w[i].ToString().PadLeft(10) + "");
            }

            //
            //  Linear signal.
            //
            n = 16;
            u = typeMethods.r8vec_linspace_new(n, 1.0, (double) n);
            v = typeMethods.r8vec_copy_new(n, u);

            Haar.haar_1d(n, ref v);

            w = typeMethods.r8vec_copy_new(n, v);
            Haar.haar_1d_inverse(n, ref w);

            Console.WriteLine("");
            Console.WriteLine("   i      U(i)        H(U)(i)  Hinv(H(U))(i)");
            Console.WriteLine("");
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                    + "  " + u[i].ToString().PadLeft(10)
                    + "  " + v[i].ToString().PadLeft(10)
                    + "  " + w[i].ToString().PadLeft(10) + "");
            }

            //
            //  Quadratic data.
            //
            n = 8;
            u = new double[n];
            u[0] = 25.0;
            u[1] = 16.0;
            u[2] = 9.0;
            u[3] = 4.0;
            u[4] = 1.0;
            u[5] = 0.0;
            u[6] = 1.0;
            u[7] = 4.0;
            v = typeMethods.r8vec_copy_new(n, u);

            Haar.haar_1d(n, ref v);

            w = typeMethods.r8vec_copy_new(n, v);
            Haar.haar_1d_inverse(n, ref w);

            Console.WriteLine("");
            Console.WriteLine("   i      U(i)        H(U)(i)  Hinv(H(U))(i)");
            Console.WriteLine("");
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                    + "  " + u[i].ToString().PadLeft(10)
                    + "  " + v[i].ToString().PadLeft(10)
                    + "  " + w[i].ToString().PadLeft(10) + "");
            }

            //
            //  N not a power of 2.
            //
            n = 99;
            seed = 123456789;
            u = UniformRNG.r8vec_uniform_01_new(n, ref seed);

            v = typeMethods.r8vec_copy_new(n, u);
            Haar.haar_1d(n, ref v);

            w = typeMethods.r8vec_copy_new(n, v);
            Haar.haar_1d_inverse(n, ref w);

            err = typeMethods.r8vec_diff_norm(n, u, w);

            Console.WriteLine("");
            Console.WriteLine("  For N = " + n
                + ", ||u-Haar.haar_1d_inverse(Haar.haar_1d(u))|| = " + err + "");
        }

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests HAAR_2D and HAAR_2D_INVERSE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 March 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double err;
            int m = 16;
            int n = 4;
            int seed;
            double[] u;
            double[] v;
            double[] w;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  HAAR_2D computes the Haar transform of an array.");
            Console.WriteLine("  HAAR_2D_INVERSE inverts the transform.");
            //
            //  Demonstrate successful inversion.
            //
            seed = 123456789;
            u = UniformRNG.r8mat_uniform_01_new(m, n, ref seed);

            typeMethods.r8mat_print(m, n, u, "  Input array U:");

            v = typeMethods.r8mat_copy_new(m, n, u);

            Haar.haar_2d(m, n, ref v);

            typeMethods.r8mat_print(m, n, v, "  Transformed array V:");

            w = typeMethods.r8mat_copy_new(m, n, v);

            Haar.haar_2d_inverse(m, n, ref w);

            typeMethods.r8mat_print(m, n, w, "  Recovered array W:");

            //
            //  M, N not powers of 2.
            //
            m = 37;
            n = 53;
            seed = 123456789;
            u = UniformRNG.r8mat_uniform_01_new(m, n, ref seed);

            v = typeMethods.r8mat_copy_new(m, n, u);
            Haar.haar_2d(m, n, ref v);

            w = typeMethods.r8mat_copy_new(m, n, v);
            Haar.haar_2d_inverse(m, n, ref w);

            err = typeMethods.r8mat_dif_fro(m, n, u, w);

            Console.WriteLine("");
            Console.WriteLine("  M = " + m
                + ", N = " + n
                + ", ||Haar.haar_2d_inverse(Haar.haar_2d(u))-u|| = " + err + "");
        }
    }
}