using System;
using System.Globalization;
using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Uniform;

namespace BLAS1ZTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for BLAS1_Z_TEST.
        //
        //  Discussion:
        //
        //    BLAS1_Z_TEST tests the BLAS1_Z library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("BLAS1_Z_TEST:");
        Console.WriteLine("  Test the BLAS1_Z library,");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();
        test07();
        test08();
        test09();

        test10();
        test11();
        test12();
        test13();
        test14();
        test15();
        test16();
        test17();

        Console.WriteLine("");
        Console.WriteLine("BLAS1_Z_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests DZASUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int MA = 5;
        int NA = 4;
        int NX = 8;

        Complex[] a =
            {
                new(-3.0, 4.0),
                new(2.0, 0.0),
                new(3.0, -4.0),
                new(2.0, 0.0),
                new(2.0, -1.0),
                new(-1.0, 1.0),
                new(0.0, 5.0),
                new(-4.0, -2.0),
                new(-4.0, 1.0),
                new(-4.0, -3.0),
                new(0.0, -2.0),
                new(1.0, 3.0),
                new(-3.0, 3.0),
                new(-3.0, 3.0),
                new(-1.0, -2.0),
                new(-1.0, 2.0),
                new(2.0, -4.0),
                new(0.0, -1.0),
                new(0.0, -1.0),
                new(-2.0, 4.0)
            }
            ;

        Complex[] x =
            {
                new(2.0, -1.0),
                new(-4.0, -2.0),
                new(3.0, 1.0),
                new(2.0, 2.0),
                new(-1.0, -1.0),
                new(-1.0, 0.0),
                new(0.0, -3.0),
                new(4.0, 0.0)
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  DZASUM adds the absolute values of elements");
        Console.WriteLine("  of a complex vector.");
        Console.WriteLine("");
        Console.WriteLine("  X =");
        Console.WriteLine("");
        for (int i = 0; i < NX; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  DZASUM ( NX,   X, 1    ) = "
                          + BLAS1Z.dzasum(NX, x, 0, 1) + "");
        Console.WriteLine("  DZASUM ( NX/2, X, 2    ) = "
                          + BLAS1Z.dzasum(NX / 2, x, 0, 2) + "");
        Console.WriteLine("  DZASUM ( 2,    X, NX/2 ) = "
                          + BLAS1Z.dzasum(2, x, 0, NX / 2) + "");

        Console.WriteLine("");
        Console.WriteLine("  Demonstrate with a matrix A:");
        Console.WriteLine("");
        for (int i = 0; i < MA; i++)
        {
            string cout = "";
            for (int j = 0; j < NA; j++)
            {
                cout += "  " + a[i + j * MA].ToString(CultureInfo.InvariantCulture).PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  DZASUM ( MA, A[1,2], 1 )   = "
                          + BLAS1Z.dzasum(MA, a, 0 + 1 * MA, 1) + "");
        Console.WriteLine("  DZASUM ( NA, A[2,1], MA ) = "
                          + BLAS1Z.dzasum(NA, a, +1 , MA) + "");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests DZNRM2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;

        Complex[] x =
            {
                new(2.0, -1.0),
                new(-4.0, -2.0),
                new(3.0, 1.0),
                new(2.0, 2.0),
                new(-1.0, -1.0)
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  DZNRM2 returns the Euclidean norm of a complex vector.");

        Console.WriteLine("");
        Console.WriteLine("  The vector X:");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(16) + "");
        }

        int incx = 1;
        double norm = BLAS1Z.dznrm2(N, x, incx);

        Console.WriteLine("");
        Console.WriteLine("  The L2 norm of X is " + norm + "");
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests IZAMAX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;

        Complex[] x =
            {
                new(2.0, -1.0),
                new(-4.0, -2.0),
                new(3.0, 1.0),
                new(2.0, 2.0),
                new(-1.0, -1.0)
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  IZAMAX returns the index of the entry of");
        Console.WriteLine("  maximum magnitude in a complex vector.");

        Console.WriteLine("");
        Console.WriteLine("  The entries and ZABS1 magnitudes:");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(16)
                                   + "  " + BLAS0.zabs1(x[i]).ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }

        int incx = 1;

        int i2 = BLAS1Z.izamax(N, x, incx);

        Console.WriteLine("");
        Console.WriteLine("  The index of maximum magnitude = " + i2 + "");
        Console.WriteLine("");
        Console.WriteLine("  Note that this is a 1-based index.");
        Console.WriteLine("  Note that the L1 norm is used.");

    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests ZABS1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  ZABS1 returns the L1 norm of a complex number.");
        Console.WriteLine("");
        Console.WriteLine("      Real      Imaginary");
        Console.WriteLine("      Part      Part           ZABS1(Z)");
        Console.WriteLine("");

        for (int i = 1; i <= 10; i++)
        {
            //
            //  The compiler seems to be unhappy with the statement
            //
            //    c = 5.0 * c8_uniform_01 ( seed );
            //
            //  Poor compiler.
            //
            Complex c = UniformRNG.c8_uniform_01(ref seed);
            c = new Complex(5.0, 0) * c;

            double c_norm = BLAS0.zabs1(c);

            Console.WriteLine("  " + c.Real.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + c.Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + c_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests ZABS2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  ZABS2 returns the L2 norm of a complex number.");
        Console.WriteLine("");
        Console.WriteLine("      Real      Imaginary");
        Console.WriteLine("      Part      Part           ZABS2(Z");
        Console.WriteLine("");

        for (int i = 1; i <= 10; i++)
        {
            Complex c = UniformRNG.c8_uniform_01(ref seed);
            c = new Complex(5.0, 0) * c;

            double c_norm = BLAS0.zabs2(c);

            Console.WriteLine("  " + c.Real.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + c.Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + c_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }
    //****************************************************************************80

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests ZAXPY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 5;

        Complex[] x =
            {
                new(2.0, -1.0),
                new(-4.0, -2.0),
                new(3.0, 1.0),
                new(2.0, 2.0),
                new(-1.0, -1.0)
            }
            ;
        Complex[] y =
            {
                new(-1.0, 0.0),
                new(0.0, -3.0),
                new(4.0, 0.0),
                new(-3.0, 4.0),
                new(-2.0, 0.0)
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  ZAXPY adds a multiple of one complex vector to another.");

        Console.WriteLine("");
        Console.WriteLine("  X =");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Y =");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + y[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + y[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }

        Complex s = new(0.50, -1.00);

        Console.WriteLine("");
        Console.WriteLine("  The scalar multiplier is: " + s + "");

        BLAS1Z.zaxpy(N, s, x, 1, ref y, 1);

        Console.WriteLine("");
        Console.WriteLine("  A * X + Y =");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + y[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + y[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests ZCOPY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N1 = 5;
        int N2 = 5;
        int N = 10;

        Complex[] a = new Complex[N1 * N2];
        Complex[] x = new Complex[N];
        Complex[] y = new Complex[N];

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  ZCOPY copies one complex vector into another.");

        for (int i = 0; i < N; i++)
        {
            x[i] = new Complex(10 * (i + 1), i + 1);
        }

        for (int i = 0; i < N; i++)
        {
            y[i] = new Complex(20 * (i + 1), 2 * (i + 1));
        }

        for (int i = 0; i < N1; i++)
        {
            for (int j = 0; j < N2; j++)
            {
                a[i + j * N1] = new Complex(10 * (i + 1), j + 1);
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  X =");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Y =");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + y[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + y[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  A =");
        Console.WriteLine("");
        for (int i = 0; i < N1; i++)
        {
            for (int j = 0; j < N2; j++)
            {
                Console.WriteLine("  " + a[i + j * N1].Real.ToString(CultureInfo.InvariantCulture).PadLeft(5)
                                       + "  " + a[i + j * N1].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(5) + "");
            }

            Console.WriteLine("");
        }

        BLAS1Z.zcopy(5, x, 1, ref y, 1);
        Console.WriteLine("");
        Console.WriteLine("  ZCOPY ( 5, X, 1, Y, 1 )");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + y[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + y[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }

        for (int i = 0; i < N; i++)
        {
            y[i] = new Complex(20 * (i + 1), 2 * (i + 1));
        }

        BLAS1Z.zcopy(3, x, 2, ref y, 3);

        Console.WriteLine("");
        Console.WriteLine("  ZCOPY ( 3, X, 2, Y, 3 )");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + y[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + y[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }

        BLAS1Z.zcopy(5, x, 1, ref a, 1);

        Console.WriteLine("");
        Console.WriteLine("  ZCOPY ( 5, X, 1, A, 1 )");
        Console.WriteLine("");
        Console.WriteLine("");
        Console.WriteLine("  A =");
        Console.WriteLine("");
        for (int i = 0; i < N1; i++)
        {
            for (int j = 0; j < N2; j++)
            {
                Console.WriteLine("  " + a[i + j * N1].Real.ToString(CultureInfo.InvariantCulture).PadLeft(5)
                                       + "  " + a[i + j * N1].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(5) + "");
            }

            Console.WriteLine("");
        }

        for (int i = 0; i < N1; i++)
        {
            for (int j = 0; j < N2; j++)
            {
                a[i + j * N1] = new Complex(10 * (i + 1), j + 1);
            }
        }

        BLAS1Z.zcopy(5, x, 2, ref a, 5);

        Console.WriteLine("");
        Console.WriteLine("  ZCOPY ( 5, X, 2, A, 5 )");
        Console.WriteLine("");
        Console.WriteLine("  A =");
        Console.WriteLine("");
        for (int i = 0; i < N1; i++)
        {
            for (int j = 0; j < N2; j++)
            {
                Console.WriteLine("  " + a[i + j * N1].Real.ToString(CultureInfo.InvariantCulture).PadLeft(5)
                                       + "  " + a[i + j * N1].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(5) + "");
            }

            Console.WriteLine("");
        }
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests ZDOTC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 5;

        Complex[] x =
            {
                new(2.0, -1.0),
                new(-4.0, -2.0),
                new(3.0, 1.0),
                new(2.0, 2.0),
                new(-1.0, -1.0)
            }
            ;
        Complex[] y =
            {
                new(-1.0, 0.0),
                new(0.0, -3.0),
                new(4.0, 0.0),
                new(-3.0, 4.0),
                new(-2.0, 0.0)
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  ZDOTC computes the conjugated dot product of");
        Console.WriteLine("  two complex vectors.");

        Console.WriteLine("");
        Console.WriteLine("  X =");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }

        Complex x_norm = BLAS1Z.zdotc(N, x, 1, x, 1);

        Console.WriteLine("");
        Console.WriteLine("  The square of the norm of X, computed as");
        Console.WriteLine("  ZDOTC(X,X) = " + x_norm + "");

        Complex xy_dot = BLAS1Z.zdotc(N, x, 1, y, 1);

        Console.WriteLine("");
        Console.WriteLine("  Y =");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + y[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + y[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  The dot product X.Y* is " + xy_dot + "");
    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests ZDOTU.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 5;

        Complex[] x =
            {
                new(2.0, -1.0),
                new(-4.0, -2.0),
                new(3.0, 1.0),
                new(2.0, 2.0),
                new(-1.0, -1.0)
            }
            ;
        Complex[] y =
            {
                new(-1.0, 0.0),
                new(0.0, -3.0),
                new(4.0, 0.0),
                new(-3.0, 4.0),
                new(-2.0, 0.0)
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  ZDOTU computes the unconjugated dot product of");
        Console.WriteLine("  two complex vectors.");

        Console.WriteLine("");
        Console.WriteLine("  X =");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }

        Complex x_norm = BLAS1Z.zdotu(N, x, 1, x, 1);

        Console.WriteLine("");
        Console.WriteLine("  The unconjugated dot product ( X dot X )");
        Console.WriteLine("  (which is NOT the square of the norm of X!):");
        Console.WriteLine("  ZDOTU(X,X) = " + x_norm + "");

        Complex xy_dot = BLAS1Z.zdotu(N, x, 1, y, 1);

        Console.WriteLine("");
        Console.WriteLine("  Y =");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + y[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + y[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  The dot product ( X dot Y ) is " + xy_dot + "");
    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests ZDROT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 6;

        Complex[] x = new Complex[N];
        Complex[] y = new Complex[N];

        for (int i = 0; i < N; i++)
        {
            x[i] = new Complex(10 * (i + 1), i + 1);
        }

        for (int i = 0; i < N; i++)
        {
            y[i] = new Complex(20 * (i + 1), 2 * (i + 1));
        }

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  ZDROT carries out a Givens rotation");
        Console.WriteLine("  on a complex vector.");
        Console.WriteLine("");
        Console.WriteLine("  X and Y");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(20) + "");
        }

        double c = 0.5;
        double s = Math.Sqrt(1.0 - c * c);
        BLAS1Z.zdrot(N, ref x, 1, ref y, 1, c, s);
        Console.WriteLine("");
        Console.WriteLine("  ZDROT ( N, X, 1, Y, 1, " + c + "," + s + " )");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(20) + "");
        }
    }

    private static void test11()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 tests ZDSCAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 6;

        Complex[] x = new Complex[N];

        for (int i = 0; i < N; i++)
        {
            x[i] = new Complex(10 * (i + 1), i + 1);
        }

        Console.WriteLine("");
        Console.WriteLine("TEST11");
        Console.WriteLine("  ZDSCAL multiplies a real scalar times a complex vector.");

        Console.WriteLine("");
        Console.WriteLine("  X =");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }

        double da = 5.0;
        BLAS1Z.zdscal(N, da, ref x, 1);
        Console.WriteLine("");
        Console.WriteLine("  ZDSCAL ( N, " + da + ", X, 1 )");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }

        for (int i = 0; i < N; i++)
        {
            x[i] = new Complex(10 * (i + 1), i + 1);
        }

        da = -2.0;
        BLAS1Z.zdscal(3, da, ref x, 2);
        Console.WriteLine("");
        Console.WriteLine("  ZDSCAL ( 3, " + da + ", X, 2 )");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }
    }

    private static void test12()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST12 tests ZMACH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TEST12");
        Console.WriteLine("  ZMACH computes several machine-dependent");
        Console.WriteLine("  complex arithmetic parameters.");

        Console.WriteLine("");
        Console.WriteLine("  ZMACH(1)  = machine epsilon = " + BLAS0.zmach(1) + "");
        Console.WriteLine("  ZMACH(2)  = a tiny value    = " + BLAS0.zmach(2) + "");
        Console.WriteLine("  ZMACH(3)  = a huge value    = " + BLAS0.zmach(3) + "");
    }

    private static void test13()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST13 tests ZROTG.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double c = 0;
        Complex s = new();
        int test_num = 5;

        Console.WriteLine("");
        Console.WriteLine("TEST13");
        Console.WriteLine("  ZROTG generates a complex Givens rotation");
        Console.WriteLine("    (  C  S ) * ( A ) = ( R )");
        Console.WriteLine("    ( -S  C )   ( B )   ( 0 )");
        Console.WriteLine("");

        int seed = 123456789;

        for (int test = 1; test <= test_num; test++)
        {
            Complex a = UniformRNG.c8_uniform_01(ref seed);
            Complex b = UniformRNG.c8_uniform_01(ref seed);

            Complex sa = a;
            Complex sb = b;

            BLAS1Z.zrotg(ref sa, sb, ref c, ref s);

            Complex r = sa;

            Console.WriteLine("");
            Console.WriteLine("  A =  " + a + "");
            Console.WriteLine("  B =  " + b + "");
            Console.WriteLine("  C =  " + c + "");
            Console.WriteLine("  S =  " + s + "");
            Console.WriteLine("  R =  " + r + "");
            Console.WriteLine("         C *A+S*B = " + (c * a + s * b) + "");
            Console.WriteLine("  -conjg(S)*A+C*B = " + (-Complex.Conjugate(s) * a + c * b) + "");
        }
    }

    private static void test14()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST14 tests ZSCAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 6;

        Complex[] x = new Complex[N];

        for (int i = 0; i < N; i++)
        {
            x[i] = new Complex(10 * (i + 1), i + 1);
        }

        Console.WriteLine("");
        Console.WriteLine("TEST14");
        Console.WriteLine("  ZSCAL multiplies a complex scalar times a vector.");

        Console.WriteLine("");
        Console.WriteLine("  X =");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }

        Complex da = new(5.0, 0.0);
        BLAS1Z.zscal(N, da, ref x, 1);
        Console.WriteLine("");
        Console.WriteLine("  ZSCAL ( N, (" + da + "), X, 1 )");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }

        for (int i = 0; i < N; i++)
        {
            x[i] = new Complex(10 * (i + 1), i + 1);
        }

        da = new Complex(-2.0, 1.0);
        BLAS1Z.zscal(3, da, ref x, 2);
        Console.WriteLine("");
        Console.WriteLine("  ZSCAL ( 3, (" + da + "), X, 2 )");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Real.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].Imaginary.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }
    }

    private static void test15()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST15 tests ZSIGN1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST15");
        Console.WriteLine("  ZSIGN1 ( C1, C2 ) transfers the sign of complex C2");
        Console.WriteLine("  to the ZABS1 magnitude of C1.");
        Console.WriteLine("");
        Console.WriteLine(
            "           C1                    C2                    C3");
        Console.WriteLine(
            "  --------------------  --------------------  --------------------");
        Console.WriteLine();

        for (int i = 1; i <= 10; i++)
        {
            Complex c1 = new Complex(5.0, 0) * UniformRNG.c8_uniform_01(ref seed);
            Complex c2 = new Complex(5.0, 0) * UniformRNG.c8_uniform_01(ref seed);
            Complex c3 = BLAS0.zsign1(c1, c2);

            Console.WriteLine("  " + c1.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + c2.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + c3.ToString(CultureInfo.InvariantCulture).PadLeft(20) + "");
        }
    }

    private static void test16()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST16 tests ZSIGN2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST16");
        Console.WriteLine("  ZSIGN2 ( C1, C2 ) transfers the sign of complex C2");
        Console.WriteLine("  to the ZABS2 magnitude of C1.");
        Console.WriteLine("");
        Console.WriteLine(
            "           C1                    C2                    C3");
        Console.WriteLine(
            "  --------------------  --------------------  --------------------");
        Console.WriteLine("");

        for (int i = 1; i <= 10; i++)
        {
            Complex c1 = new Complex(5.0, 0) * UniformRNG.c8_uniform_01(ref seed);
            Complex c2 = new Complex(5.0, 0) * UniformRNG.c8_uniform_01(ref seed);
            Complex c3 = BLAS0.zsign2(c1, c2);

            Console.WriteLine("  " + c1.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + c2.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + c3.ToString(CultureInfo.InvariantCulture).PadLeft(20) + "");
        }
    }

    private static void test17()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST17 tests ZSWAP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 5;

        Complex[] x = new Complex[N];
        Complex[] y = new Complex[N];

        for (int i = 0; i < N; i++)
        {
            x[i] = new Complex(10 * (i + 1), i + 1);
        }

        for (int i = 0; i < N; i++)
        {
            y[i] = new Complex(20 * (i + 1), 2 * (i + 1));
        }

        Console.WriteLine("");
        Console.WriteLine("TEST17");
        Console.WriteLine("  ZSWAP swaps two complex vectors.");
        Console.WriteLine("");
        Console.WriteLine("  X and Y");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(20) + "");
        }

        BLAS1Z.zswap(N, ref x, 1, ref y, 1);
        Console.WriteLine("");
        Console.WriteLine("  ZSWAP ( N, X, 1, Y, 1 )");
        Console.WriteLine("");
        Console.WriteLine("  X and Y");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(20) + "");
        }

        for (int i = 0; i < N; i++)
        {
            x[i] = new Complex(10 * (i + 1), i + 1);
        }

        for (int i = 0; i < N; i++)
        {
            y[i] = new Complex(20 * (i + 1), 2 * (i + 1));
        }

        BLAS1Z.zswap(3, ref x, 2, ref y, 1);
        Console.WriteLine("");
        Console.WriteLine("  ZSWAP ( 3, X, 2, Y, 1 )");
        Console.WriteLine("");
        Console.WriteLine("  X and Y");
        Console.WriteLine("");
        for (int i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(20) + "");
        }
    }
}