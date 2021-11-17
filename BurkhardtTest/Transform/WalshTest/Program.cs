using System;
using Burkardt.Transform;
using Burkardt.Types;
using Burkardt.Uniform;

namespace WalshTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for WALSH_TEST.
        //
        //  Discussion:
        //
        //    WALSH_TEST tests the WALSH library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("WALSH_TEST");
        Console.WriteLine("  Test the WALSH library.");

        test01();
        test02();
        test03();
        test04();

        Console.WriteLine("");
        Console.WriteLine("WALSH_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests FWT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int n = 16;
        int seed;
        double[] w;
        double[] x;
        double[] y;
        double[] z;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  FWT computes a fast Walsh transform.");

        for (j = 1; j <= 2; j++)
        {
            switch (j)
            {
                case 1:
                    seed = 123456789;
                    w = UniformRNG.r8vec_uniform_01_new(n, ref seed);
                    break;
                default:
                {
                    w = new double[n];
                    for (i = 0; i < n; i++)
                    {
                        w[i] = i + 1;
                    }

                    break;
                }
            }

            x = typeMethods.r8vec_copy_new(n, w);
            Walsh.fwt(n, ref w);
            y = typeMethods.r8vec_copy_new(n, w);
            for (i = 0; i < n; i++)
            {
                y[i] /= n;
            }

            Walsh.fwt(n, ref w);
            z = typeMethods.r8vec_copy_new(n, w);
            for (i = 0; i < n; i++)
            {
                z[i] /= n;
            }

            Console.WriteLine("");
            Console.WriteLine("     I        X(I)    Y=FWT(X)/N   Z=FWT(Y)/N");
            Console.WriteLine("");
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                                       + "  " + x[i].ToString().PadLeft(10)
                                       + "  " + y[i].ToString().PadLeft(10)
                                       + "  " + z[i].ToString().PadLeft(10) + "");
            }
        }

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests WALSH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int n = 16;
        int seed;
        double[] w;
        double[] x;
        double[] y;
        double[] z;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  WALSH computes a fast Walsh transform.");

        for (j = 1; j <= 2; j++)
        {
            switch (j)
            {
                case 1:
                    seed = 123456789;
                    w = UniformRNG.r8vec_uniform_01_new(n, ref seed);
                    break;
                default:
                {
                    w = new double[n];
                    for (i = 0; i < n; i++)
                    {
                        w[i] = i + 1;
                    }

                    break;
                }
            }

            x = typeMethods.r8vec_copy_new(n, w);
            Walsh.walsh(n, ref w);
            y = typeMethods.r8vec_copy_new(n, w);
            for (i = 0; i < n; i++)
            {
                y[i] /= n;
            }

            Walsh.walsh(n, ref w);
            z = typeMethods.r8vec_copy_new(n, w);
            for (i = 0; i < n; i++)
            {
                z[i] /= n;
            }

            Console.WriteLine("");
            Console.WriteLine("     I        X(I)    Y=FWT(X)/N   Z=FWT(Y)/N");
            Console.WriteLine("");
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                                       + "  " + x[i].ToString().PadLeft(10)
                                       + "  " + y[i].ToString().PadLeft(10)
                                       + "  " + z[i].ToString().PadLeft(10) + "");
            }
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests HAAR, HAARIN and HNORM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int n = 16;
        int seed;
        double[] w;
        double[] x;
        double[] y;
        double[] z;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  HAAR computes a Haar transform.");
        Console.WriteLine("  HNORM normalizes the transformed data.");
        Console.WriteLine("  HAARIN computes an inverse Haar transform.");

        for (j = 1; j <= 2; j++)
        {
            switch (j)
            {
                case 1:
                    seed = 123456789;
                    w = UniformRNG.r8vec_uniform_01_new(n, ref seed);
                    break;
                default:
                {
                    w = new double[n];
                    for (i = 0; i < n; i++)
                    {
                        w[i] = i + 1;
                    }

                    break;
                }
            }

            x = typeMethods.r8vec_copy_new(n, w);

            Haar.haar(n, ref w);

            y = typeMethods.r8vec_copy_new(n, w);

            Haar.hnorm(n, ref w);

            z = typeMethods.r8vec_copy_new(n, w);

            Haar.haarin(n, ref w);

            Console.WriteLine("");
            Console.WriteLine("     I        X(I)    Y=HAAR(X)  Z=HNORM(Y)  W=HAARIN(Z)");
            Console.WriteLine("");
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                                       + "  " + x[i].ToString().PadLeft(10)
                                       + "  " + y[i].ToString().PadLeft(10)
                                       + "  " + z[i].ToString().PadLeft(10)
                                       + "  " + w[i].ToString().PadLeft(10) + "");
            }
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests FFWT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int n = 16;
        int seed;
        double[] w;
        double[] x;
        double[] y;
        double[] z;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  FFWT computes a fast Walsh transform.");

        for (j = 1; j <= 2; j++)
        {
            switch (j)
            {
                case 1:
                    seed = 123456789;
                    w = UniformRNG.r8vec_uniform_01_new(n, ref seed);
                    break;
                default:
                {
                    w = new double[n];
                    for (i = 0; i < n; i++)
                    {
                        w[i] = i + 1;
                    }

                    break;
                }
            }

            x = typeMethods.r8vec_copy_new(n, w);
            Walsh.ffwt(n, ref w);
            y = typeMethods.r8vec_copy_new(n, w);
            for (i = 0; i < n; i++)
            {
                y[i] /= n;
            }

            Walsh.ffwt(n, ref w);
            z = typeMethods.r8vec_copy_new(n, w);
            for (i = 0; i < n; i++)
            {
                z[i] /= n;
            }

            Console.WriteLine("");
            Console.WriteLine("     I        X(I)   Y=FFWT(X)/N  Z=FFWT(Y)/N");
            Console.WriteLine("");
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(4)
                                       + "  " + x[i].ToString().PadLeft(10)
                                       + "  " + y[i].ToString().PadLeft(10)
                                       + "  " + z[i].ToString().PadLeft(10) + "");
            }
        }
    }
}