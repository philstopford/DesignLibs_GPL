using System;
using Burkardt.IHSNS;
using Burkardt.Types;

namespace IHSTest;

internal static class Program
{
    private static void Main()
    {
        Console.WriteLine();
            
        Console.WriteLine("IHS_TEST");
        Console.WriteLine("  C# version");
        Console.WriteLine("  Test the IHS library.");

        test01();
        test02();
        test03();
        test04();

        Console.WriteLine();
        Console.WriteLine("IHS_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine();
            
    }

    private static void test01( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests the improved distributed hypercube sampling algorithm.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 April 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim_num;
        const int duplication = 5;
        const int point_num = 10;

        Console.WriteLine();
        Console.WriteLine("TEST01");
        Console.WriteLine("  IHS implements the IHS Algorithm");
        Console.WriteLine("  (Improved Distributed Hypercube Sampling)");
        Console.WriteLine();
        Console.WriteLine("  Demonstrate the code for a fixed number of points");
        Console.WriteLine("  and an increasing dimension.");

        for ( dim_num = 1; dim_num <= 4; dim_num++ )
        {
            int seed = 17;

            double opt = point_num /
                         Math.Pow ( point_num,
                             1.0E+00 / dim_num );

            Console.WriteLine();
            Console.WriteLine("  Random number seed =       " + seed);
            Console.WriteLine("  Spatial dimension =        " + dim_num);
            Console.WriteLine("  Number of points =         " + point_num);
            Console.WriteLine("  Duplication factor =       " + duplication);
            Console.WriteLine("  Desired minimum distance = " + opt);
            //
            //  Get the points.
            //
            int[] x = IHS.ihs( dim_num, point_num, duplication, ref seed );
            //
            //  Compute the covariance.
            //
            covariance c = new( dim_num, point_num, x );

            Console.WriteLine();
            Console.WriteLine("  Average minimum distance " + c.average);
            Console.WriteLine("  Standard deviation:      " + c.std);
            Console.WriteLine("  Covariance:              " + c.covc);

            typeMethods.i4mat_transpose_print ( dim_num, point_num, x, "  X:" );
        }

    }

    private static void test02 ()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests the improved distributed hypercube sampling algorithm.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 April 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int dim_num = 2;
        const int point_num = 10;

        Console.WriteLine();
        Console.WriteLine("TEST02");
        Console.WriteLine("  IHS implements the IHS Algorithm");
        Console.WriteLine("  (Improved Distributed Hypercube Sampling)");
        Console.WriteLine();
        Console.WriteLine("  Demonstrate the code for a fixed number of points");
        Console.WriteLine("  and dimension, but vary the duplication value.");

        double opt = point_num /
                     Math.Pow ( point_num,
                         1.0E+00 / dim_num );

        Console.WriteLine();
        Console.WriteLine("  Spatial dimension =        " + dim_num);
        Console.WriteLine("  Number of points =         " + point_num);
        Console.WriteLine("  Desired minimum distance = " + opt);

        for (int  duplication = 1; duplication <= 5; duplication++ )
        {
            int seed = 17;

            Console.WriteLine();
            Console.WriteLine("  Random number seed =       " + seed);
            Console.WriteLine("  Duplication factor =       " + duplication);
            //
            //  Get the points.
            //
            int[] x = IHS.ihs( dim_num, point_num, duplication, ref seed );
            //
            //  Compute the covariance.
            //
            covariance c = new( dim_num, point_num, x );

            Console.WriteLine();
            Console.WriteLine("  Average minimum distance " + c.average);
            Console.WriteLine("  Standard deviation:      " + c.std);
            Console.WriteLine("  Covariance:              " + c.covc);

            typeMethods.i4mat_transpose_print ( dim_num, point_num, x, "  X:" );
        }
    }

    private static void test03 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests the improved distributed hypercube sampling algorithm.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 April 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int dim_num = 2;
        const int duplication = 5;

        Console.WriteLine();
        Console.WriteLine("TEST03");
        Console.WriteLine("  IHS implements the IHS Algorithm");
        Console.WriteLine("  (Improved Distributed Hypercube Sampling)");
        Console.WriteLine();
        Console.WriteLine("  Demonstrate the code for a fixed dimension");
        Console.WriteLine("  and duplication value, and increasing number of points.");

        Console.WriteLine();
        Console.WriteLine("  Spatial dimension =        " + dim_num);
        Console.WriteLine("  Duplication factor =       " + duplication);

        int point_num = 5;

        for (int k = 1; k <= 5; k++ )
        {
            point_num = 2 * point_num;

            double opt = point_num /
                         Math.Pow ( point_num,
                             1.0E+00 / dim_num );

            int seed = 17;

            Console.WriteLine();
            Console.WriteLine("  Random number seed =       " + seed);
            Console.WriteLine("  Number of points =         " + point_num);
            Console.WriteLine("  Desired minimum distance = " + opt);
            //
            //  Get the points.
            //
            int[] x = IHS.ihs( dim_num, point_num, duplication, ref seed );
            //
            //  Compute the covariance.
            //
            covariance c = new( dim_num, point_num, x );

            Console.WriteLine();
            Console.WriteLine("  Average minimum distance " + c.average);
            Console.WriteLine("  Standard deviation:      " + c.std);
            Console.WriteLine("  Covariance:              " + c.covc);

            Console.WriteLine();

            for (int j = 0; j < point_num; j++ )
            {
                if (j <= 10 || point_num-10 <= j )
                {
                    string t = j+1 + "    ";
                    string line = t.PadLeft(4);
                    for (int i = 0; i < dim_num; i++ )
                    {
                        t = x[i+j*dim_num] + "  ";
                        line += t.PadLeft(4);
                    }
                    Console.WriteLine(line);
                    Console.WriteLine();
                }
                else
                {
                    switch (j)
                    {
                        case 11:
                            Console.WriteLine("....    ........");
                            break;
                    }
                }
            }
        }
    }

    private static void test04 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests the improved distributed hypercube sampling algorithm.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 April 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int dim_num = 17;
        const int point_num = 1000;

        const int duplication = 5;

        Console.WriteLine();
        Console.WriteLine("TEST04");
        Console.WriteLine("  IHS implements the IHS Algorithm");
        Console.WriteLine("  (Improved Distributed Hypercube Sampling)");
        Console.WriteLine();
        Console.WriteLine("  Demonstrate the code for a fixed number of points,");
        Console.WriteLine("  dimension, and duplication factor, but with a");
        Console.WriteLine("  varying random number seed.");

        double opt = point_num /
                     Math.Pow ( point_num,
                         1.0 / dim_num );

        Console.WriteLine();
        Console.WriteLine("  Spatial dimension =        " + dim_num);
        Console.WriteLine("  Number of points =         " + point_num);
        Console.WriteLine("  Duplication factor =       " + duplication);
        Console.WriteLine("  Desired minimum distance = " + opt);

        int seed = 17;

        for (int k = 1; k <= 4; k++ )
        {
            Console.WriteLine();
            Console.WriteLine("  Random number seed =       " + seed);
            //
            //  Get the points.
            //
            int[] x = IHS.ihs( dim_num, point_num, duplication, ref seed );
            //
            //  Compute the covariance.
            //
            covariance c = new( dim_num, point_num, x );

            Console.WriteLine();
            Console.WriteLine("  Average minimum distance " + c.average);
            Console.WriteLine("  Standard deviation:      " + c.std);
            Console.WriteLine("  Covariance:              " + c.covc);

            typeMethods.i4mat_transpose_print ( dim_num, point_num, x, "  X:" );
        }
    }
}