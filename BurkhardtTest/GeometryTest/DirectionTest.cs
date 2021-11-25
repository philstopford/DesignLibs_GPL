﻿using System;
using System.Globalization;
using Burkardt.Geometry;

namespace GeometryTest;

public static class DirectionTest
{
    public static void test021()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST021 tests DIRECTION_PERT_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int itest;
        double[] vbase = {1.0, 0.0, 0.0};

        Console.WriteLine("");
        Console.WriteLine("TEST021");
        Console.WriteLine("  DIRECTION_PERT_3D perturbs a direction vector.");
        Console.WriteLine("");

        int seed = entropyRNG.RNG.nextint();

        Console.WriteLine("");
        Console.WriteLine("  We use SEED = " + seed + "");

        Console.WriteLine("");
        Console.WriteLine("  Base vector:");
        Console.WriteLine("  " + vbase[0]
                               + "  " + vbase[1]
                               + "  " + vbase[2] + "");

        for (itest = 0; itest < 3; itest++)
        {
            double sigma = itest switch
            {
                0 => 0.99,
                1 => 0.5,
                _ => 0.1
            };

            Console.WriteLine("");
            Console.WriteLine("  Using sigma = " + sigma + "");
            Console.WriteLine("");

            int i;
            for (i = 0; i < 20; i++)
            {
                double[] vran = Direction.direction_pert_3d(sigma, vbase, ref seed);
                Console.WriteLine("  " + vran[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + vran[1].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + vran[2].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            }

        }

    }

    public static void test022()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST022 tests DIRECTION_UNIFORM_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST022");
        Console.WriteLine("  DIRECTION_UNIFORM_3D picks a random direction vector.");

        int seed = entropyRNG.RNG.nextint();

        Console.WriteLine("");
        Console.WriteLine("  We use SEED = " + seed + "");

        Console.WriteLine("");

        for (i = 0; i < 10; i++)
        {
            double[] vran = Direction.direction_uniform_3d(ref seed);
            Console.WriteLine("  " + vran[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + vran[1].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + vran[2].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

    }

    public static void direction_uniform_nd_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRECTION_UNIFORM_ND_TEST tests DIRECTION_UNIFORM_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 4;

        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("DIRECTION_UNIFORM_ND_TEST");
        Console.WriteLine("  DIRECTION_UNIFORM_ND picks a random direction vector.");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double[] vran = Direction.direction_uniform_nd(DIM_NUM, ref seed);
            string cout = "";
            int j;
            for (j = 0; j < DIM_NUM; j++)
            {
                cout += "  " + vran[j].ToString(CultureInfo.InvariantCulture).PadLeft(8);
            }

            Console.WriteLine(cout);
        }

    }

}