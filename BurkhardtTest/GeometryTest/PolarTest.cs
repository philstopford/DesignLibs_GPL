using System;
using System.Globalization;
using Burkardt;
using Burkardt.Geometry;
using Burkardt.Uniform;

namespace GeometryTest;

public static class PolarTest
{
    public static void test0685 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0685 tests POLAR_TO_XY and XY_TO_POLAR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int TEST_NUM = 10;

        double r = 0;
        double t = 0;
        int test;
        double[] xy1 = new double[2];
        double[] xy2 = new double[2];

        Console.WriteLine("");
        Console.WriteLine("TEST0685");
        Console.WriteLine("  POLAR_TO_XY converts (R,Theta) to (X,Y);");
        Console.WriteLine("  XY_TO_POLAR converts (X,Y) to (R,Theta).");
        Console.WriteLine("");
        Console.WriteLine("         X           Y     ===>  R           T   =>      X           Y");
        Console.WriteLine("");

        const double b = -1.0;
        const double c = +1.0;
        int seed = 123456789;

        for ( test = 1; test <= TEST_NUM; test++ )
        {
            xy1[0] = UniformRNG.r8_uniform_ab ( b, c, ref seed );
            xy1[1] = UniformRNG.r8_uniform_ab ( b, c, ref seed );

            XY.xy_to_polar ( xy1, ref r, ref t );
            Helpers.polar_to_xy ( r, t, ref xy2 );

            Console.WriteLine("  " + xy1[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + xy1[1].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + t.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + xy2[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + xy2[1].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

}