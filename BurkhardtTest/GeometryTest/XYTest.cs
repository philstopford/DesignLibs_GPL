using System;
using Burkardt.Geometry;
using Burkardt.Uniform;

namespace GeometryTest;

public static class XYTest
{
    public static void test1893 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //   TEST1893 tests RTP_TO_XYZ and XYZ_TO_RTP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    221 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = -2.0;
        double b =  3.0;
        double phi = 0;
        double r = 0;
        int seed;
        int test;
        int test_num = 5;
        double theta = 0;
        double[] xyz1;
        double[] xyz2 = new double[3];

        Console.WriteLine("");
        Console.WriteLine("TEST1893");
        Console.WriteLine("  RTP_TO_XYZ converts XYZ to (R,Theta,Phi) coordinates.");
        Console.WriteLine("  XYZ_TO_RTP converts (R,Theta,Phi) to XYZ coordinates.");
        Console.WriteLine("");
        Console.WriteLine("      X1     Y1     Z1     R    THETA    PHI    X2     Y2     Z2");
        Console.WriteLine("");

        seed = 123456789;

        for ( test = 1; test <= test_num; test++ )
        {
            xyz1 = UniformRNG.r8vec_uniform_ab_new ( 3, a, b, ref seed );

            XY.xyz_to_rtp ( xyz1, ref r, ref theta, ref phi );
            XY.rtp_to_xyz ( r, theta, phi, ref xyz2 );

            Console.WriteLine("  " + xyz1[0].ToString().PadLeft(7)
                                   + "  " + xyz1[1].ToString().PadLeft(7)
                                   + "  " + xyz1[1].ToString().PadLeft(7)
                                   + "  " + r.ToString().PadLeft(7)
                                   + "  " + theta.ToString().PadLeft(7)
                                   + "  " + phi.ToString().PadLeft(7)
                                   + "  " + xyz2[0].ToString().PadLeft(7)
                                   + "  " + xyz2[1].ToString().PadLeft(7)
                                   + "  " + xyz2[2].ToString().PadLeft(7) + "");

        }
    }


}