using System;
using Burkardt;
using Burkardt.AnnulusNS;

namespace GeometryTest;

public static class AnnulusTest
{
    public static void annulus_sector_centroid_2d_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ANNULUS_SECTOR_CENTROID_2D_TEST tests ANNULUS_SECTOR_CENTROID_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] centroid;
        double[] pc = { 5.0, 3.0 };
        double r1 = 2.0;
        double r2 = 3.0;
        double theta1;
        double theta2;

        theta1 = Helpers.degrees_to_radians ( 30.0 );
        theta2 = Helpers.degrees_to_radians ( 60.0 );

        Console.WriteLine("");
        Console.WriteLine("ANNULUS_SECTOR_CENTROID_2D_TEST");
        Console.WriteLine("  ANNULUS_SECTOR_CENTROID_2D computes the centroid of a");
        Console.WriteLine("  circular annulus.");
        Console.WriteLine("");
        Console.WriteLine("  The circle has center        " + pc[0]
                                                            + "  " + pc[1] + "");
        Console.WriteLine("  The inner radius is R1 =     " + r1     + "");
        Console.WriteLine("  The outer radius is R2 =     " + r2     + "");
        Console.WriteLine("  The first angle is THETA1 =  " + theta1 + "");
        Console.WriteLine("  The second angle is THETA2 = " + theta2 + "");

        centroid = Annulus.annulus_sector_centroid_2d ( pc, r1, r2, theta1, theta2 );

        Console.WriteLine("");
        Console.WriteLine("  Centroid: "
                          + centroid[0].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                          + centroid[1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            
    }

}