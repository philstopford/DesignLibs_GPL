using System;
using System.Linq;
using Burkardt;
using Burkardt.Geometry;

namespace GeometryTest;

public static class RadecTest
{
    public static void test173()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST173 tests RADEC_DISTANCE_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 6;

        double dec1 = 0;
        double dec2 = 0;
        double[] p1;
        double[] p2;
        double[] p_test =
        {
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
            1.0, 1.0, 1.0,
            5.0, -2.0, -1.0,
            -2.0, -2.0, -2.0
        };
        double ra1 = 0;
        double ra2 = 0;
        int test1;
        int test2;
        double theta;
        double theta_deg;

        Console.WriteLine("");
        Console.WriteLine("TEST173");
        Console.WriteLine("  RADEC_DISTANCE_3D computes the angular separation");
        Console.WriteLine("  between two points on a sphere described in terms of");
        Console.WriteLine("  right ascension and declination.");
        Console.WriteLine("");
        Console.WriteLine("     RA1       DEC1      RA2       DEC2    Radians  Degrees");
        Console.WriteLine("");

        for (test1 = 0; test1 < TEST_NUM; test1++)
        {
            p1 = p_test.Skip(+test1 * DIM_NUM).ToArray();

            XY.xyz_to_radec(p1, ref ra1, ref dec1);

            for (test2 = test1 + 1; test2 < TEST_NUM; test2++)
            {
                p2 = p_test.Skip(+test2 * DIM_NUM).ToArray();

                XY.xyz_to_radec(p2, ref ra2, ref dec2);
                theta = Radec.radec_distance_3d(ra1, dec1, ra2, dec2);
                theta_deg = Helpers.radians_to_degrees(theta);

                Console.WriteLine("  " + ra1.ToString().PadLeft(8)
                                       + "  " + dec1.ToString().PadLeft(8)
                                       + "  " + ra2.ToString().PadLeft(8)
                                       + "  " + dec2.ToString().PadLeft(8)
                                       + "  " + theta.ToString().PadLeft(8)
                                       + "  " + theta_deg.ToString().PadLeft(8) + "");
            }
        }
    }

    public static void test174()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST174 tests RADEC_TO_XYZ and XYZ_TO_RADEC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 6;

        double dec = 0;
        int i;
        double[] p1;
        double[] p2;
        double[] p_test =
        {
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
            1.0, 1.0, 1.0,
            5.0, -2.0, -1.0,
            -2.0, -2.0, -2.0
        };
        double ra = 0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST174");
        Console.WriteLine("  RADEC_TO_XYZ converts XYZ to RADEC coordinates.");
        Console.WriteLine("  XYZ_TO_RADEC converts RADEC to XYZ coordinates.");
        Console.WriteLine("");
        Console.WriteLine("          P1          RA     DEC           P2");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p1 = p_test.Skip(+test * DIM_NUM).ToArray();

            string cout = "";

            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + p1[i].ToString().PadLeft(7);
            }

            XY.xyz_to_radec(p1, ref ra, ref dec);

            cout += ra.ToString().PadLeft(7);
            cout += dec.ToString().PadLeft(7);

            p2 = Radec.radec_to_xyz(ra, dec);

            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + p2[i].ToString().PadLeft(7);
            }

            Console.WriteLine(cout);

        }

    }

}