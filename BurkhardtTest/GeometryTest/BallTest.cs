﻿using System;
using System.Globalization;
using Burkardt.Ball;
using Burkardt.Types;

namespace GeometryTest;

public static class BallTest
{
    public static void ball01_sample_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BALL01_SAMPLE_2D_TEST tests BALL01_SAMPLE_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 2;

        double[] average = new double[DIM_NUM];
        int i;
        int j;
        const int sample_num = 1000;
        int seed = 123456789;
        double[] x;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("BALL01_SAMPLE_2D_TEST");
        Console.WriteLine("  For the unit ball in 2 dimensions (the disk):");
        Console.WriteLine("  BALL01_SAMPLE_2D samples;");

        Console.WriteLine("");
        Console.WriteLine("  A few sample values:");
        Console.WriteLine("");

        for (i = 1; i <= 5; i++)
        {
            cout = "";
            x = Geometry.ball01_sample_2d(ref seed);
            for (j = 0; j < DIM_NUM; j++)
            {
                cout += "  " + x[j].ToString(CultureInfo.InvariantCulture).PadLeft(10);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of sample points = " + sample_num + "");

        typeMethods.r8vec_zero(DIM_NUM, ref average);

        for (i = 1; i <= sample_num; i++)
        {
            x = Geometry.ball01_sample_2d(ref seed);
            for (j = 0; j < DIM_NUM; j++)
            {
                average[j] += x[j];
            }
        }

        for (j = 0; j < DIM_NUM; j++)
        {
            average[j] /= sample_num;
        }

        Console.WriteLine("");
        Console.WriteLine("  Now average the points, which should get a value");
        Console.WriteLine("  close to zero, and closer as N increases.");
        Console.WriteLine("");
        cout = "  Average:        ";
        for (j = 0; j < DIM_NUM; j++)
        {
            cout += "  " + average[j].ToString(CultureInfo.InvariantCulture).PadLeft(10);
        }

        Console.WriteLine(cout);
        cout = "";

        double average_r = 0.0;

        for (i = 1; i <= sample_num; i++)
        {
            x = Geometry.ball01_sample_2d(ref seed);
            double temp = 0.0;
            for (j = 0; j < DIM_NUM; j++)
            {
                temp += x[j] * x[j];
            }

            average_r += Math.Sqrt(temp);
        }

        average_r /= sample_num;

        Console.WriteLine("");
        Console.WriteLine("  Now average the distance of the points from the center,");
        Console.WriteLine("  which should be 1/sqrt(2) = "
                          + 1.0 / Math.Sqrt(2.0) + "");
        Console.WriteLine("");
        Console.WriteLine("  Average:        " + average_r + "");

        double average_theta = 0.0;

        for (i = 1; i <= sample_num; i++)
        {
            x = Geometry.ball01_sample_2d(ref seed);
            double theta = typeMethods.r8_atan(x[1], x[0]);
            average_theta += theta;
        }

        average_theta /= sample_num;

        Console.WriteLine("");
        Console.WriteLine("  Now average the angle THETA,");
        Console.WriteLine("  which should be PI.");
        Console.WriteLine("");
        Console.WriteLine("  Average:        " + average_theta + "");

    }

    public static void ball01_sample_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BALL01_SAMPLE_3D_TEST tests BALL01_SAMPLE_3D.
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
        const int DIM_NUM = 3;

        double[] average = new double[DIM_NUM];
        int i;
        int j;
        const int sample_num = 1000;
        double r;
        int seed = 123456789;
        double[] x;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("BALL01_SAMPLE_3D_TEST");
        Console.WriteLine("  For the unit ball in 3 dimensions:");
        Console.WriteLine("  BALL01_SAMPLE_3D samples;");

        Console.WriteLine("");
        Console.WriteLine("  A few sample values:");
        Console.WriteLine("");

        for (i = 1; i <= 5; i++)
        {
            cout = "";
            x = Geometry.ball01_sample_3d(ref seed);
            for (j = 0; j < DIM_NUM; j++)
            {
                cout += "  " + x[j].ToString(CultureInfo.InvariantCulture).PadLeft(8);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of sample points = " + sample_num + "");

        typeMethods.r8vec_zero(DIM_NUM, ref average);

        for (i = 1; i <= sample_num; i++)
        {
            x = Geometry.ball01_sample_3d(ref seed);
            for (j = 0; j < DIM_NUM; j++)
            {
                average[j] += x[j];
            }
        }

        for (j = 0; j < DIM_NUM; j++)
        {
            average[j] /= sample_num;
        }

        Console.WriteLine("");
        Console.WriteLine("  Now average the points, which should get a value");
        Console.WriteLine("  close to zero, and closer as sample_num increases.");
        Console.WriteLine("");
        cout = "  Average:";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + average[i].ToString(CultureInfo.InvariantCulture).PadLeft(8);
        }

        Console.WriteLine(cout);
        cout = "";

        seed = 123456789;

        double average_r = 0.0;

        for (i = 1; i <= sample_num; i++)
        {
            x = Geometry.ball01_sample_3d(ref seed);
            r = 0.0;
            for (j = 0; j < DIM_NUM; j++)
            {
                r += x[j] * x[j];
            }

            r = Math.Sqrt(r);
            average_r += r;
        }

        average_r /= sample_num;

        Console.WriteLine("");
        Console.WriteLine("  Now average the distance of the points from");
        Console.WriteLine("  the center, which should be the ");
        Console.WriteLine("  1/2^(1/n) = " + Math.Pow(0.5, 1.0 / DIM_NUM) + "");
        Console.WriteLine("");
        Console.WriteLine("  Average:        " + average_r + "");

        double average_theta = 0.0;

        for (i = 1; i <= sample_num; i++)
        {
            x = Geometry.ball01_sample_3d(ref seed);
            double theta = typeMethods.r8_atan(x[1], x[0]);
            average_theta += theta;
        }

        average_theta /= sample_num;

        Console.WriteLine("");
        Console.WriteLine("  Now average the angle THETA,");
        Console.WriteLine("  which should be PI.");
        Console.WriteLine("");
        Console.WriteLine("  Average:        " + average_theta + "");

        double average_phi = 0.0;

        for (i = 1; i <= sample_num; i++)
        {
            x = Geometry.ball01_sample_3d(ref seed);
            r = 0.0;
            for (j = 0; j < DIM_NUM; j++)
            {
                r += x[j] * x[j];
            }

            r = Math.Sqrt(r);
            double phi = Math.Acos(x[2] / r);
            average_phi += phi;
        }

        average_phi /= sample_num;

        Console.WriteLine("");
        Console.WriteLine("  Now average the angle PHI,");
        Console.WriteLine("  which should be PI/2.");
        Console.WriteLine("");
        Console.WriteLine("  Average:        " + average_phi + "");

    }

    public static void ball01_sample_nd_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BALL01_SAMPLE_ND_TEST tests BALL01_SAMPLE_ND.
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
        const int DIM_NUM = 3;

        double[] average = new double[DIM_NUM];
        int i;
        int j;
        const int sample_num = 1000;
        double r;
        int seed = 123456789;
        double[] x;
        string cout;

        Console.WriteLine("");
        Console.WriteLine("BALL01_SAMPLE_ND_TEST");
        Console.WriteLine("  For the unit ball in N dimensions:");
        Console.WriteLine("  BALL01_SAMPLE_ND samples;");

        Console.WriteLine("");
        Console.WriteLine("  A few sample values:");
        Console.WriteLine("");

        for (i = 1; i <= 5; i++)
        {
            cout = "";
            x = Geometry.ball01_sample_nd(DIM_NUM, ref seed);
            for (j = 0; j < DIM_NUM; j++)
            {
                cout += "  " + x[j].ToString(CultureInfo.InvariantCulture).PadLeft(8);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of sample points = " + sample_num + "");

        typeMethods.r8vec_zero(DIM_NUM, ref average);

        for (i = 1; i <= sample_num; i++)
        {
            x = Geometry.ball01_sample_nd(DIM_NUM, ref seed);
            for (j = 0; j < DIM_NUM; j++)
            {
                average[j] += x[j];
            }
        }

        for (j = 0; j < DIM_NUM; j++)
        {
            average[j] /= sample_num;
        }

        Console.WriteLine("");
        Console.WriteLine("  Now average the points, which should get a value");
        Console.WriteLine("  close to zero, and closer as N increases.");
        Console.WriteLine("");
        cout = "  Average:        ";
        for (j = 0; j < DIM_NUM; j++)
        {
            cout += "  " + average[j].ToString(CultureInfo.InvariantCulture).PadLeft(8);
        }

        Console.WriteLine(cout);
        cout = "";

        seed = 123456789;

        double average_r = 0.0;

        for (i = 1; i <= sample_num; i++)
        {
            x = Geometry.ball01_sample_nd(DIM_NUM, ref seed);
            r = 0.0;
            for (j = 0; j < DIM_NUM; j++)
            {
                r += x[j] * x[j];
            }

            r = Math.Sqrt(r);
            average_r += r;
        }

        average_r /= sample_num;

        Console.WriteLine("");
        Console.WriteLine("  Now average the distance of the points from");
        Console.WriteLine("  the center, which should be the");
        Console.WriteLine("  1/2^(1/dim_num) = "
                          + Math.Pow(0.5, 1.0 / DIM_NUM) + "");
        Console.WriteLine("");
        Console.WriteLine("  Average:        " + average_r + "");

        double average_theta = 0.0;

        for (i = 1; i <= sample_num; i++)
        {
            x = Geometry.ball01_sample_nd(DIM_NUM, ref seed);
            double theta = typeMethods.r8_atan(x[1], x[0]);
            average_theta += theta;
        }

        average_theta /= sample_num;

        Console.WriteLine("");
        Console.WriteLine("  Now average the angle THETA,");
        Console.WriteLine("  which should be PI.");
        Console.WriteLine("");
        Console.WriteLine("  Average:        " + average_theta + "");

        double average_phi = 0.0;

        for (i = 1; i <= sample_num; i++)
        {
            x = Geometry.ball01_sample_nd(DIM_NUM, ref seed);
            r = 0.0;
            for (j = 0; j < DIM_NUM; j++)
            {
                r += x[j] * x[j];
            }

            r = Math.Sqrt(r);
            double phi = Math.Acos(x[2] / r);
            average_phi += phi;
        }

        average_phi /= sample_num;

        Console.WriteLine("");
        Console.WriteLine("  Now average the angle PHI,");
        Console.WriteLine("  which should be PI/2.");
        Console.WriteLine("");
        Console.WriteLine("  Average:        " + average_phi + "");

    }

    public static void ball01_volume_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BALL01_VOLUME_TEST tests BALL01_VOLUME.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("BALL01_VOLUME_TEST");
        Console.WriteLine("  BALL01_VOLUME returns the volume of the unit ball.");

        double volume = Geometry.ball01_volume();

        Console.WriteLine("");
        Console.WriteLine("  Volume = " + volume + "");

    }

}