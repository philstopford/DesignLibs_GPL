using System;
using System.Globalization;
using Burkardt.Composition;
using Burkardt.SphereNS;
using Burkardt.Table;
using Burkardt.Types;

namespace SphereExactnessTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SPHERE_EXACTNESS.
        //
        //  Discussion:
        //
        //    This program investigates the polynomial exactness of a quadrature
        //    rule for the unit sphere
        //
        //  Usage:
        //
        //    sphere_exactness files prefix degree_max
        //
        //    where
        //
        //    * files explains how the quadrature rule is stored:
        //      'XYZW'  for file 'prefix.xyzw' containing (X,Y,Z,Weight);
        //      'RTPW'  for file 'prefix.rtpw' containing  (Theta, Phi, Weight) (radians);
        //      'DTPW'  for file 'prefix.dtpw' containing  (Theta, Phi, Weight) (degrees);
        //      'XYZ+W' for file 'prefix.xyz' containing (X,Y,Z)
        //              and file 'prefix.w' containing Weight;
        //      'RTP+W' for file 'prefix.rtp' containing (Theta, Phi ) in radians,
        //              and file 'prefix.w' containing Weight;
        //      'DTP+W' for file 'prefix.dtp' containing (Theta, Phi ) in degrees,
        //              and file 'prefix.w' containing Weight;
        //      'XYZ1'  for file 'prefix.xyz' containing (X,Y,Z), 
        //              and equal weights, which do not need to be read in.
        //      'RTP1'  for file 'prefix.rtp' containing (Theta, Phi ) in radians,
        //              and equal weights, which do not need to be read in.
        //      'DTP1'  for file 'prefix.dtp' containing (Theta, Phi ) in degrees,'
        //              and equal weights, which do not need to be read in.
        //    * prefix is the common file prefix;
        //    * degree_max is the maximum monomial degree to check.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int degree;
        int degree_max;
        int dim_num;
        double[] dtp;
        string filename;
        string files;
        int i;
        int j;
        int point_num;
        string prefix;
        double[] rtp;
        double[] w;
        double[] xyz;

        Console.WriteLine("");
        Console.WriteLine("SPHERE_EXACTNESS");
        Console.WriteLine("");
        Console.WriteLine("  Investigate the polynomial exactness of a quadrature");
        Console.WriteLine("  rule for the unit sphere by integrating all monomials");
        Console.WriteLine("  of a given degree.");
        //
        //  Get the file structure;
        //
        try
        {
            files = args[0];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("SPHERE_EXACTNESS:");
            Console.WriteLine("  Describe the files to be read:");
            Console.WriteLine("");
            Console.WriteLine("  For coordinates and weights in one file:");
            Console.WriteLine("    XYZW     (X,Y,Z,Weight)");
            Console.WriteLine("    RTPW     (Theta, Phi, Weight) (radians)");
            Console.WriteLine("    DTPW     (Theta, Phi, Weight) (degrees)");
            Console.WriteLine("  For coordinates in one file and weights in another:");
            Console.WriteLine("    XYZ+W    (X,Y,Z)       + Weight");
            Console.WriteLine("    RTP+W    (Theta, Phi ) + Weight");
            Console.WriteLine("    DTP+W    (Theta, Phi ) + Weight");
            Console.WriteLine("  For coordinates in one file, and equal weights:");
            Console.WriteLine("    XYZ1     (X,Y,Z)");
            Console.WriteLine("    RTP1     (Theta, Phi ) (radians)");
            Console.WriteLine("    DTP1     (Theta, Phi ) (degrees)");

            files = Console.ReadLine();
        }

        //
        //  Get the file prefix.
        //
        try
        {
            prefix = args[1];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("SPHERE_EXACTNESS:");
            Console.WriteLine("  Enter the filename prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Get the maximum degree.
        //
        try
        {
            degree_max = Convert.ToInt32(args[2]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("SPHERE_EXACTNESS:");
            Console.WriteLine("  Please enter the maximum total degree to check.");

            degree_max = Convert.ToInt32(Console.ReadLine());
        }

        //
        //  Summarize the input.
        //
        Console.WriteLine("");
        Console.WriteLine("SPHERE_EXACTNESS: User input:");
        Console.WriteLine("  File structure = \"" + files + "\".");
        Console.WriteLine("  Filename prefix = \"" + prefix + "\".");
        Console.WriteLine("  Maximum total degree to check = " + degree_max + "");
        //
        //  Read data needed to create XYZ and W arrays.
        //
        if (typeMethods.s_eqi(files, "xyzw"))
        {
            filename = prefix + ".xyzw";

            TableHeader hd = typeMethods.r8mat_header_read(filename);
            dim_num = hd.m;
            point_num = hd.n;

            double[] xyzw = typeMethods.r8mat_data_read(filename, 4, point_num);

            xyz = new double[3 * point_num];
            w = new double[point_num];

            for (j = 0; j < point_num; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    xyz[i + j * 3] = xyzw[i + j * 4];
                }

                w[j] = xyzw[3 + j * 4];
            }
        }
        else if (typeMethods.s_eqi(files, "xyz+w"))
        {
            filename = prefix + ".xyz";

            TableHeader hd = typeMethods.r8mat_header_read(filename);
            dim_num = hd.m;
            point_num = hd.n;

            xyz = typeMethods.r8mat_data_read(filename, 3, point_num);

            filename = prefix + ".w";

            w = typeMethods.r8mat_data_read(filename, 1, point_num);
        }
        else if (typeMethods.s_eqi(files, "xyz1"))
        {
            filename = prefix + ".xyz";

            TableHeader hd = typeMethods.r8mat_header_read(filename);
            dim_num = hd.m;
            point_num = hd.n;

            xyz = typeMethods.r8mat_data_read(filename, 3, point_num);

            w = new double[point_num];
            for (i = 0; i < point_num; i++)
            {
                w[i] = 4.0 * Math.PI / point_num;
            }
        }
        else if (typeMethods.s_eqi(files, "rtpw"))
        {
            filename = prefix + ".rtpw";

            TableHeader hd = typeMethods.r8mat_header_read(filename);
            dim_num = hd.m;
            point_num = hd.n;

            double[] rtpw = typeMethods.r8mat_data_read(filename, 3, point_num);

            xyz = new double[3 * point_num];
            w = new double[point_num];

            for (j = 0; j < point_num; j++)
            {
                xyz[0 + j * 3] = Math.Cos(rtpw[0 + j * 3]) * Math.Sin(rtpw[1 + j * 3]);
                xyz[1 + j * 3] = Math.Sin(rtpw[0 + j * 3]) * Math.Sin(rtpw[1 + j * 3]);
                xyz[2 + j * 3] = Math.Cos(rtpw[1 + j * 3]);
                w[j] = rtpw[2 + j * 3];
            }
        }
        else if (typeMethods.s_eqi(files, "rtp+w"))
        {
            filename = prefix + ".rtp";

            TableHeader hd = typeMethods.r8mat_header_read(filename);
            dim_num = hd.m;
            point_num = hd.n;

            rtp = typeMethods.r8mat_data_read(filename, 2, point_num);

            filename = prefix + ".w";

            w = typeMethods.r8mat_data_read(filename, 1, point_num);

            xyz = new double[3 * point_num];

            for (j = 0; j < point_num; j++)
            {
                xyz[0 + j * 3] = Math.Cos(rtp[0 + j * 2]) * Math.Sin(rtp[1 + j * 2]);
                xyz[1 + j * 3] = Math.Sin(rtp[0 + j * 2]) * Math.Sin(rtp[1 + j * 2]);
                xyz[2 + j * 3] = Math.Cos(rtp[1 + j * 2]);
            }
        }
        else if (typeMethods.s_eqi(files, "rtp1"))
        {
            filename = prefix + ".rtp";

            TableHeader hd = typeMethods.r8mat_header_read(filename);
            dim_num = hd.m;
            point_num = hd.n;

            rtp = typeMethods.r8mat_data_read(filename, 2, point_num);

            xyz = new double[3 * point_num];

            for (j = 0; j < point_num; j++)
            {
                xyz[0 + j * 3] = Math.Cos(rtp[0 + j * 2]) * Math.Sin(rtp[1 + j * 2]);
                xyz[1 + j * 3] = Math.Sin(rtp[0 + j * 2]) * Math.Sin(rtp[1 + j * 2]);
                xyz[2 + j * 3] = Math.Cos(rtp[1 + j * 2]);
            }

            w = new double[point_num];
            for (i = 0; i < point_num; i++)
            {
                w[i] = 4.0 * Math.PI / point_num;
            }
        }
        else if (typeMethods.s_eqi(files, "dtpw"))
        {
            filename = prefix + ".dtpw";

            TableHeader hd = typeMethods.r8mat_header_read(filename);

            dim_num = hd.m;
            point_num = hd.n;

            double[] dtpw = typeMethods.r8mat_data_read(filename, 3, point_num);

            xyz = new double[3 * point_num];
            w = new double[point_num];

            for (j = 0; j < point_num; j++)
            {
                dtpw[0 + j * 3] = dtpw[0 + j * 3] * Math.PI / 180.0;
                dtpw[1 + j * 3] = dtpw[1 + j * 3] * Math.PI / 180.0;
            }

            for (j = 0; j < point_num; j++)
            {
                xyz[0 + j * 3] = Math.Cos(dtpw[0 + j * 3]) * Math.Sin(dtpw[1 + j * 3]);
                xyz[1 + j * 3] = Math.Sin(dtpw[0 + j * 3]) * Math.Sin(dtpw[1 + j * 3]);
                xyz[2 + j * 3] = Math.Cos(dtpw[1 + j * 3]);
                w[j] = dtpw[2 + j * 3];
            }

        }
        else if (typeMethods.s_eqi(files, "dtp+w"))
        {
            filename = prefix + ".dtp";

            TableHeader hd = typeMethods.r8mat_header_read(filename);
            dim_num = hd.m;
            point_num = hd.n;

            dtp = typeMethods.r8mat_data_read(filename, 2, point_num);

            filename = prefix + ".w";

            w = typeMethods.r8mat_data_read(filename, 1, point_num);

            xyz = new double[3 * point_num];

            for (j = 0; j < point_num; j++)
            {
                dtp[0 + j * 2] = dtp[0 + j * 2] * Math.PI / 180.0;
                dtp[1 + j * 2] = dtp[1 + j * 2] * Math.PI / 180.0;
            }

            for (j = 0; j < point_num; j++)
            {
                xyz[0 + j * 3] = Math.Cos(dtp[0 + j * 2]) * Math.Sin(dtp[1 + j * 2]);
                xyz[1 + j * 3] = Math.Sin(dtp[0 + j * 2]) * Math.Sin(dtp[1 + j * 2]);
                xyz[2 + j * 3] = Math.Cos(dtp[1 + j * 2]);
            }
        }
        else if (typeMethods.s_eqi(files, "dtp1"))
        {
            filename = prefix + ".dtp";

            TableHeader hd = typeMethods.r8mat_header_read(filename);
            dim_num = hd.m;
            point_num = hd.n;

            dtp = typeMethods.r8mat_data_read(filename, 2, point_num);

            xyz = new double[3 * point_num];

            for (j = 0; j < point_num; j++)
            {
                dtp[0 + j * 2] = dtp[0 + j * 2] * Math.PI / 180.0;
                dtp[1 + j * 2] = dtp[1 + j * 2] * Math.PI / 180.0;
            }

            for (j = 0; j < point_num; j++)
            {
                xyz[0 + j * 3] = Math.Cos(dtp[0 + j * 2]) * Math.Sin(dtp[1 + j * 2]);
                xyz[1 + j * 3] = Math.Sin(dtp[0 + j * 2]) * Math.Sin(dtp[1 + j * 2]);
                xyz[2 + j * 3] = Math.Cos(dtp[1 + j * 2]);
            }

            w = new double[point_num];
            for (i = 0; i < point_num; i++)
            {
                w[i] = 4.0 * Math.PI / point_num;
            }
        }
        else
        {
            Console.WriteLine("");
            Console.WriteLine("SPHERE_EXACTNESS - Fatal error!");
            Console.WriteLine("  Unrecognized file structure choice");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of points  = " + point_num + "");
        //
        //  The W's should sum to 4 * PI.
        //
        double w_sum = typeMethods.r8vec_sum(point_num, w);

        for (i = 0; i < point_num; i++)
        {
            w[i] = 4.0 * Math.PI * w[i] / w_sum;
        }

        //
        //  Explore the monomials.
        //
        int[] expon = new int[dim_num];

        Console.WriteLine("");
        Console.WriteLine("      Error    Degree  Exponents");

        for (degree = 0; degree <= degree_max; degree++)
        {
            Console.WriteLine("");
            bool more = false;
            int h = 0;
            int t = 0;

            for (;;)
            {
                Comp.comp_next(degree, dim_num, ref expon, ref more, ref h, ref t);

                double quad_error = Sphere.sphere01_monomial_quadrature(expon, point_num, xyz, w);

                string cout = "  " + quad_error.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "     " + degree.ToString().PadLeft(2)
                                   + "  ";

                int dim;
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += expon[dim].ToString().PadLeft(3);
                }

                Console.WriteLine(cout);

                if (!more)
                {
                    break;
                }
            }

        }

        Console.WriteLine("");
        Console.WriteLine("SPHERE_EXACTNESS:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}