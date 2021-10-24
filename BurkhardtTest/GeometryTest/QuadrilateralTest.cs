using System;
using System.Linq;
using Burkardt.Quadrilateral;
using Burkardt.Types;

namespace GeometryTest
{
    public static class QuadrilateralTest
    {
        public static void test171()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST171 tests QUAD_AREA_2D and QUAD_AREA2_2D;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 May 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;

            double area;
            double[] quad =
            {
                0.0, 0.0,
                1.0, 0.0,
                1.0, 1.0,
                0.0, 1.0
            };

            Console.WriteLine("");
            Console.WriteLine("TEST171");
            Console.WriteLine("  For a quadrilateral in 2D:");
            Console.WriteLine("  QUAD_AREA_2D finds the area;");
            Console.WriteLine("  QUAD_AREA2_2D finds the area;");

            typeMethods.r8mat_transpose_print(DIM_NUM, 4, quad, "  The vertices:");

            area = Geometry.quad_area_2d(quad);

            Console.WriteLine("");
            Console.WriteLine("  QUAD_AREA_2D area is  " + area + "");

            area = Geometry.quad_area2_2d(quad);

            Console.WriteLine("  QUAD_AREA2_2D area is " + area + "");

        }

        public static void test1712()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST1712 tests QUAD_AREA_3D;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 May 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double area;
            double area1;
            double area2;
            int i;
            int j;
            double[] q =
            {
                2.0, 2.0, 0.0,
                0.0, 0.0, 0.0,
                1.0, 1.0, 1.0,
                3.0, 3.0, 1.0
            };
            double[] t = new double[3 * 3];

            Console.WriteLine("");
            Console.WriteLine("TEST1712");
            Console.WriteLine("  For a quadrilateral in 3D:");
            Console.WriteLine("  QUAD_AREA_3D finds the area.");

            typeMethods.r8mat_transpose_print(3, 4, q, "  The vertices:");

            area = Geometry.quad_area_3d(q);

            Console.WriteLine("");
            Console.WriteLine("  QUAD_AREA_3D area is     " + area + "");

            for (j = 0; j < 3; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    t[i + j * 3] = q[i + j * 3];
                }
            }

            area1 = typeMethods.triangle_area_3d(t);
            for (j = 0; j < 2; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    t[i + j * 3] = q[i + (j + 2) * 3];
                }
            }

            for (i = 0; i < 3; i++)
            {
                t[i + 2 * 3] = q[i + 0 * 3];
            }

            area2 = typeMethods.triangle_area_3d(t);
            Console.WriteLine("  Sum of TRIANGLE_AREA_3D: " + area1 + area2 + "");

        }

        public static void test1715()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST1715 tests QUAD_CONTAINS_POINT_2D, QUAD_POINT_DIST_2D, QUAD_POINT_DIST_SIGNED_2D.
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
            int DIM_NUM = 2;
            int TEST_NUM = 7;

            double dist;
            double dist_signed;
            bool inside;
            double[] p;
            double[] p_test =
            {
                0.25, 0.25,
                0.75, 0.25,
                1.00, 1.00,
                11.00, 0.50,
                0.00, 0.50,
                0.50, -10.00,
                2.00, 2.00
            };
            double[] q =
            {
                0.0, 0.0,
                1.0, 0.0,
                1.0, 1.0,
                0.0, 1.0
            };
            int test;

            Console.WriteLine("");
            Console.WriteLine("TEST1715");
            Console.WriteLine("  For a quadrilateral in 2D:");
            Console.WriteLine("  QUAD_AREA_2D finds the area;");
            Console.WriteLine("  QUAD_CONTAINS_POINT_2D tells if a point is inside;");
            Console.WriteLine("  QUAD_POINT_DIST_2D computes the distance.");
            Console.WriteLine("  QUAD_POINT_DIST_SIGNED_2D computes signed distance.");

            typeMethods.r8mat_transpose_print(DIM_NUM, 4, q, "  The vertices:");

            Console.WriteLine("");
            Console.WriteLine("        P        Contains  Dist    Dist");
            Console.WriteLine("                          Signed  Unsigned");
            Console.WriteLine("");

            for (test = 0; test < TEST_NUM; test++)
            {
                p = p_test.Skip(+test * DIM_NUM).ToArray();

                inside = Geometry.quad_contains_point_2d(q, p);

                dist_signed = Geometry.quad_point_dist_signed_2d(q, p);

                dist = Geometry.quad_point_dist_2d(q, p);

                Console.WriteLine("  " + p[0].ToString().PadLeft(12)
                                       + "  " + p[1].ToString().PadLeft(12)
                                       + "  " + inside.ToString().PadLeft(1)
                                       + "  " + dist_signed.ToString().PadLeft(12)
                                       + "  " + dist.ToString().PadLeft(12) + "");
            }
        }

    }
}