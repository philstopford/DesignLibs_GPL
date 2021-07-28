using System;
using Burkardt.TriangleNS;
using Burkardt.Types;

namespace TrianglePropertiesTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TRIANGLE_PROPERTIES_TEST.
            //
            //  Discussion:
            //
            //    TRIANGLE_PROPERTIES_TEST tests the TRIANGLE_PROPERTIES library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_PROPERTIES_TEST");
            Console.WriteLine("  Test the TRIANGLE_PROPERTIES library.");

            triangle_angles_test();
            triangle_area_test();
            triangle_centroid_test();
            triangle_circumcircle_test();
            triangle_contains_point_test();
            triangle_diameter_test();
            triangle_edge_length_test();
            triangle_incircle_test();
            triangle_orientation_test();
            triangle_orthocenter_test();
            triangle_point_dist_test();
            triangle_point_near_test();
            triangle_quality_test();
            triangle_reference_sample_test();
            triangle_sample_test();
            triangle_xsi_to_xy_test();
            triangle_xy_to_xsi_test();
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_PROPERTIES_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void triangle_angles_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_ANGLES_TEST tests TRIANGLE_ANGLES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] angle;
            int i;
            const double r8_pi = 3.141592653589793;
            double[] t =
            {
                0.0, 1.0,
                0.0, 0.0,
                1.0, 0.0
            };

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_ANGLES_TEST");
            Console.WriteLine("  TRIANGLE_ANGLES computes the angles in a triangle.");

            typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

            angle = typeMethods.triangle_angles_2d_new(t);

            Console.WriteLine("");
            Console.WriteLine("      Radians      Degrees");
            Console.WriteLine("");
            for (i = 0; i < 3; i++)
            {
                Console.WriteLine("  " + angle[i].ToString().PadLeft(14)
                                       + "  " + (angle[i] * 180.0 / r8_pi) + "");
            }
        }

        static void triangle_area_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_AREA_TEST tests TRIANGLE_AREA.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double area;
            double[] t =
            {
                0.0, 1.0,
                0.0, 0.0,
                1.0, 0.0
            };

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_AREA_TEST");
            Console.WriteLine("  TRIANGLE_AREA computes the area of a triangle.");

            typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

            area = Integrals.triangle_area(t);

            Console.WriteLine("");
            Console.WriteLine("  Triangle area is " + area + "");

            return;
        }

        static void triangle_centroid_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_CENTROID_TEST tests TRIANGLE_CENTROID;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] centroid;
            int i;
            int j;
            double[] t = new double[2 * 3];
            double[] t_test =
            {
                0.0, 0.0,
                1.0, 0.0,
                0.0, 1.0,
                0.0, 0.0,
                1.0, 0.0,
                0.5, 0.86602539,
                0.0, 0.0,
                1.0, 0.0,
                0.5, 10.0,
                0.0, 0.0,
                1.0, 0.0,
                10.0, 2.0
            };
            int test;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_CENTROID_TEST");
            Console.WriteLine("  TRIANGLE_CENTROID computes the centroid of a triangle.");

            for (test = 0; test < 4; test++)
            {
                for (j = 0; j < 3; j++)
                {
                    for (i = 0; i < 2; i++)
                    {
                        t[i + j * 2] = t_test[i + j * 2 + test * 2 * 3];
                    }
                }

                typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

                centroid = typeMethods.triangle_centroid_2d(t);

                typeMethods.r8vec_print(2, centroid, "  Centroid:");

            }
        }

        static void triangle_circumcircle_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_CIRCUMCIRCLE_TEST tests TRIANGLE_CIRCUMCIRCLE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int j;
            double[] pc = new double[2];
            double r = 0;
            double[] t = new double[2 * 3];
            double[] t_test =
            {
                0.0, 0.0,
                1.0, 0.0,
                0.0, 1.0,
                0.0, 0.0,
                1.0, 0.0,
                0.5, 0.86602539,
                0.0, 0.0,
                1.0, 0.0,
                0.5, 10.0,
                0.0, 0.0,
                1.0, 0.0,
                10.0, 2.0
            };
            int test;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_CIRCUMCIRCLE_TEST");
            Console.WriteLine("  TRIANGLE_CIRCUMCIRCLE computes the circumcircle of a triangle.");

            for (test = 0; test < 4; test++)
            {
                for (j = 0; j < 3; j++)
                {
                    for (i = 0; i < 2; i++)
                    {
                        t[i + j * 2] = t_test[i + j * 2 + test * 2 * 3];
                    }
                }

                typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

                typeMethods.triangle_circumcircle_2d(t, ref r, ref pc);

                typeMethods.r8vec_print(2, pc, "  Circumcenter");
                Console.WriteLine("  Circumradius: " + r + "");

            }

            return;
        }

        static void triangle_contains_point_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_CONTAINS_POINT_TEST tests TRIANGLE_CONTAINS_POINT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            bool inside;
            int j;
            double[] p = new double[2];
            double[] p_test =
            {
                0.25, 0.25,
                0.75, 0.25,
                1.00, 1.00,
                11.00, 0.50,
                0.00, 1.00,
                0.50, -10.00,
                0.60, 0.60
            };
            double[] t =
            {
                0.0, 1.0,
                0.0, 0.0,
                1.0, 0.0
            };
            double[] t2 = new double[2 * 3];
            int test;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_CONTAINS_POINT_TEST");
            Console.WriteLine("  TRIANGLE_CONTAINS_POINT reports if a point");
            Console.WriteLine("  is inside a triangle.");

            typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

            Console.WriteLine("");
            Console.WriteLine("       X       Y     Inside");
            Console.WriteLine("");

            for (test = 0; test < 7; test++)
            {
                p[0] = p_test[0 + test * 2];
                p[1] = p_test[1 + test * 2];

                inside = typeMethods.triangle_contains_point_2d(t, p);

                Console.WriteLine("  " + p[0].ToString().PadLeft(8)
                                       + "  " + p[1].ToString().PadLeft(8)
                                       + "  " + inside + "");
            }

            //
            //  Make a copy of the triangle with vertices in reverse order.
            //
            Console.WriteLine("");
            Console.WriteLine("  Repeat the test, but reverse the triangle vertex");
            Console.WriteLine("  ordering.");

            for (j = 0; j < 3; j++)
            {
                t2[0 + j * 2] = t[0 + (2 - j) * 3];
                t2[1 + j * 2] = t[1 + (2 - j) * 3];
            }

            typeMethods.r8mat_transpose_print(2, 3, t2, "  Triangle vertices (reversed):");

            Console.WriteLine("");
            Console.WriteLine("       X       Y     Inside");
            Console.WriteLine("");

            for (test = 0; test < 7; test++)
            {
                p[0] = p_test[0 + test * 2];
                p[1] = p_test[1 + test * 2];

                inside = typeMethods.triangle_contains_point_2d(t2, p);

                Console.WriteLine("  " + p[0].ToString().PadLeft(8)
                                       + "  " + p[1].ToString().PadLeft(8)
                                       + "  " + inside + "");
            }

            return;
        }

        static void triangle_diameter_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_DIAMETER_TEST tests TRIANGLE_DIAMETER.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double diameter;
            int i;
            int j;
            double[] t = new double[2 * 3];
            double[] t_test =
            {
                4.0, 2.0,
                1.0, 5.0,
                -2.0, 2.0,
                4.0, 2.0,
                5.0, 4.0,
                6.0, 6.0,
                4.0, 2.0,
                1.0, 5.0,
                4.0, 2.0
            };
            int test;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_DIAMETER_TEST");
            Console.WriteLine("  TRIANGLE_DIAMETER computes the diameter of");
            Console.WriteLine("  the SMALLEST circle around the triangle.");

            for (test = 0; test < 3; test++)
            {
                for (j = 0; j < 3; j++)
                {
                    for (i = 0; i < 2; i++)
                    {
                        t[i + j * 2] = t_test[i + j * 2 + test * 6];
                    }
                }

                typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

                diameter = typeMethods.triangle_diameter_2d(t);

                Console.WriteLine("");
                Console.WriteLine("  Diameter = " + diameter + "");
            }

        }

        static void triangle_edge_length_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_EDGE_LENGTH_TEST tests TRIANGLE_EDGE_LENGTH.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] edge_length;
            int i;
            int j;
            double[] t = new double[2 * 3];
            double[] t_test =
            {
                4.0, 2.0,
                1.0, 5.0,
                -2.0, 2.0,
                4.0, 2.0,
                5.0, 4.0,
                6.0, 6.0,
                4.0, 2.0,
                1.0, 5.0,
                4.0, 2.0
            };
            int test;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_EDGE_LENGTH_TEST");
            Console.WriteLine("  TRIANGLE_EDGE_LENGTH computes the edge lengths of a triangle.");

            for (test = 0; test < 3; test++)
            {
                for (j = 0; j < 3; j++)
                {
                    for (i = 0; i < 2; i++)
                    {
                        t[i + j * 2] = t_test[i + j * 2 + test * 6];
                    }
                }

                typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

                edge_length = typeMethods.triangle_edge_length_2d(t);

                typeMethods.r8vec_print(3, edge_length, "  EDGE_LENGTHS:");

            }
        }

        static void triangle_incircle_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_INCIRCLE_TEST tests TRIANGLE_INCIRCLE;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] pc = new double[2];
            double r = 0;
            double[] t =
            {
                0.0, 1.0,
                0.0, 0.0,
                1.0, 0.0
            };

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_INCIRCLE_TEST");
            Console.WriteLine("  TRIANGLE_INCIRCLE_2D computes the incircle of a triangle.");

            typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

            typeMethods.triangle_incircle_2d(t, ref pc, ref r);

            typeMethods.r8vec_print(2, pc, "  Incenter");

            Console.WriteLine("");
            Console.WriteLine("  Incircle radius is " + r + "");

        }

        static void triangle_orientation_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_ORIENTATION_TEST tests TRIANGLE_ORIENTATION.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int j;
            double[] t = new double[2 * 3];
            double[] t_test =
            {
                4.0, 2.0,
                1.0, 5.0,
                -2.0, 2.0,
                1.0, 5.0,
                4.0, 2.0,
                1.0, -1.0,
                1.0, 5.0,
                2.0, 7.0,
                3.0, 9.0,
                1.0, 5.0,
                4.0, 2.0,
                1.0, 5.0
            };
            int test;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_ORIENTATION_TEST");
            Console.WriteLine("  TRIANGLE_ORIENTATION determines the orientation of a triangle.");

            for (test = 0; test < 4; test++)
            {
                for (j = 0; j < 3; j++)
                {
                    for (i = 0; i < 2; i++)
                    {
                        t[i + j * 2] = t_test[i + j * 2 + test * 6];
                    }
                }

                i = typeMethods.triangle_orientation_2d(t);

                typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

                if (i == 0)
                {
                    Console.WriteLine("  The points are counterclockwise.");
                }
                else if (i == 1)
                {
                    Console.WriteLine("  The points are clockwise.");
                }
                else if (i == 2)
                {
                    Console.WriteLine("  The points are colinear.");
                }
                else if (i == 3)
                {
                    Console.WriteLine("  The points are not distinct.");
                }
                else
                {
                    Console.WriteLine("  The return value makes no sense.");
                }

            }

        }

        static void triangle_orthocenter_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_ORTHOCENTER_TEST tests TRIANGLE_ORTHOCENTER;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            bool flag = false;
            int i;
            int j;
            double[] pc = new double[2];
            double[] t = new double[2 * 3];
            double[] t_test =
            {
                0.0, 0.0,
                1.0, 0.0,
                0.0, 1.0,
                0.0, 0.0,
                1.0, 0.0,
                0.5, 0.86602539,
                0.0, 0.0,
                1.0, 0.0,
                0.5, 10.0,
                0.0, 0.0,
                1.0, 0.0,
                10.0, 2.0
            };
            int test;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_ORTHOCENTER_TEST");
            Console.WriteLine("  TRIANGLE_ORTHOCENTER_2D computes the orthocenter of a triangle.");

            for (test = 0; test < 4; test++)
            {
                for (j = 0; j < 3; j++)
                {
                    for (i = 0; i < 2; i++)
                    {
                        t[i + j * 2] = t_test[i + j * 2 + test * 6];
                    }
                }

                typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

                typeMethods.triangle_orthocenter_2d(t, pc, ref flag);

                typeMethods.r8vec_print(2, pc, "  Orthocenter");
            }

        }

        static void triangle_point_dist_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_POINT_DIST_TEST tests TRIANGLE_POINT_DIST;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double dist;
            double[] p = new double[2];
            double[] p_test =
            {
                0.25, 0.25,
                0.75, 0.25,
                1.00, 1.00,
                11.00, 0.50,
                0.00, 1.00,
                0.50, -10.00,
                0.60, 0.60
            };
            double[] t =
            {
                0.0, 1.0,
                0.0, 0.0,
                1.0, 0.0
            };
            int test;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_POINT_DIST_TEST");
            Console.WriteLine("  TRIANGLE_POINT_DIST computes the distance");
            Console.WriteLine("  from a point to a triangle.");

            typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

            Console.WriteLine("");
            Console.WriteLine("           P            DIST");
            Console.WriteLine("");

            for (test = 0; test < 7; test++)
            {
                p[0] = p_test[0 + test * 2];
                p[1] = p_test[1 + test * 2];

                dist = typeMethods.triangle_point_dist(t, p);

                Console.WriteLine("  " + p[0].ToString().PadLeft(8)
                                       + "  " + p[1].ToString().PadLeft(8)
                                       + "  "
                                       + "  " + dist.ToString().PadLeft(8) + "");
            }

        }

        static void triangle_point_near_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_POINT_NEAR_TEST tests TRIANGLE_POINT_NEAR.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double dist = 0;
            double[] p = new double[2];
            double[] p_test =
            {
                0.25, 0.25,
                0.75, 0.25,
                1.00, 1.00,
                11.00, 0.50,
                0.00, 1.00,
                0.50, -10.00,
                0.60, 0.60
            };
            double[] pn = new double[2];
            double[] t =
            {
                0.0, 1.0,
                0.0, 0.0,
                1.0, 0.0
            };
            int test;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_POINT_NEAR_TEST");
            Console.WriteLine("  TRIANGLE_POINT_NEAR computes the nearest");
            Console.WriteLine("  point on a triangle to a given point.");

            typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

            Console.WriteLine("");
            Console.WriteLine("           P                PN");
            Console.WriteLine("");

            for (test = 0; test < 7; test++)
            {
                p[0] = p_test[0 + test * 2];
                p[1] = p_test[1 + test * 2];


                typeMethods.triangle_point_near(t, p, ref pn, ref dist);

                Console.WriteLine("  " + p[0].ToString().PadLeft(8)
                                       + "  " + p[1].ToString().PadLeft(8)
                                       + "  "
                                       + "  " + pn[0].ToString().PadLeft(8)
                                       + "  " + pn[1].ToString().PadLeft(8) + "");
            }
        }

        static void triangle_quality_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_QUALITY_TEST tests TRIANGLE_QUALITY.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int j;
            double quality;
            double[] t = new double[2 * 3];
            double[] t_test =
            {
                0.0, 0.0,
                1.0, 0.0,
                0.0, 1.0,
                0.0, 0.0,
                1.0, 0.0,
                0.5, 0.86602539,
                0.0, 0.0,
                1.0, 0.0,
                0.5, 10.0,
                0.0, 0.0,
                1.0, 0.0,
                10.0, 2.0
            };
            int test;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_QUALITY_TEST");
            Console.WriteLine("  TRIANGLE_QUALITY computes the quality of a triangle.");

            for (test = 0; test < 4; test++)
            {
                for (j = 0; j < 3; j++)
                {
                    for (i = 0; i < 2; i++)
                    {
                        t[i + j * 2] = t_test[i + j * 2 + test * 6];
                    }
                }

                typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

                quality = typeMethods.triangle_quality_2d(t);

                Console.WriteLine("");
                Console.WriteLine("  Quality = " + quality + "");
            }

        }

        static void triangle_reference_sample_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_REFERENCE_SAMPLE_TEST tests TRIANGLE_REFERENCE_SAMPLE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] p = null;
            int seed = 123456789;
            double[] t =
            {
                0.0, 0.0,
                1.0, 0.0,
                0.0, 1.0
            };
            int test;
            double[] xsi;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_REFERENCE_SAMPLE_TEST");
            Console.WriteLine("  TRIANGLE_REFERENCE_SAMPLE samples the reference triangle.");

            Console.WriteLine("");
            Console.WriteLine("  Sample points (X,Y) and (XSI1,XSI2,XSI3) coordinates:");
            Console.WriteLine("");

            for (test = 0; test < 10; test++)
            {
                typeMethods.triangle_reference_sample(1, ref seed, ref p);
                xsi = typeMethods.triangle_xy_to_xsi(t, p);
                Console.WriteLine("  " + p[0].ToString().PadLeft(8)
                                       + "  " + p[1].ToString().PadLeft(8)
                                       + "  "
                                       + "  " + xsi[0].ToString().PadLeft(8)
                                       + "  " + xsi[1].ToString().PadLeft(8)
                                       + "  " + xsi[2].ToString().PadLeft(8) + "");
            }
        }

        static void triangle_sample_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_SAMPLE_TEST tests TRIANGLE_SAMPLE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] p = null;
            int seed = 123456789;
            double[] t =
            {
                4.0, 2.0,
                1.0, 5.0,
                -2.0, 2.0
            };
            int test;
            double[] xsi;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_SAMPLE_TEST");
            Console.WriteLine("  TRIANGLE_SAMPLE samples points from a triangle.");

            typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

            Console.WriteLine("");
            Console.WriteLine("  Sample points (X,Y) and (XSI1,XSI2,XSI3) coordinates:");
            Console.WriteLine("");

            for (test = 0; test < 10; test++)
            {
                typeMethods.triangle_sample(t, 1, ref seed, ref p);
                xsi = typeMethods.triangle_xy_to_xsi(t, p);
                Console.WriteLine("  " + p[0].ToString().PadLeft(8)
                                       + "  " + p[1].ToString().PadLeft(8)
                                       + "  "
                                       + "  " + xsi[0].ToString().PadLeft(8)
                                       + "  " + xsi[1].ToString().PadLeft(8)
                                       + "  " + xsi[2].ToString().PadLeft(8) + "");
            }
        }

        static void triangle_xsi_to_xy_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_XSI_TO_XY_TEST tests TRIANGLE_XSI_TO_XY.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] p = null;
            double[] p2;
            int seed = 123456789;
            double[] t =
            {
                4.0, 2.0,
                1.0, 5.0,
                -2.0, 2.0
            };
            int test;
            double[] xsi;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_XSI_TO_XY_TEST");
            Console.WriteLine("  TRIANGLE_XSI_TO_XY converts XSI to XY coordinates.");
            Console.WriteLine("");
            Console.WriteLine("  We verify that (X,Y) -> (XSI1,XSI2,XSI3) -> (X,Y)");
            Console.WriteLine("  works properly.");

            typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

            Console.WriteLine("");
            Console.WriteLine("  Sample points:");
            Console.WriteLine("");

            for (test = 0; test < 10; test++)
            {

                if (test == 0)
                {
                    p = typeMethods.triangle_centroid_2d(t);
                }
                else if (test == 1)
                {
                    p = new double[2];
                    p[0] = 3.0;
                    p[1] = 0.0;
                }
                else
                {
                    typeMethods.triangle_sample(t, 1, ref seed, ref p);
                }

                xsi = typeMethods.triangle_xy_to_xsi(t, p);

                p2 = typeMethods.triangle_xsi_to_xy(t, xsi);

                Console.WriteLine("  " + p[0].ToString().PadLeft(8)
                                       + "  " + p[1].ToString().PadLeft(8)
                                       + "  "
                                       + "  " + xsi[0].ToString().PadLeft(8)
                                       + "  " + xsi[1].ToString().PadLeft(8)
                                       + "  " + xsi[2].ToString().PadLeft(8)
                                       + "  "
                                       + "  " + p2[0].ToString().PadLeft(8)
                                       + "  " + p2[1].ToString().PadLeft(8) + "");

            }
        }

        static void triangle_xy_to_xsi_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_XY_TO_XSI_TEST tests TRIANGLE_XY_TO_XSI.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] p = null;
            double[] p2;
            int seed = 123456789;
            double[] t =
            {
                4.0, 2.0,
                1.0, 5.0,
                -2.0, 2.0
            };
            int test;
            double[] xsi;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_XY_TO_XSI_TEST");
            Console.WriteLine("  TRIANGLE_XY_TO_XSI converts XY to XSI coordinates.");
            Console.WriteLine("");
            Console.WriteLine("  We verify that (X,Y) -> (XSI1,XSI2,XSI3) -> (X,Y)");
            Console.WriteLine("  works properly.");

            typeMethods.r8mat_transpose_print(2, 3, t, "  Triangle vertices:");

            Console.WriteLine("");
            Console.WriteLine("  Sample points:");
            Console.WriteLine("");

            for (test = 0; test < 10; test++)
            {

                if (test == 0)
                {
                    p = typeMethods.triangle_centroid_2d(t);
                }
                else if (test == 1)
                {
                    p = new double[2];
                    p[0] = 3.0;
                    p[1] = 0.0;
                }
                else
                {
                    typeMethods.triangle_sample(t, 1, ref seed, ref p);
                }

                xsi = typeMethods.triangle_xy_to_xsi(t, p);

                p2 = typeMethods.triangle_xsi_to_xy(t, xsi);

                Console.WriteLine("  " + p[0].ToString().PadLeft(8)
                                       + "  " + p[1].ToString().PadLeft(8)
                                       + "  "
                                       + "  " + xsi[0].ToString().PadLeft(8)
                                       + "  " + xsi[1].ToString().PadLeft(8)
                                       + "  " + xsi[2].ToString().PadLeft(8)
                                       + "  "
                                       + "  " + p2[0].ToString().PadLeft(8)
                                       + "  " + p2[1].ToString().PadLeft(8) + "");
            }
        }
    }
}