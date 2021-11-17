using System;
using System.Linq;
using Burkardt.Plane;
using Burkardt.Types;
using Burkardt.Uniform;

namespace GeometryTest;

public static class PlaneTest
{
      public static void plane_exp_normal_3d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_EXP_NORMAL_3D_TEST tests PLANE_EXP_NORMAL_3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
      {
            int DIM_NUM = 3;

            double[] p1 = {-10.56, -10.56, 78.09};
            double[] p2 = {44.66, -65.77, 0.0};
            double[] p3 = {44.66, 44.66, 0.0};
            double[] pn = new double[DIM_NUM];

            Console.WriteLine("");
            Console.WriteLine("PLANE_EXP_NORMAL_3D_TEST");
            Console.WriteLine("  PLANE_EXP_NORMAL_3D finds the normal to a plane.");
            Console.WriteLine("");
            Console.WriteLine("  Coordinates of 3 points:");
            Console.WriteLine("");
            Console.WriteLine("  P1 = " + p1[0].ToString().PadLeft(12)
                                        + "  " + p1[1].ToString().PadLeft(12)
                                        + "  " + p1[2].ToString().PadLeft(12) + "");
            Console.WriteLine("  P2 = " + p2[0].ToString().PadLeft(12)
                                        + "  " + p2[1].ToString().PadLeft(12)
                                        + "  " + p2[2].ToString().PadLeft(12) + "");
            Console.WriteLine("  P3 = " + p3[0].ToString().PadLeft(12)
                                        + "  " + p3[1].ToString().PadLeft(12)
                                        + "  " + p3[2].ToString().PadLeft(12) + "");

            Geometry.plane_exp_normal_3d(p1, p2, p3, ref pn);

            Console.WriteLine("");
            Console.WriteLine("  Unit normal vector:");
            Console.WriteLine("");
            Console.WriteLine("  PN = " + pn[0].ToString().PadLeft(12)
                                        + "  " + pn[1].ToString().PadLeft(12)
                                        + "  " + pn[2].ToString().PadLeft(12) + "");

      }

      public static void test051()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST051 tests PLANE_EXP2IMP_3D.
            //
            //  Modified:
            //
            //    06 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
      {
            int DIM_NUM = 3;

            double a = 0;
            double b = 0;
            double c = 0;
            double d = 0;
            double[] p1 = new double[DIM_NUM];
            double[] p2 = new double[DIM_NUM];
            double[] p3 = new double[DIM_NUM];

            Console.WriteLine("");
            Console.WriteLine("TEST051");
            Console.WriteLine("  PLANE_EXP2IMP_3D puts a plane defined by");
            Console.WriteLine("  3 points into A*X+B*Y+C*Z+D = 0 form.");
            Console.WriteLine("");

            p1[0] = -1.0;
            p1[1] = 0.0;
            p1[2] = -1.0;

            p2[0] = -4.0;
            p2[1] = 0.0;
            p2[2] = 0.0;

            p3[0] = -20.0;
            p3[1] = 2.0;
            p3[2] = 4.0;

            Geometry.plane_exp2imp_3d(p1, p2, p3, ref a, ref b, ref c, ref d);

            Console.WriteLine("");
            Console.WriteLine("  Input:");
            Console.WriteLine("");
            Console.WriteLine("  P1 = " + p1[0].ToString().PadLeft(12)
                                        + "  " + p1[1].ToString().PadLeft(12)
                                        + "  " + p1[2].ToString().PadLeft(12) + "");
            Console.WriteLine("  P2 = " + p2[0].ToString().PadLeft(12)
                                        + "  " + p2[1].ToString().PadLeft(12)
                                        + "  " + p2[2].ToString().PadLeft(12) + "");
            Console.WriteLine("  P3 = " + p3[0].ToString().PadLeft(12)
                                        + "  " + p3[1].ToString().PadLeft(12)
                                        + "  " + p3[2].ToString().PadLeft(12) + "");
            Console.WriteLine("");
            Console.WriteLine("  Output:");
            Console.WriteLine("");
            Console.WriteLine("  (A,B,C,D)= " + a + "  " + b + "  " + c + "  " + d + "");
            Console.WriteLine("  Correct answer is a multiple of 1, 2, 3, 4.");

            p1[0] = -16.0;
            p1[1] = 2.0;
            p1[2] = 4.0;

            p2[0] = 0.0;
            p2[1] = 0.0;
            p2[2] = 0.0;

            p3[0] = 4.0;
            p3[1] = -2.0;
            p3[2] = 0.0;

            Geometry.plane_exp2imp_3d(p1, p2, p3, ref a, ref b, ref c, ref d);

            Console.WriteLine("");
            Console.WriteLine("  Input:");
            Console.WriteLine("");
            Console.WriteLine("  P1 = " + p1[0].ToString().PadLeft(12)
                                        + "  " + p1[1].ToString().PadLeft(12)
                                        + "  " + p1[2].ToString().PadLeft(12) + "");
            Console.WriteLine("  P2 = " + p2[0].ToString().PadLeft(12)
                                        + "  " + p2[1].ToString().PadLeft(12)
                                        + "  " + p2[2].ToString().PadLeft(12) + "");
            Console.WriteLine("  P3 = " + p3[0].ToString().PadLeft(12)
                                        + "  " + p3[1].ToString().PadLeft(12)
                                        + "  " + p3[2].ToString().PadLeft(12) + "");
            Console.WriteLine("");
            Console.WriteLine("  Output:");
            Console.WriteLine("");
            Console.WriteLine("  (A,B,C,D)= " + a + "  " + b + "  " + c + "  " + d + "");
            Console.WriteLine("  Correct answer is a multiple of 1, 2, 3, 0.");

      }

      public static void test052()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST052 tests PLANE_EXP2NORMAL_3D.
            //
            //  Modified:
            //
            //    06 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
      {
            int DIM_NUM = 3;

            int i;
            double[] p1 = new double[DIM_NUM];
            double[] p2 = new double[DIM_NUM];
            double[] p3 = new double[DIM_NUM];
            double[] pn = new double[DIM_NUM];
            double[] pp = new double[DIM_NUM];

            Console.WriteLine("");
            Console.WriteLine("TEST052");
            Console.WriteLine("  PLANE_EXP2NORMAL_3D puts a plane defined by");
            Console.WriteLine("  3 points into point, normal form");
            Console.WriteLine("");

            p1[0] = -1.0;
            p1[1] = 0.0;
            p1[2] = -1.0;

            p2[0] = -4.0;
            p2[1] = 0.0;
            p2[2] = 0.0;

            p3[0] = -20.0;
            p3[1] = 2.0;
            p3[2] = 4.0;

            Geometry.plane_exp2normal_3d(p1, p2, p3, ref pp, ref pn);

            Console.WriteLine("");
            Console.WriteLine("  Input:");
            Console.WriteLine("");
            Console.WriteLine("  P1 = " + p1[0].ToString().PadLeft(12)
                                        + "  " + p1[1].ToString().PadLeft(12)
                                        + "  " + p1[2].ToString().PadLeft(12) + "");
            Console.WriteLine("  P2 = " + p2[0].ToString().PadLeft(12)
                                        + "  " + p2[1].ToString().PadLeft(12)
                                        + "  " + p2[2].ToString().PadLeft(12) + "");
            Console.WriteLine("  P3 = " + p3[0].ToString().PadLeft(12)
                                        + "  " + p3[1].ToString().PadLeft(12)
                                        + "  " + p3[2].ToString().PadLeft(12) + "");
            Console.WriteLine("");
            Console.WriteLine("  Output:");
            Console.WriteLine("");
            Console.WriteLine("  PP = " + pp[0] + "  " + pp[1] + "  " + pp[2] + "");
            Console.WriteLine("  PN = " + pn[0] + "  " + pn[1] + "  " + pn[2] + "");

            p1[0] = -16.0;
            p1[1] = 2.0;
            p1[2] = 4.0;

            p2[0] = 0.0;
            p2[1] = 0.0;
            p2[2] = 0.0;

            p3[0] = 4.0;
            p3[1] = -2.0;
            p3[2] = 0.0;

            string cout = "  P1: ";
            for (i = 0; i < DIM_NUM; i++)
            {
                  cout += "  " + p1[i].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
            cout = "  P2: ";
            for (i = 0; i < DIM_NUM; i++)
            {
                  cout += "  " + p2[i].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
            cout = "  P3: ";
            for (i = 0; i < DIM_NUM; i++)
            {
                  cout += "  " + p3[i].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);

            Geometry.plane_exp2normal_3d(p1, p2, p3, ref pp, ref pn);

            Console.WriteLine("");
            Console.WriteLine("  Output:");
            Console.WriteLine("");
            Console.WriteLine("  PP = " + pp[0] + "  " + pp[1] + "  " + pp[2] + "");
            Console.WriteLine("  PN = " + pn[0] + "  " + pn[1] + "  " + pn[2] + "");

      }

      public static void test053()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST053 tests PLANE_EXP_PROJECT_3D.
            //
            //  Discussion:
            //
            //    1: Projection   is ( 0, 0.5, 0.5 ), IVIS is 3.
            //    2: Projection   is ( 4, 5, -8 ), IVIS is 2.
            //    3: Projection   is ( 0.33, 0.33, 0.33), IVIS is 1.
            //    4: "Projection" is ( 0, 0, 0 ), IVIS is 0.
            //    5: Projection   is ( 1, 0, 0 ), IVIS is -1.
            //
            //  Modified:
            //
            //    18 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
      {
            int DIM_NUM = 3;
            int N = 5;

            int i;
            int[] ivis = new int[N];
            int j;
            double[] p1 = {1.0, 0.0, 0.0};
            double[] p2 = {0.0, 1.0, 0.0};
            double[] p3 = {0.0, 0.0, 1.0};
            double[] pf = {0.0, 0.0, 0.0};
            double[] po =
            {
                  0.00, 2.00, 2.00,
                  4.00, 5.00, -8.00,
                  0.25, 0.25, 0.25,
                  5.00, -2.00, -3.00,
                  -2.00, 0.00, 0.00
            };
            double[] pp = new double[DIM_NUM * N];

            Console.WriteLine("");
            Console.WriteLine("TEST053");
            Console.WriteLine("  PLANE_EXP_PROJECT_3D projects a point through");
            Console.WriteLine("  a focus point into a plane.");
            Console.WriteLine("");
            Console.WriteLine("      PO      PP      IVIS");
            Console.WriteLine("");

            Geometry.plane_exp_project_3d(p1, p2, p3, pf, N, po, ref pp, ref ivis);

            string cout = "";
            for (j = 0; j < N; j++)
            {
                  for (i = 0; i < DIM_NUM; i++)
                  {
                        cout += "  " + po[i + j * DIM_NUM].ToString().PadLeft(10);
                  }

                  for (i = 0; i < DIM_NUM; i++)
                  {
                        cout += "  " + pp[i + j * DIM_NUM].ToString().PadLeft(10);
                  }

                  cout += "  " + ivis[j].ToString().PadLeft(4);
                  Console.WriteLine(cout);
            }

      }

      public static void test054()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST054 tests PLANE_IMP2EXP_3D.
            //
            //  Modified:
            //
            //    06 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
      {
            int DIM_NUM = 3;

            double a;
            double b;
            double c;
            double d;
            int i;
            double[] p1 = new double[DIM_NUM];
            double[] p2 = new double[DIM_NUM];
            double[] p3 = new double[DIM_NUM];

            Console.WriteLine("");
            Console.WriteLine("TEST054");
            Console.WriteLine("  PLANE_IMP2EXP_3D converts a plane in implicit");
            Console.WriteLine("  (A,B,C,D) form to explicit form.");

            a = 1.0;
            b = -2.0;
            c = -3.0;
            d = 6.0;

            Console.WriteLine("");
            Console.WriteLine("  A = " + a
                                       + "  B = " + b
                                       + "  C = " + c
                                       + "  D = " + d + "");

            Geometry.plane_imp2exp_3d(a, b, c, d, ref p1, ref p2, ref p3);

            string cout = "  P1: ";
            for (i = 0; i < DIM_NUM; i++)
            {
                  cout += "  " + p1[i].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
            cout = "  P2: ";
            for (i = 0; i < DIM_NUM; i++)
            {
                  cout += "  " + p2[i].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
            cout = "  P3: ";
            for (i = 0; i < DIM_NUM; i++)
            {
                  cout += "  " + p3[i].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);

      }

      public static void plane_imp2normal_3d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_IMP2NORMAL_3D_TEST tests PLANE_IMP2NORMAL_3D.
            //
            //  Modified:
            //
            //    06 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
      {
            int DIM_NUM = 3;

            double a;
            double b;
            double c;
            double d;
            double[] pn = new double[DIM_NUM];
            double[] pp = new double[DIM_NUM];

            Console.WriteLine("");
            Console.WriteLine("PLANE_IMP2NORMAL_3D_TEST");
            Console.WriteLine("  PLANE_IMP2NORMAL_3D converts a plane in implicit");
            Console.WriteLine("  (A,B,C,D) form to point, normal form.");

            a = 1.0;
            b = -2.0;
            c = -3.0;
            d = 6.0;

            Console.WriteLine("");
            Console.WriteLine("  A = " + a
                                       + "  B = " + b
                                       + "  C = " + c
                                       + "  D = " + d + "");

            Geometry.plane_imp2normal_3d(a, b, c, d, ref pp, ref pn);

            Console.WriteLine("");
            Console.WriteLine("  PP = " + pp[0] + "  " + pp[1] + "  " + pp[2] + "");
            Console.WriteLine("  PN = " + pn[0] + "  " + pn[1] + "  " + pn[2] + "");

      }

      public static void plane_imp_line_par_int_3d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_IMP_LINE_PAR_INT_3D_TEST tests PLANE_IMP_LINE_PAR_INT_3D.
            //
            //  Modified:
            //
            //    18 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
      {
            int DIM_NUM = 3;

            double a;
            double b;
            double c;
            double d;
            double f;
            double g;
            double h;
            int i;
            bool intersect;
            double[] p = new double[DIM_NUM];
            double x0;
            double y0;
            double z0;

            Console.WriteLine("");
            Console.WriteLine("PLANE_IMP_LINE_PAR_INT_3D_TEST");
            Console.WriteLine("  PLANE_IMP_LINE_PAR_INT_3D finds the");
            Console.WriteLine("  intersection of an implicit plane and");
            Console.WriteLine("  a parametric line, in 3D.");

            a = 1.0;
            b = -2.0;
            c = -3.0;
            d = 6.0;

            f = 2.0;
            g = 1.0;
            h = 5.0;
            x0 = 3.0;
            y0 = 0.0;
            z0 = -7.0;
            string cout = "";

            intersect = Geometry.plane_imp_line_par_int_3d(a, b, c, d, x0, y0, z0, f, g, h, ref p);

            switch (intersect)
            {
                  case true:
                  {
                        Console.WriteLine("");
                        Console.WriteLine("  The plane and line intersect at");
                        for (i = 0; i < DIM_NUM; i++)
                        {
                              cout += "  " + p[i].ToString().PadLeft(12);
                        }

                        Console.WriteLine(cout);
                        break;
                  }
                  default:
                        Console.WriteLine("");
                        Console.WriteLine("  The plane and the line do not intersect.");
                        break;
            }

            Console.WriteLine("");
            Console.WriteLine("  Expected answer:");
            Console.WriteLine("    The plane and line intersect at ");
            Console.WriteLine("    7, 2, 3.");

      }

      public static void test057()

            //****************************************************************************80*
            //
            //  Purpose:
            //
            //    TEST057 tests PLANE_IMP_SEGMENT_NEAR_3D.
            //
            //  Modified:
            //
            //    06 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
      {
            int DIM_NUM = 3;
            int TEST_NUM = 2;

            double a;
            double b;
            double c;
            double d;
            double dist = 0;
            double[] p1 = new double[DIM_NUM];
            double[] p2 = new double[DIM_NUM];
            double[] pls = new double[DIM_NUM];
            double[] pp = new double[DIM_NUM];
            int test;

            Console.WriteLine("");
            Console.WriteLine("TEST057");
            Console.WriteLine("  PLANE_IMP_SEGMENT_NEAR_3D finds the point");
            Console.WriteLine("  on a line segment nearest a plane.");

            for (test = 0; test < TEST_NUM; test++)
            {
                  p1[0] = 3.0;
                  p1[1] = 0.0;
                  p1[2] = -7.0;

                  switch (test)
                  {
                        case 0:
                              p2[0] = 9.0;
                              p2[1] = 3.0;
                              p2[2] = 8.0;
                              break;
                        case 1:
                              p2[0] = 5.0;
                              p2[1] = 1.0;
                              p2[2] = -2.0;
                              break;
                  }

                  a = 1.0;
                  b = -2.0;
                  c = -3.0;
                  d = 6.0;

                  Geometry.plane_imp_segment_near_3d(p1, p2, a, b, c, d, ref dist, ref pp, ref pls);

                  Console.WriteLine("");
                  Console.WriteLine("  The distance between the plane and the");
                  Console.WriteLine("  line segment is " + dist + "");
                  Console.WriteLine("");
                  Console.WriteLine("  A nearest point on the line segment is");
                  Console.WriteLine("  PLS = " + pls[0] + "  " + pls[1] + "  " + pls[2] + "");
                  Console.WriteLine("");
                  Console.WriteLine("  A nearest point on the plane is");
                  Console.WriteLine("  PP = " + pp[0] + "  " + pp[1] + "  " + pp[2] + "");
            }
      }

      public static void test058()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST058 tests PLANE_IMP_POINT_DIST_3D, PLANE_IMP_POINT_DIST_SIGNED_3D;
            //
            //  Modified:
            //
            //    06 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
      {
            int TEST_NUM = 4;

            double a;
            double b;
            double c;
            double d;
            double dist;
            double dist_signed;
            double[] p = new double[3];
            int test;
            double[] xtest = new double[TEST_NUM];
            double[] ytest = new double[TEST_NUM];
            double[] ztest = new double[TEST_NUM];

            a = 0.0;
            b = 0.0;
            c = 1.0;
            d = -10.0;

            xtest[0] = -12.0;
            ytest[0] = 14.0;
            ztest[0] = 0.0;

            xtest[1] = 7.0;
            ytest[1] = 8.0;
            ztest[1] = 9.0;

            xtest[2] = 1.0;
            ytest[2] = 2.0;
            ztest[2] = 10.0;

            xtest[3] = 0.0;
            ytest[3] = 0.0;
            ztest[3] = 12.0;

            Console.WriteLine("");
            Console.WriteLine("TEST058");
            Console.WriteLine("  PLANE_IMP_POINT_DIST_3D computes the distance");
            Console.WriteLine("  between an implicit plane and a point in 3D;");
            Console.WriteLine("  PLANE_IMP_POINT_DIST_SIGNED 3D computes the");
            Console.WriteLine("  signed distance between an implicit plane");
            Console.WriteLine("  and a point in 3D.");
            Console.WriteLine("");
            Console.WriteLine("  For all tests, we use the implicit plane with");
            Console.WriteLine("  (A,B,C,D) = "
                              + a + "  " + b + "  " + c + "  " + d + "");
            Console.WriteLine("");
            Console.WriteLine("  P     DISTANCE   SIGNED_DISTANCE");
            Console.WriteLine("");

            for (test = 0; test < TEST_NUM; test++)
            {
                  p[0] = xtest[test];
                  p[1] = ytest[test];
                  p[2] = ztest[test];

                  dist = Geometry.plane_imp_point_dist_3d(a, b, c, d, p);

                  dist_signed = Geometry.plane_imp_point_dist_signed_3d(a, b, c, d, p);

                  Console.WriteLine("  " + p[0].ToString().PadLeft(10)
                                         + "  " + p[1].ToString().PadLeft(10)
                                         + "  " + p[2].ToString().PadLeft(10)
                                         + "  " + dist.ToString().PadLeft(10)
                                         + "  " + dist_signed.ToString().PadLeft(10) + "");
            }

      }

      public static void test059()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST059 tests PLANE_IMP_TRIANGLE_NEAR_3D.
            //
            //  Modified:
            //
            //    18 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
      {
            int DIM_NUM = 3;
            int TEST_NUM = 2;

            double a;
            double b;
            double c;
            double d;
            double dist = 0;
            int near_num;
            double[] pn = new double[DIM_NUM * 6];
            double[] t =
            {
                  3.0, 0.0, -7.0,
                  13.0, -4.0, -1.0,
                  0.0, 0.0, 0.0
            };
            int test;

            Console.WriteLine("");
            Console.WriteLine("TEST059");
            Console.WriteLine("  PLANE_IMP_TRIANGLE_NEAR_3D finds the nearest");
            Console.WriteLine("  points on an implicit plane and a triangle.");

            a = 1.0;
            b = -2.0;
            c = -3.0;
            d = 6.0;

            Console.WriteLine("");
            Console.WriteLine("  Implicit plane: A*X + B*Y + C*Z + D = 0.");
            Console.WriteLine("  A = " + a
                                       + "  B = " + b
                                       + "  C = " + c
                                       + "  D = " + d + "");

            for (test = 0; test < TEST_NUM; test++)
            {
                  switch (test)
                  {
                        case 1:
                              t[0 + 2 * DIM_NUM] = 5.0;
                              t[1 + 2 * DIM_NUM] = 1.0;
                              t[2 + 2 * DIM_NUM] = -2.0;
                              break;
                        case 2:
                              t[0 + 2 * DIM_NUM] = 9.0;
                              t[1 + 2 * DIM_NUM] = 3.0;
                              t[2 + 2 * DIM_NUM] = 8.0;
                              break;
                  }

                  typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

                  near_num = Geometry.plane_imp_triangle_near_3d(t, a, b, c, d, ref dist, ref pn);

                  Console.WriteLine("");
                  Console.WriteLine("  Triangle to plane distance is " + dist + "");

                  typeMethods.r8mat_transpose_print(DIM_NUM, near_num, pn, "  Nearest points:");
            }
      }

      public static void test060()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST060 tests PLANE_IMP_TRIANGLE_INT_3D.
            //
            //  Modified:
            //
            //    06 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
      {
            int DIM_NUM = 3;
            int TEST_NUM = 4;

            double a;
            double b;
            double c;
            double d;
            int int_num = 0;
            int j;
            double[] p = new double[DIM_NUM * 3];
            double[] t = new double[DIM_NUM * 3];
            int test;

            Console.WriteLine("");
            Console.WriteLine("TEST060");
            Console.WriteLine("  PLANE_IMP_TRIANGLE_INT_3D finds the");
            Console.WriteLine("  intersection points of an implicit plane");
            Console.WriteLine("  and a triangle.");

            a = 1.0;
            b = -2.0;
            c = -3.0;
            d = 6.0;

            Console.WriteLine("");
            Console.WriteLine("  The implicit plane: A*X + B*Y + C*Z + D = 0.");
            Console.WriteLine("  (A,B,C,D) = " + a + "  " + b + "  " + c + "  " + d + "");

            for (test = 0; test < TEST_NUM; test++)
            {
                  switch (test)
                  {
                        case 0:
                              t[0 + 0 * 3] = 3.0;
                              t[1 + 0 * 3] = 0.0;
                              t[2 + 0 * 3] = -7.0;
                              t[0 + 1 * 3] = 13.0;
                              t[1 + 1 * 3] = -4.0;
                              t[2 + 1 * 3] = -1.0;
                              t[0 + 2 * 3] = 5.0;
                              t[1 + 2 * 3] = 1.0;
                              t[2 + 2 * 3] = -2.0;
                              break;
                        case 1:
                              t[0 + 0 * 3] = 3.0;
                              t[1 + 0 * 3] = 0.0;
                              t[2 + 0 * 3] = -7.0;
                              t[0 + 1 * 3] = 13.0;
                              t[1 + 1 * 3] = -4.0;
                              t[2 + 1 * 3] = -1.0;
                              t[0 + 2 * 3] = 9.0;
                              t[1 + 2 * 3] = 3.0;
                              t[2 + 2 * 3] = 8.0;
                              break;
                        case 2:
                              t[0 + 0 * 3] = -6.0;
                              t[1 + 0 * 3] = 0.0;
                              t[2 + 0 * 3] = 0.0;
                              t[0 + 1 * 3] = 0.0;
                              t[1 + 1 * 3] = 3.0;
                              t[2 + 1 * 3] = 0.0;
                              t[0 + 2 * 3] = 0.0;
                              t[1 + 2 * 3] = 0.0;
                              t[2 + 2 * 3] = 2.0;
                              break;
                        case 3:
                              t[0 + 0 * 3] = -4.0;
                              t[1 + 0 * 3] = +1.0;
                              t[2 + 0 * 3] = 0.0;
                              t[0 + 1 * 3] = 0.0;
                              t[1 + 1 * 3] = 6.0;
                              t[2 + 1 * 3] = -2.0;
                              t[0 + 2 * 3] = 0.0;
                              t[1 + 2 * 3] = 0.0;
                              t[2 + 2 * 3] = 1.0;
                              break;
                  }

                  Console.WriteLine("");
                  Console.WriteLine("  Case " + test + "");

                  typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

                  Geometry.plane_imp_triangle_int_3d(a, b, c, d, t, ref int_num, ref p);

                  Console.WriteLine("");
                  Console.WriteLine("  Number of intersection points is " + int_num + "");
                  Console.WriteLine("");

                  switch (int_num)
                  {
                        case > 0:
                        {
                              for (j = 0; j < int_num; j++)
                              {
                                    Console.WriteLine("  " + (j + 1).ToString().PadLeft(4)
                                                           + "  " + p[0 + j * 3].ToString().PadLeft(10)
                                                           + "  " + p[1 + j * 3].ToString().PadLeft(10)
                                                           + "  " + p[2 + j * 3].ToString().PadLeft(10) + "");
                              }

                              break;
                        }
                  }

            }

      }

      public static void test061()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST061 tests PLANE_NORMAL_BASIS_3D.
            //
            //  Modified:
            //
            //    27 August 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
      {
            int DIM_NUM = 3;
            int TEST_NUM = 5;

            double[] b = new double[DIM_NUM * DIM_NUM];
            double[] pn;
            double[] pp;
            double[] pq = new double[DIM_NUM];
            double[] pr = new double[DIM_NUM];
            int seed = 123456789;
            int test;

            Console.WriteLine("");
            Console.WriteLine("TEST061");
            Console.WriteLine("  PLANE_NORMAL_BASIS_3D, given a plane in");
            Console.WriteLine("  point, normal form (P,N), finds two unit");
            Console.WriteLine("  vectors Q and R that lie in the plane");
            Console.WriteLine("  and are mutually orthogonal.");

            for (test = 0; test < TEST_NUM; test++)
            {
                  pp = UniformRNG.r8vec_uniform_01_new(DIM_NUM, ref seed);
                  pn = UniformRNG.r8vec_uniform_01_new(DIM_NUM, ref seed);

                  Geometry.plane_normal_basis_3d(pp, pn, ref pq, ref pr);

                  switch (test)
                  {
                        case 0:
                              Console.WriteLine("");
                              Console.WriteLine("  Data for first case:");
                              Console.WriteLine("");
                              Console.WriteLine("  PP = " + pp[0] + "  " + pp[1] + "  " + pp[2] + "");
                              Console.WriteLine("  PN = " + pn[0] + "  " + pn[1] + "  " + pn[2] + "");
                              Console.WriteLine("  PQ = " + pq[0] + "  " + pq[1] + "  " + pq[2] + "");
                              Console.WriteLine("  PR = " + pr[0] + "  " + pr[1] + "  " + pr[2] + "");
                              Console.WriteLine("");
                              break;
                  }

                  b[0 + 0 * 3] = typeMethods.r8vec_dot_product(DIM_NUM, pn, pn);
                  b[1 + 0 * 3] = typeMethods.r8vec_dot_product(DIM_NUM, pn, pq);
                  b[2 + 0 * 3] = typeMethods.r8vec_dot_product(DIM_NUM, pn, pr);
                  b[0 + 1 * 3] = typeMethods.r8vec_dot_product(DIM_NUM, pq, pn);
                  b[1 + 1 * 3] = typeMethods.r8vec_dot_product(DIM_NUM, pq, pq);
                  b[2 + 1 * 3] = typeMethods.r8vec_dot_product(DIM_NUM, pq, pr);
                  b[0 + 2 * 3] = typeMethods.r8vec_dot_product(DIM_NUM, pr, pn);
                  b[1 + 2 * 3] = typeMethods.r8vec_dot_product(DIM_NUM, pr, pq);
                  b[2 + 2 * 3] = typeMethods.r8vec_dot_product(DIM_NUM, pr, pr);

                  typeMethods.r8mat_print(DIM_NUM, DIM_NUM, b, "  Dot product matrix:");

            }
      }

      public static void test0615()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST0615 tests PLANE_NORMAL_LINE_EXP_INT_3D.
            //
            //  Modified:
            //
            //    18 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
      {
            int DIM_NUM = 3;

            int i;
            int ival;
            double[] normal = {1.0, -2.0, -3.0};
            double[] p1 = {3.0, 0.0, -7.0};
            double[] p2 = {5.0, 1.0, -2.0};
            double[] pint = new double[DIM_NUM];
            double[] pp = {-1.0, +1.0, +1.0};
            double temp;

            Console.WriteLine("");
            Console.WriteLine("TEST0615");
            Console.WriteLine("  PLANE_NORMAL_LINE_EXP_INT_3D finds the");
            Console.WriteLine("  intersection of a normal plane and");
            Console.WriteLine("  an explicit line, in 3D.");

            temp = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                  temp += Math.Pow(normal[i], 2);
            }

            temp = Math.Sqrt(temp);
            for (i = 0; i < DIM_NUM; i++)
            {
                  normal[i] /= temp;
            }

            typeMethods.r8vec_print(DIM_NUM, pp, "  Plane point PP:");
            typeMethods.r8vec_print(DIM_NUM, normal, "  Plane Normal:");

            typeMethods.r8vec_print(DIM_NUM, p1, "  Line point P1:");
            typeMethods.r8vec_print(DIM_NUM, p2, "  Line point P2:");

            ival = Geometry.plane_normal_line_exp_int_3d(pp, normal, p1, p2, ref pint);

            Console.WriteLine("");

            switch (ival)
            {
                  case 0:
                        Console.WriteLine("  The plane and line do not intersect.");
                        break;
                  case 1:
                        typeMethods.r8vec_print(DIM_NUM, pint, "  The unique intersection point:");
                        break;
                  case 2:
                        typeMethods.r8vec_print(DIM_NUM, pint, "  One of infinitely many intersections:");
                        break;
                  default:
                        Console.WriteLine("  The plane and the line do not intersect.");
                        break;
            }

            Console.WriteLine("");
            Console.WriteLine("  Expected answer:");
            Console.WriteLine("    The unique intersection point:");
            Console.WriteLine("    7, 2, 3.");

      }

      public static void test0616()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST0616 tests PLANE_NORMAL_QR_TO_XYZ and PLANE_NORMAL_XYZ_TO_QR.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 November 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
      {
            double dif;
            int i;
            int j;
            int m = 3;
            int n = 5;
            double[] normal;
            double[] pp;
            double[] pq = new double[3];
            double[] pr = new double[3];
            double[] qr1;
            double[] qr2;
            int seed;
            double t;
            double[] xyz;

            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TEST0616");
            Console.WriteLine("  For a normal plane, with point PP and NORMAL vector,");
            Console.WriteLine("  and in-plane basis vectors PQ and PR,");
            Console.WriteLine("  PLANE_NORMAL_QR_TO_XYZ converts QR to XYZ coordinates;");
            Console.WriteLine("  PLANE_NORMAL_XYZ_TO_QR converts XYZ to QR coordinates.");
            //
            //  Choose PP and NORMAL at random.
            //
            pp = UniformRNG.r8vec_uniform_01_new(m, ref seed);

            normal = UniformRNG.r8vec_uniform_01_new(m, ref seed);
            //
            //  Compute in-plane basis vectors PQ and PR.
            //
            Geometry.plane_normal_basis_3d(pp, normal, ref pq, ref pr);
            //
            //  Choose random Q, R coordinates.
            //
            qr1 = UniformRNG.r8mat_uniform_01_new(m - 1, n, ref seed);
            //
            //  Convert to XYZ.
            //
            xyz = Geometry.plane_normal_qr_to_xyz(pp, normal, pq, pr, n, qr1);
            //
            //  Convert XYZ to QR.
            //
            qr2 = Geometry.plane_normal_xyz_to_qr(pp, normal, pq, pr, n, xyz);

            dif = 0.0;
            for (j = 0; j < n; j++)
            {
                  t = 0.0;
                  for (i = 0; i < m - 1; i++)
                  {
                        t += Math.Pow(qr1[0 + j * 2] - qr2[0 + j * 2], 2);
                  }

                  t = Math.Sqrt(t);
                  dif = Math.Max(dif, t);
            }

            Console.WriteLine("");
            Console.WriteLine("  Maximum difference was " + dif + "");

      }

      public static void test0617()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST0617 tests PLANE_NORMAL_TETRAHEDRON_INTERSECT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 June 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
      {
            int i;
            int j;
            int k;
            int l;
            int int_num = 0;
            double[] normal = new double[3];
            double[] pint = new double[3 * 4];
            double[] pp = new double[3];
            double[] t =
            {
                  0.0, 0.0, 0.0,
                  1.0, 0.0, 0.0,
                  0.0, 1.0, 0.0,
                  0.0, 0.0, 1.0
            };

            string cout = "";

            Console.WriteLine("");
            Console.WriteLine("TEST0617");
            Console.WriteLine("  PLANE_NORMAL_TETRAHEDRON_INTERSECT determines");
            Console.WriteLine("  the intersection of a plane and tetrahedron.");

            for (k = 1; k <= 2; k++)
            {
                  switch (k)
                  {
                        case 1:
                              normal[0] = 0.0;
                              normal[1] = 0.0;
                              normal[2] = 1.0;
                              break;
                        default:
                              normal[0] = 1.0 / Math.Sqrt(2.0);
                              normal[1] = 1.0 / Math.Sqrt(2.0);
                              normal[2] = 0.0;
                              break;
                  }

                  Console.WriteLine("");
                  Console.WriteLine("  Plane normal vector number " + k + "");
                  Console.WriteLine("");
                  cout = "";
                  for (i = 0; i < 3; i++)
                  {
                        cout += "  " + normal[i];
                  }

                  Console.WriteLine(cout);

                  for (l = 0; l <= 6; l++)
                  {
                        for (i = 0; i < 3; i++)
                        {
                              pp[i] = normal[i] * l / 5.0;
                        }

                        Geometry.plane_normal_tetrahedron_intersect(pp, normal, t, ref int_num, ref pint);

                        Console.WriteLine("");
                        Console.WriteLine("  Point on plane:");
                        Console.WriteLine("");
                        cout = "";
                        for (i = 0; i < 3; i++)
                        {
                              cout += "  " + pp[i];
                        }

                        Console.WriteLine(cout);
                        Console.WriteLine("");
                        Console.WriteLine("  Number of intersection points = " + int_num + "");
                        Console.WriteLine("");
                        for (j = 0; j < int_num; j++)
                        {
                              cout = "  " + j.ToString().PadLeft(4);
                              for (i = 0; i < 3; i++)
                              {
                                    cout += "  " + pint[i + j * 3];
                              }

                              Console.WriteLine(cout);
                        }
                  }
            }
      }

      public static void test062()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST062 tests PLANE_NORMAL_TRIANGLE_INT_3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
      {
            int DIM_NUM = 3;
            int TEST_NUM = 4;

            int int_num;
            double[] normal = {1.0, -2.0, -3.0};
            double[] pint = new double[DIM_NUM * 3];
            double[] pp = {0.0, 0.0, 2.0};
            double[] t;
            double[] t_test =
            {
                  3.0, 0.0, -7.0,
                  13.0, -4.0, -1.0,
                  5.0, 1.0, -2.0,
                  3.0, 0.0, -7.0,
                  13.0, -4.0, -1.0,
                  9.0, 3.0, 8.0,
                  -6.0, 0.0, 0.0,
                  0.0, 3.0, 0.0,
                  0.0, 0.0, 2.0,
                  -4.0, 1.0, 0.0,
                  0.0, 6.0, -2.0,
                  0.0, 0.0, 1.0
            };
            int test;

            Console.WriteLine("");
            Console.WriteLine("TEST062");
            Console.WriteLine("  PLANE_NORMAL_TRIANGLE_INT_3D finds the");
            Console.WriteLine("  intersection points of a normal form plane");
            Console.WriteLine("  and a triangle.");

            typeMethods.r8vec_print(DIM_NUM, pp, "  The point PP:");
            typeMethods.r8vec_print(DIM_NUM, normal, "  The normal vector N:");

            for (test = 0; test < TEST_NUM; test++)
            {
                  t = t_test.Skip(+DIM_NUM * 3 * test).ToArray();

                  typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

                  int_num = Geometry.plane_normal_triangle_int_3d(pp, normal, t, pint);

                  Console.WriteLine("");
                  Console.WriteLine("  Number of intersection points is " + int_num + "");
                  Console.WriteLine("");

                  typeMethods.r8mat_transpose_print(DIM_NUM, int_num, pint, "  Intersection points:");
            }
      }

      public static void test063()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST063 tests PLANE_NORMAL2EXP_3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
      {
            int DIM_NUM = 3;

            int i;
            double[] normal = {-0.2672612, -0.5345225, -0.8017837};
            double[] p1 = new double[DIM_NUM];
            double[] p2 = new double[DIM_NUM];
            double[] p3 = new double[DIM_NUM];
            double[] pp = {-1.0, 0.0, -1.0};

            Console.WriteLine("");
            Console.WriteLine("TEST063");
            Console.WriteLine("  PLANE_NORMAL2EXP_3D puts a plane defined by");
            Console.WriteLine("  point, normal form into explicit form.");

            typeMethods.r8vec_print(DIM_NUM, pp, "  The point PP:");

            typeMethods.r8vec_print(DIM_NUM, normal, "  Normal vector:");

            Geometry.plane_normal2exp_3d(pp, normal, ref p1, ref p2, ref p3);

            string cout = "";
            cout = "  P1: ";
            for (i = 0; i < DIM_NUM; i++)
            {
                  cout += "  " + p1[i].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
            cout = "  P2: ";
            for (i = 0; i < DIM_NUM; i++)
            {
                  cout += "  " + p2[i].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
            cout = "  P3: ";
            for (i = 0; i < DIM_NUM; i++)
            {
                  cout += "  " + p3[i].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);

      }

      public static void test064()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST064 tests PLANE_NORMAL2IMP_3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
      {
            int DIM_NUM = 3;

            double a = 0;
            double b = 0;
            double c = 0;
            double d = 0;
            double[] pn = new double[DIM_NUM];
            double[] pp = new double[DIM_NUM];

            Console.WriteLine("");
            Console.WriteLine("TEST064");
            Console.WriteLine("  PLANE_NORMAL2IMP_3D puts a plane defined by");
            Console.WriteLine("  point, normal form into implicit ABCD form.");
            Console.WriteLine("");

            pp[0] = -1.0;
            pp[1] = 0.0;
            pp[2] = -1.0;

            pn[0] = -0.2672612;
            pn[1] = -0.5345225;
            pn[2] = -0.8017837;

            Console.WriteLine("");
            Console.WriteLine("  Input:");
            Console.WriteLine("");
            Console.WriteLine("  PP = " + pp[0] + "  " + pp[1] + "  " + pp[2] + "");
            Console.WriteLine("  PN = " + pn[0] + "  " + pn[1] + "  " + pn[2] + "");

            Geometry.plane_normal2imp_3d(pp, pn, ref a, ref b, ref c, ref d);

            Console.WriteLine("");
            Console.WriteLine("  Output:");
            Console.WriteLine("");
            Console.WriteLine("  (A,B,C,D) = " + a + "  " + b + "  " + c + "  " + d + "");

            pp[0] = -16.0;
            pp[1] = 2.0;
            pp[2] = 4.0;

            pn[0] = -0.2672612;
            pn[1] = -0.5345225;
            pn[2] = -0.8017837;

            Geometry.plane_normal2imp_3d(pp, pn, ref a, ref b, ref c, ref d);

            Console.WriteLine("");
            Console.WriteLine("  Input:");
            Console.WriteLine("");
            Console.WriteLine("  PP = " + pp[0] + "  " + pp[1] + "  " + pp[2] + "");
            Console.WriteLine("  PN = " + pn[0] + "  " + pn[1] + "  " + pn[2] + "");
            Console.WriteLine("");
            Console.WriteLine("  Output:");
            Console.WriteLine("");
            Console.WriteLine("  (A,B,C,D) = " + a + "  " + b + "  " + c + "  " + d + "");

      }

      public static void plane_exp_pro3_test ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_EXP_PRO3_TEST tests PLANE_EXP_PRO3.
            //
            //  Discussion:
            //
            //    Projection is ( -1, 1, 1 ).
            //    Projection is ( 4, 5, -8 ).
            //    Projection is ( 0.33, 0.33, 0.33).
            //    Projection is ( 5.33, -1.66, -2.66 ).
            //    Projection is ( -1, 1, 1 ).
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
            int TEST_NUM = 5;

            int i;
            double[] p1 = { 1.0, 0.0, 0.0 };
            double[] p2 = { 0.0, 1.0, 0.0 };
            double[] p3 = { 0.0, 0.0, 1.0 };
            double[] po = {
                  0.0,  2.0,  2.0,
                  4.0,  5.0, -8.0,
                  0.25, 0.25, 0.25,
                  5.0, -2.0, -3.0,
                  -2.0,  0.0,  0.0 };
            double[] pp = new double[DIM_NUM*TEST_NUM];
            int test;

            Console.WriteLine("");
            Console.WriteLine("PLANE_EXP_PRO3_TEST");
            Console.WriteLine("  PLANE_EXP_PRO3 projects an object point");
            Console.WriteLine("  orthographically into a plane.");
            Console.WriteLine("");
            Console.WriteLine("            PO                PP");
            Console.WriteLine("");

            Geometry.plane_exp_pro3 ( p1, p2, p3, TEST_NUM, po, ref pp );

            for ( test = 0; test < TEST_NUM; test++ )
            {
                  string cout = "";
                  for ( i = 0; i < DIM_NUM; i++ )
                  {
                        cout += "  " + po[i+test*DIM_NUM].ToString().PadLeft(10);
                  }
                  for ( i = 0; i < DIM_NUM; i++ )
                  {
                        cout += "  " + pp[i+test*DIM_NUM].ToString().PadLeft(10);
                  }
                  Console.WriteLine(cout);
            }
      }

}