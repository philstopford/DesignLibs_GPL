using System;
using Burkardt.Types;
using Burkardt.Voxels;

namespace GeometryTest
{
    public static class VoxelTest
    {
        public static void voxels_dist_l1_nd_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VOXELS_DIST_L1_ND_TEST tests VOXELS_DIST_L1_ND.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;

            int dist;
            int[] p1 = {1, 1, 5};
            int[] p2 = {9, 4, 4};

            Console.WriteLine("");
            Console.WriteLine("VOXELS_DIST_L1_ND_TEST");
            Console.WriteLine("  VOXELS_DIST_L1_ND prints the voxels on a line in 3D.");

            Console.WriteLine("");
            Console.WriteLine("  P1:");
            Console.WriteLine("  " + p1[0].ToString().PadLeft(6)
                                   + "  " + p1[1].ToString().PadLeft(6)
                                   + "  " + p1[2].ToString().PadLeft(6) + "");


            Console.WriteLine("");
            Console.WriteLine("  P2:");
            Console.WriteLine("  " + p2[0].ToString().PadLeft(6)
                                   + "  " + p2[1].ToString().PadLeft(6)
                                   + "  " + p2[2].ToString().PadLeft(6) + "");

            dist = Geometry.voxels_dist_l1_nd(DIM_NUM, p1, p2);

            Console.WriteLine("");
            Console.WriteLine("  L1 distance = " + dist + "");

        }

        public static void voxels_line_3d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VOXELS_LINE_3D_TEST tests VOXELS_LINE_3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;

            int n;
            int[] v;
            int[] p1 = {1, 1, 5};
            int[] p2 = {9, 4, 4};

            Console.WriteLine("");
            Console.WriteLine("VOXELS_LINE_3D_TEST");
            Console.WriteLine("  VOXELS_LINE_3D computes the voxels on a line in 3D");
            Console.WriteLine("  starting at the first voxel, and heading towards");
            Console.WriteLine("  the second one.");

            Console.WriteLine("");
            Console.WriteLine("  Starting voxel:");
            Console.WriteLine("  " + p1[0].ToString().PadLeft(6)
                                   + "  " + p1[1].ToString().PadLeft(6)
                                   + "  " + p1[2].ToString().PadLeft(6) + "");

            Console.WriteLine("");
            Console.WriteLine("  Heading voxel:");
            Console.WriteLine("  " + p2[0].ToString().PadLeft(6)
                                   + "  " + p2[1].ToString().PadLeft(6)
                                   + "  " + p2[2].ToString().PadLeft(6) + "");

            n = Geometry.voxels_dist_l1_nd(DIM_NUM, p1, p2) + 1;

            Console.WriteLine("");
            Console.WriteLine("  Number of voxels we will compute is " + n + "");

            v = new int[3 * n];

            Geometry.voxels_line_3d(p1, p2, n, ref v);

            typeMethods.i4mat_transpose_print(3, n, v, "  The voxels:");

        }

        public static void voxels_region_3d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VOXELS_REGION_3D_TEST tests VOXELS_REGION_3D.
            //
            //  Discussion:
            //
            //    The test region is 8 by 9 by 1 voxels:
            //
            //    123456789
            //  1 .........
            //  2 ...11.1..
            //  3 ..11111..
            //  4 ...11.1..
            //  5 ......1..
            //  6 .11..11..
            //  7 ..1......
            //  8 .......1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 February 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;
            int LIST_MAX = 100;
            int NX = 8;
            int NY = 9;
            int NZ = 1;

            int i;
            int[] ishow = new int[NX * NY * NZ];
            int j;
            int k;
            int l;
            int[] list = new int[LIST_MAX];
            int list_num = 0;
            int nelements;
            int region;
            int region_num = 0;

            Console.WriteLine("");
            Console.WriteLine("VOXELS_REGION_3D_TEST");
            Console.WriteLine("  VOXELS_REGION_3D groups voxels into regions.");

            typeMethods.i4vec_zero(NX * NY * NZ, ref ishow);

            ishow[1 + 3 * NX + 0 * NX * NY] = 1;
            ishow[1 + 4 * NX + 0 * NX * NY] = 1;
            ishow[1 + 6 * NX + 0 * NX * NY] = 1;

            ishow[2 + 2 * NX + 0 * NX * NY] = 1;
            ishow[2 + 3 * NX + 0 * NX * NY] = 1;
            ishow[2 + 4 * NX + 0 * NX * NY] = 1;
            ishow[2 + 5 * NX + 0 * NX * NY] = 1;
            ishow[2 + 6 * NX + 0 * NX * NY] = 1;

            ishow[3 + 3 * NX + 0 * NX * NY] = 1;
            ishow[3 + 4 * NX + 0 * NX * NY] = 1;
            ishow[3 + 6 * NX + 0 * NX * NY] = 1;

            ishow[4 + 6 * NX + 0 * NX * NY] = 1;

            ishow[5 + 1 * NX + 0 * NX * NY] = 1;
            ishow[5 + 2 * NX + 0 * NX * NY] = 1;
            ishow[5 + 5 * NX + 0 * NX * NY] = 1;
            ishow[5 + 6 * NX + 0 * NX * NY] = 1;

            ishow[6 + 2 * NX + 0 * NX * NY] = 1;

            ishow[7 + 7 * NX + 0 * NX * NY] = 1;

            Geometry.voxels_region_3d(LIST_MAX, NX, NY, NZ, ref ishow, ref list_num, ref list,
                ref region_num);

            Console.WriteLine("");
            Console.WriteLine("  Number of regions found = " + region_num + "");
            Console.WriteLine("");
            Console.WriteLine("  The nonzero ISHOW array elements are:");
            Console.WriteLine("");

            for (i = 0; i < NX; i++)
            {
                for (j = 0; j < NY; j++)
                {
                    for (k = 0; k < NZ; k++)
                    {
                        l = ishow[i + j * NX + k * NX * NY];
                        if (l != 0)
                        {
                            Console.WriteLine("  " + (i + 1).ToString().PadLeft(6)
                                                   + "  " + (j + 1).ToString().PadLeft(6)
                                                   + "  " + (k + 1).ToString().PadLeft(6)
                                                   + "  " + l.ToString().PadLeft(6) + "");
                        }
                    }
                }
            }

            if (LIST_MAX < list_num)
            {
                Console.WriteLine("");
                Console.WriteLine("  The stack-based list of regions is unusable.");
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("  The stack-based list of regions is:");
                Console.WriteLine("");

                region = region_num;

                while (0 < list_num)
                {
                    nelements = list[list_num - 1];
                    list_num = list_num - 1;

                    Console.WriteLine("");
                    Console.WriteLine("  Region " + region
                                                  + " includes " + nelements + " voxels:");
                    Console.WriteLine("");

                    for (l = 1; l <= nelements; l++)
                    {
                        k = list[list_num - 1];
                        list_num = list_num - 1;
                        j = list[list_num - 1];
                        list_num = list_num - 1;
                        i = list[list_num - 1];
                        list_num = list_num - 1;
                        Console.WriteLine("  " + i.ToString().PadLeft(6)
                                               + "  " + j.ToString().PadLeft(6)
                                               + "  " + k.ToString().PadLeft(6) + "");
                    }

                    region = region - 1;
                }
            }

        }

        public static void voxels_step_3d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VOXELS_STEP_3D_TEST tests VOXELS_STEP_3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;

            int i;
            int inc;
            int jnc;
            int knc;
            int[] v1 = {1, 1, 5};
            int[] v2 = new int[DIM_NUM];
            int[] v3 = new int[DIM_NUM];

            Console.WriteLine("");
            Console.WriteLine("VOXELS_STEP_3D_TEST");
            Console.WriteLine("  VOXELS_STEP_3D steps along a line from");
            Console.WriteLine("  one voxel to another.");

            typeMethods.i4vec_copy(DIM_NUM, v1, ref v2);

            inc = 7;
            jnc = 3;
            knc = -1;

            Console.WriteLine("");
            Console.WriteLine("  " + 0.ToString().PadLeft(4)
                                   + "  " + v2[0].ToString().PadLeft(6)
                                   + "  " + v2[1].ToString().PadLeft(6)
                                   + "  " + v2[2].ToString().PadLeft(6) + "");

            for (i = 1; i <= 10; i++)
            {
                Geometry.voxels_step_3d(v1, v2, inc, jnc, knc, ref v3);
                Console.WriteLine("  " + i.ToString().PadLeft(4)
                                       + "  " + v3[0].ToString().PadLeft(6)
                                       + "  " + v3[1].ToString().PadLeft(6)
                                       + "  " + v3[2].ToString().PadLeft(6) + "");
                typeMethods.i4vec_copy(DIM_NUM, v3, ref v2);
            }

            Console.WriteLine("");
            Console.WriteLine("  Now, as a check, reverse direction and return.");
            Console.WriteLine("");

            typeMethods.i4vec_copy(DIM_NUM, v2, ref v1);

            inc = -inc;
            jnc = -jnc;
            knc = -knc;

            typeMethods.i4vec_copy(DIM_NUM, v1, ref v2);

            Console.WriteLine("  " + 0.ToString().PadLeft(4)
                                   + "  " + v2[0].ToString().PadLeft(6)
                                   + "  " + v2[1].ToString().PadLeft(6)
                                   + "  " + v2[2].ToString().PadLeft(6) + "");

            for (i = 1; i <= 10; i++)
            {
                Geometry.voxels_step_3d(v1, v2, inc, jnc, knc, ref v3);
                Console.WriteLine("  " + i.ToString().PadLeft(4)
                                       + "  " + v3[0].ToString().PadLeft(6)
                                       + "  " + v3[1].ToString().PadLeft(6)
                                       + "  " + v3[2].ToString().PadLeft(6) + "");
                typeMethods.i4vec_copy(DIM_NUM, v3, ref v2);
            }

        }

    }
}