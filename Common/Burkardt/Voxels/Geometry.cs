using System;
using Burkardt.Types;

namespace Burkardt.Voxels;

public static class Geometry
{
    public static int voxels_dist_l1_3d(int[] v1, int[] v2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VOXELS_DIST_L1_3D computes the L1 distance between voxels in 3D.
        //
        //  Discussion:
        //
        //    We can imagine that, in traveling from (X1,Y1,Z1) to (X2,Y2,Z2),
        //    we are allowed to increment or decrement just one coordinate at
        //    at time.  The minimum number of such changes required is the
        //    L1 distance.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int V1[3], the voxel that begins the line.
        //
        //    Input, int V2[3], the voxel that ends the line.
        //
        //    Output, int VOXELS_DIST_L1_3D, the L1 distance between the voxels.
        //
    {
        int value = Math.Abs(v1[0] - v2[0])
                    + Math.Abs(v1[1] - v2[1])
                    + Math.Abs(v1[2] - v2[2]);

        return value;
    }

    public static int voxels_dist_l1_nd(int dim_num, int[] v1, int[] v2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VOXELS_DIST_L1_ND computes the L1 distance between voxels in ND.
        //
        //  Discussion:
        //
        //    A voxel is generally a point in 3D space with integer coordinates.
        //    There's no reason to stick with 3D, so this routine will handle
        //    any dimension.
        //
        //    We can imagine that, in traveling from V1 to V2, we are allowed to
        //    increment or decrement just one coordinate at a time.  The minimum number
        //    of such changes required is the L1 distance.
        //
        //    More formally,
        //
        //      DIST_L1 ( V1, V2 ) = sum ( 1 <= I <= N ) | V1(I) - V2(I) |
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int V1[DIM_NUM], the voxel that begins the line.
        //
        //    Input, int V2[DIM_NUM], the voxel that ends the line.
        //
        //    Output, int VOXELS_DIST_L1_ND, the L1 distance between the voxels.
        //
    {
        int i;

        int value = 0;
        for (i = 0; i < dim_num; i++)
        {
            value += Math.Abs(v1[i] - v2[i]);
        }

        return value;
    }

    public static void voxels_line_3d(int[] v1, int[] v2, int n, ref int[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VOXELS_LINE_3D computes voxels along a line in 3D.
        //
        //  Discussion:
        //
        //    The line itself is defined by two voxels.  The line will begin
        //    at the first voxel, and move towards the second.  If the value of
        //    N is equal to the L1 distance between the two voxels, then the
        //    line will "almost" reach the second voxel.  Depending on the
        //    direction, 1, 2 or 3 more steps may be needed.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Daniel Cohen,
        //    Voxel Traversal along a 3D Line,
        //    Graphics Gems IV,
        //    edited by Paul Heckbert,
        //    AP Professional, 1994, T385.G6974.
        //
        //  Parameters:
        //
        //    Input, int V1[3], the voxel that begins the line.
        //
        //    Input, int V2[3], the voxel that ends the line.
        //
        //    Input, int N, the number of voxels to compute.
        //
        //    Output, int V[3*N], a sequence of voxels, whose
        //    first value is V1 and which proceeds towards V2.
        //
    {
        const int DIM_NUM = 3;

        int[] a = new int[DIM_NUM];
        int i;
        int j;
        int[] s = new int[DIM_NUM];

        switch (n)
        {
            case <= 0:
                return;
        }

        //
        //  Determine the number of voxels on the line.
        //
        for (i = 0; i < DIM_NUM; i++)
        {
            if (v2[i] < v1[i])
            {
                s[i] = -1;
            }
            else
            {
                s[i] = +1;
            }

            a[i] = Math.Abs(v2[i] - v1[i]);
        }

        int exy = a[1] - a[0];
        int exz = a[2] - a[0];
        int ezy = a[1] - a[2];
        //
        //  We start at the starting point.
        //
        for (i = 0; i < DIM_NUM; i++)
        {
            v[i + 0 * DIM_NUM] = v1[i];
        }

        for (j = 1; j < n; j++)
        {
            for (i = 0; i < DIM_NUM; i++)
            {
                v[i + j * DIM_NUM] = v[i + (j - 1) * DIM_NUM];
            }

            switch (exy)
            {
                case < 0 when exz < 0:
                    v[0 + j * DIM_NUM] += s[0];
                    exy += 2 * a[1];
                    exz += 2 * a[2];
                    break;
                case < 0:
                    v[2 + j * DIM_NUM] += s[2];
                    exz -= 2 * a[0];
                    ezy += 2 * a[1];
                    break;
                default:
                {
                    switch (ezy)
                    {
                        case < 0:
                            v[2 + j * DIM_NUM] += s[2];
                            exz -= 2 * a[0];
                            ezy += 2 * a[1];
                            break;
                        default:
                            v[1 + j * DIM_NUM] += s[1];
                            exy -= 2 * a[0];
                            ezy -= 2 * a[2];
                            break;
                    }

                    break;
                }
            }
        }
    }

    public static void voxels_region_3d(int list_max, int nx, int ny, int nz, ref int[] ishow,
            ref int list_num, ref int[] list, ref int region_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VOXELS_REGION_3D arranges a set of voxels into contiguous regions in 3D.
        //
        //  Discussion:
        //
        //    On input, the ISHOW array contains zero and nonzero values.  The nonzero
        //    values are taken to be active voxels.  On output, the zero voxels remain
        //    zero, and all the active voxels have been assigned a value which now
        //    indicates membership in a region, or group of contiguous voxels.
        //
        //    On output, the array LIST contains information about the regions.
        //    The last used element of LIST is LIST_NUM.
        //
        //    The number of elements in region REGION_NUM is NELEM = LIST(LIST_NUM).
        //    The (I,J,K) indices of the last element in this region are in
        //    LIST(LIST_NUM-3) through LIST(LIST_NUM-1), and the first element is
        //    listed in LIST(LIST_NUM-3*NELEM), LIST(LIST_NUM-3*NELEM+1),
        //    LIST(LIST_NUM-3*NELEM+2).
        //
        //    The number of elements in REGION_NUM-1 is listed in LIST(LIST_NUM-3*NELEM-1),
        //    and the (I,J,K) indices of the these elements are listed there.
        //
        //  Picture:
        //
        //    Input:
        //
        //      0  2  0  0 17  0  3
        //      0  0  3  0  1  0  4
        //      1  0  4  8  8  0  7
        //      3  0  6 45  0  0  0
        //      3 17  0  5  9  2  5
        //
        //    Output:
        //
        //      0  1  0  0  2  0  3
        //      0  0  2  0  2  0  3
        //      4  0  2  2  2  0  3
        //      4  0  2  2  0  0  0
        //      4  4  0  2  2  2  2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int ISHOW[NX*NY*NZ].  On input, the only significance to
        //    the entries is whether they are zero or nonzero.  On output, the nonzero
        //    entries have now been revalued so that contiguous entries have the same
        //    value, indicating a grouping into a region.
        //
        //    Output, int LIST[LIST_MAX], contains, in stack form, a list
        //    of the indices of the elements in each region.
        //
        //    Input, int LIST_MAX, the maximum length of the array used to
        //    list the elements of the regions.
        //
        //    Output, int *LIST_NUM, the number of entries of LIST that were used.
        //    However, if LIST_MAX < LIST_NUM, then there was not enough space in
        //    LIST to store the data properly, and LIST should not be used,
        //    although the data in ISHOW should be correct.
        //
        //    Output, int *REGION_NUM, the number of regions discovered.
        //
        //    Input, int NX, NY, NZ, the number of voxels in the X, Y and
        //    Z directions.
        //
    {
        const int STACK_MAX = 100;

        int i;
        int j;
        int k;
        int[] stack = new int[STACK_MAX];
        //
        //  Reset all nonzero entries of ISHOW to -1.
        //
        for (k = 0; k < nz; k++)
        {
            for (j = 0; j < ny; j++)
            {
                for (i = 0; i < nx; i++)
                {
                    if (ishow[i + nx * (j + ny * k)] != 0)
                    {
                        ishow[i + nx * (j + ny * k)] = -1;
                    }
                }
            }
        }

        //
        //  Start the number of items in the region list at 0.
        //
        list_num = 0;
        //
        //  Start the number of regions at 0.
        //
        region_num = 0;
        //
        //  The stack begins empty.
        //
        int nstack = 0;
        //
        //  Search for an unused "ON" voxel from which we can "grow" a new region.
        //
        for (k = 1; k <= nz; k++)
        {
            for (j = 1; j <= ny; j++)
            {
                for (i = 1; i <= nx; i++)
                {
                    switch (ishow[i - 1 + nx * (j - 1 + ny * (k - 1))])
                    {
                        //
                        //  We found a voxel that is "ON", and does not belong to any region.
                        //
                        case -1:
                        {
                            //
                            //  Increase the number of regions.
                            //
                            region_num += 1;
                            //
                            //  Add this voxel to the region.
                            //
                            ishow[i - 1 + nx * (j - 1 + ny * (k - 1))] = region_num;
                            //
                            //  Add this voxel to the stack.
                            //
                            if (STACK_MAX < nstack + 4)
                            {
                                Console.WriteLine("");
                                Console.WriteLine("VOXELS_REGION - Fatal error!");
                                Console.WriteLine("  The internal stack overflowed.");
                                Console.WriteLine("  The algorithm has failed.");
                                return;
                            }

                            stack[nstack + 1 - 1] = i;
                            stack[nstack + 2 - 1] = j;
                            stack[nstack + 3 - 1] = k;
                            stack[nstack + 4 - 1] = 1;

                            nstack += 4;
                            //
                            //  Add this voxel to the description of the region.
                            //
                            int nelements = 1;

                            if (list_num + 3 <= list_max)
                            {
                                list[list_num + 1 - 1] = i;
                                list[list_num + 2 - 1] = j;
                                list[list_num + 3 - 1] = k;
                            }

                            list_num += 3;

                            for (;;)
                            {
                                //
                                //  Find all neighbors of BASE that are "ON" but unused.
                                //  Mark them as belonging to this region, and stack their indices.
                                //
                                int ibase = stack[nstack - 3 - 1];
                                int jbase = stack[nstack - 2 - 1];
                                int kbase = stack[nstack - 1 - 1];

                                int ilo = Math.Max(ibase - 1, 1);
                                int ihi = Math.Min(ibase + 1, nx);
                                int jlo = Math.Max(jbase - 1, 1);
                                int jhi = Math.Min(jbase + 1, ny);
                                int klo = Math.Max(kbase - 1, 1);
                                int khi = Math.Min(kbase + 1, nz);

                                int nabes = 0;

                                int k2;
                                for (k2 = klo; k2 <= khi; k2++)
                                {
                                    int j2;
                                    for (j2 = jlo; j2 <= jhi; j2++)
                                    {
                                        int i2;
                                        for (i2 = ilo; i2 <= ihi; i2++)
                                        {
                                            switch (ishow[i2 - 1 + nx * (j2 - 1 + ny * (k2 - 1))])
                                            {
                                                //
                                                //  We found a neighbor to our current search point, which is "ON" and unused.
                                                //
                                                case -1:
                                                {
                                                    //
                                                    //  Increase the number of neighbors.
                                                    //
                                                    nabes += 1;
                                                    //
                                                    //  Mark the neighbor as belonging to the region.
                                                    //
                                                    ishow[i2 - 1 + nx * (j2 - 1 + ny * (k2 - 1))] = region_num;
                                                    //
                                                    //  Add the neighbor to the stack.
                                                    //
                                                    if (STACK_MAX < nstack + 3)
                                                    {
                                                        Console.WriteLine("");
                                                        Console.WriteLine("VOXELS_REGION - Fatal error!");
                                                        Console.WriteLine("  The internal stack overflowed.");
                                                        Console.WriteLine("  The algorithm has failed.");
                                                        return;
                                                    }

                                                    stack[nstack + 1 - 1] = i2;
                                                    stack[nstack + 2 - 1] = j2;
                                                    stack[nstack + 3 - 1] = k2;

                                                    nstack += 3;
                                                    //
                                                    //  Add the neighbor to the description of the region.
                                                    //
                                                    nelements += 1;

                                                    if (list_num + 3 <= list_max)
                                                    {
                                                        list[list_num + 1 - 1] = i2;
                                                        list[list_num + 2 - 1] = j2;
                                                        list[list_num + 3 - 1] = k2;
                                                    }

                                                    list_num += 3;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }

                                switch (nabes)
                                {
                                    //
                                    //  If any new neighbors were found, take the last one as the basis
                                    //  for a deeper search.
                                    //
                                    case > 0 when STACK_MAX < nstack + 1:
                                        Console.WriteLine("");
                                        Console.WriteLine("VOXELS_REGION - Fatal error!");
                                        Console.WriteLine("  The internal stack overflowed.");
                                        Console.WriteLine("  The algorithm has failed.");
                                        return;
                                    case > 0:
                                        stack[nstack + 1 - 1] = nabes;
                                        nstack += 1;
                                        continue;
                                }

                                //
                                //  If the current search point had no new neighbors, drop it from the stack.
                                //
                                int ncan = stack[nstack - 1] - 1;
                                nstack -= 3;
                                stack[nstack - 1] = ncan;
                                switch (stack[nstack - 1])
                                {
                                    //
                                    //  If there are still any unused candidates at this level, take the
                                    //  last one as the basis for a deeper search.
                                    //
                                    case > 0:
                                        continue;
                                }

                                //
                                //  If there are no more unused candidates at this level, then we need
                                //  to back up a level in the stack.  If there are any candidates at
                                //  that earlier level, then we can still do more searching.
                                //
                                nstack -= 1;

                                if (nstack <= 0)
                                {
                                    break;
                                }
                            }

                            //
                            //  If we have exhausted the stack, we have completed this region.
                            //  Tag the number of elements to the end of the region description list.
                            //
                            list_num += 1;
                            if (list_num <= list_max)
                            {
                                list[list_num - 1] = nelements;
                            }

                            break;
                        }
                    }
                }
            }
        }

        //
        //  Print some warnings.
        //
        if (list_max >= list_num)
        {
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("VOXELS_REGION - Warning!");
        Console.WriteLine("  LIST_MAX was too small to list the regions.");
        Console.WriteLine("  Do not try to use the LIST array!");
        Console.WriteLine("  The ISHOW data is OK, however.");
    }

    public static void voxels_step_3d(int[] v1, int[] v2, int inc, int jnc, int knc,
            ref int[] v3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VOXELS_STEP_3D computes voxels along a line from a given point in 3D.
        //
        //  Discussion:
        //
        //    If you input INC = JNC = KNC, then no movement is possible,
        //    and none is made.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int V1[3], the coordinates of the base voxel from
        //    which the line begins.
        //
        //    Input, int V2[3], the coordinates of the current voxel on the
        //    line.  For the first call, these might be equal to V1.
        //
        //    Input, int INC, JNC, KNC, the increments to the voxels.
        //    These values define the direction along which the line proceeds.
        //    However, the voxels on the line will typically be incremented
        //    by a fractional value of the vector (INC,JNC,KNC), and the
        //    result is essentially rounded.
        //
        //    Output, int V3[3], the coordinates of the next voxel along
        //    the line.
        //
    {
        const int DIM_NUM = 3;

        double alpha = 0.0;
        double alphai = 0.0;
        double alphaj = 0.0;
        double alphak = 0.0;

        typeMethods.i4vec_copy(DIM_NUM, v2, ref v3);
        switch (inc)
        {
            //
            //  Assuming for the moment that (I,J,K) can take on real values,
            //  points on the line have the form:
            //
            //    I = V1[0] + alpha * inc
            //    J = V1[1] + alpha * jnc
            //    K = V1[2] + alpha * knc
            //
            case 0 when jnc == 0 && knc == 0:
                return;
        }

        alpha = 0.0;
        alphai = inc switch
        {
            //
            //  Compute the smallest ALPHA that will change I2, J2 or K2 by +-0.5.
            //
            > 0 => (v2[0] - v1[0] + 0.5) / inc,
            < 0 => (v2[0] - v1[0] - 0.5) / inc,
            _ => typeMethods.r8_huge()
        };

        alphaj = jnc switch
        {
            > 0 => (v2[1] - v1[1] + 0.5) / jnc,
            < 0 => (v2[1] - v1[1] - 0.5) / jnc,
            _ => typeMethods.r8_huge()
        };

        switch (knc)
        {
            case > 0:
                alphak = (v2[2] - v1[2] + 0.5) / knc;
                break;
            case < 0:
                alphak = (v2[2] - v1[2] - 0.5) / knc;
                break;
            default:
                alphaj = typeMethods.r8_huge();
                break;
        }

        //
        //  The ALPHA of smallest positive magnitude represents the closest next voxel.
        //
        alpha = typeMethods.r8_huge();

        switch (alphai)
        {
            case > 0.0:
                alpha = Math.Min(alpha, alphai);
                break;
        }

        alpha = alphak switch
        {
            > 0.0 => Math.Min(alpha, alphak),
            _ => alphaj switch
            {
                > 0.0 => Math.Min(alpha, alphaj),
                _ => alpha
            }
        };

        //
        //  Move to the new voxel.  Whichever index just made the half
        //  step must be forced to take a whole step.
        //
        if (Math.Abs(alpha - alphai) <= double.Epsilon)
        {
            v3[0] = v2[0] + typeMethods.i4_sign(inc);
            v3[1] = v1[1] + (int) typeMethods.r8_nint(alpha * jnc);
            v3[2] = v1[2] + (int) typeMethods.r8_nint(alpha * knc);
        }
        else if (Math.Abs(alpha - alphaj) <= double.Epsilon)
        {
            v3[0] = v1[0] + (int) typeMethods.r8_nint(alpha * inc);
            v3[1] = v2[1] + typeMethods.i4_sign(jnc);
            v3[2] = v1[2] + (int) typeMethods.r8_nint(alpha * knc);
        }
        else if (Math.Abs(alpha - alphak) <= double.Epsilon)
        {
            v3[0] = v1[0] + (int) typeMethods.r8_nint(alpha * inc);
            v3[1] = v1[1] + (int) typeMethods.r8_nint(alpha * jnc);
            v3[2] = v2[2] + typeMethods.i4_sign(knc);
        }
    }
}