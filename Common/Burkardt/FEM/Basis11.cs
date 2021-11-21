using System;
using System.Globalization;
using Burkardt.Types;

namespace Burkardt.FEM;

public static class Basis11
{
    public static void basis_one_t3(double[] t, int i, double[] p, ref double qi,
            ref double dqidx, ref double dqidy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_ONE_T3 evaluates basis functions for a linear triangular element.
        //
        //  Discussion:
        //
        //    The routine is given the coordinates of the nodes of a triangle.
        //
        //           3
        //          / .
        //         /   .
        //        /     .
        //       1-------2
        //
        //    It evaluates the linear basis function Q(I)(X,Y) associated with
        //    node I, which has the property that it is a linear function
        //    which is 1 at node I and zero at the other two nodes.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 January 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the coordinates of the nodes.
        //
        //    Input, int I, the index of the desired basis function.
        //    I should be between 1 and 3.
        //
        //    Input, double P[2], the coordinates of the point where
        //    the basis function is to be evaluated.
        //
        //    Output, double *QI, *DQIDX, *DQIDY, the value of the I-th basis function
        //    and its X and Y derivatives.
        //
    {
        double area = t[0 + 0 * 2] * (t[1 + 1 * 2] - t[1 + 2 * 2])
                      + t[0 + 1 * 2] * (t[1 + 2 * 2] - t[1 + 0 * 2])
                      + t[0 + 2 * 2] * (t[1 + 0 * 2] - t[1 + 1 * 2]);

        switch (area)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("BASIS_ONE_T3 - Fatal error!");
                Console.WriteLine("  Element has zero area.");
                Console.WriteLine("  Area = " + area + "");
                Console.WriteLine("");
                Console.WriteLine("  Node 1: ( " + t[0 + 0 * 2] + ", " + t[1 + 0 * 2] + " )");
                Console.WriteLine("  Node 2: ( " + t[0 + 1 * 2] + ", " + t[1 + 1 * 2] + " )");
                Console.WriteLine("  Node 3: ( " + t[0 + 2 * 2] + ", " + t[1 + 2 * 2] + " )");
                return;
        }

        switch (i)
        {
            case < 1:
            case > 3:
                Console.WriteLine("");
                Console.WriteLine("BASIS_ONE_T3 - Fatal error!");
                Console.WriteLine("  Basis index I is not between 1 and 3.");
                Console.WriteLine("  I = " + i + "");
                return;
        }

        int ip1 = typeMethods.i4_wrap(i + 1, 1, 3);
        int ip2 = typeMethods.i4_wrap(i + 2, 1, 3);

        qi = ((t[0 + (ip2 - 1) * 2] - t[0 + (ip1 - 1) * 2])
              * (p[1] - t[1 + (ip1 - 1) * 2])
              - (t[1 + (ip2 - 1) * 2] - t[1 + (ip1 - 1) * 2])
              * (p[0] - t[0 + (ip1 - 1) * 2])) / area;

        dqidx = -(t[1 + (ip2 - 1) * 2] - t[1 + (ip1 - 1) * 2]) / area;
        dqidy = (t[0 + (ip2 - 1) * 2] - t[0 + (ip1 - 1) * 2]) / area;
    }

    public static void basis_11_t3(double[] t, int i, double[] p, ref double qi,
            ref double dqidx, ref double dqidy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_11_T3: one basis at one point for a T3 element.
        //
        //  Discussion:
        //
        //    The routine is given the coordinates of the nodes of a triangle.
        //
        //           3
        //          . .
        //         .   .
        //        .     .
        //       1-------2
        //
        //    It evaluates the linear basis function Q(I)(X,Y) associated with
        //    node I, which has the property that it is a linear function
        //    which is 1 at node I and zero at the other two nodes.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 January 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double[] t, the coordinates of the nodes.
        //
        //    Input, int I, the index of the desired basis function.
        //    I should be between 1 and 3.
        //
        //    Input, double[] p, the coordinates of the point where 
        //    the basis function is to be evaluated.
        //
        //    Output, ref double QI, *DQIDX, *DQIDY, the value of the I-th basis function
        //    and its X and Y derivatives.
        //
    {
        double area = t[0 + 0 * 2] * (t[1 + 1 * 2] - t[1 + 2 * 2])
                      + t[0 + 1 * 2] * (t[1 + 2 * 2] - t[1 + 0 * 2])
                      + t[0 + 2 * 2] * (t[1 + 0 * 2] - t[1 + 1 * 2]);

        switch (area)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("BASIS_11_T3 - Fatal error!");
                Console.WriteLine("  Element has zero area.");
                Console.WriteLine("  Area = " + area + "");
                return;
        }

        switch (i)
        {
            case < 1:
            case > 3:
                Console.WriteLine("");
                Console.WriteLine("BASIS_11_T3 - Fatal error!");
                Console.WriteLine("  Basis index I is not between 1 and 3.");
                Console.WriteLine("  I = " + i + "");
                return;
        }

        int ip1 = typeMethods.i4_wrap(i + 1, 1, 3);
        int ip2 = typeMethods.i4_wrap(i + 2, 1, 3);

        qi = ((t[0 + (ip2 - 1) * 2] - t[0 + (ip1 - 1) * 2])
              * (p[1] - t[1 + (ip1 - 1) * 2])
              - (t[1 + (ip2 - 1) * 2] - t[1 + (ip1 - 1) * 2])
              * (p[0] - t[0 + (ip1 - 1) * 2])) / area;

        dqidx = -(t[1 + (ip2 - 1) * 2] - t[1 + (ip1 - 1) * 2]) / area;
        dqidy = (t[0 + (ip2 - 1) * 2] - t[0 + (ip1 - 1) * 2]) / area;

    }

    public static void basis_11_t3_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_11_T3_TEST verifies BASIS_11_T3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 January 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    None
        //
    {
        const int NODE_NUM = 3;

        double dqjdx = 0;
        double dqjdy = 0;
        int i;
        int j;
        double qj = 0;
        double[] p = new double[2];
        double[] t =
        {
            2.0, 0.0,
            4.0, 3.0,
            0.0, 4.0
        };

        Console.WriteLine("");
        Console.WriteLine("BASIS_11_T3_TEST:");
        Console.WriteLine("  Verify basis functions for element T3.");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + NODE_NUM + "");

        Console.WriteLine("");
        Console.WriteLine("  Physical Nodes:");
        Console.WriteLine("");
        for (j = 0; j < NODE_NUM; j++)
        {
            Console.WriteLine("  "
                              + t[0 + j * 2].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + t[1 + j * 2].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  The basis function values at basis nodes");
        Console.WriteLine("  should form the identity matrix.");
        Console.WriteLine("");

        for (i = 0; i < NODE_NUM; i++)
        {
            p[0] = t[0 + i * 2];
            p[1] = t[1 + i * 2];

            string cout = "";
            for (j = 0; j < NODE_NUM; j++)
            {
                basis_11_t3(t, j + 1, p, ref qj, ref dqjdx, ref dqjdy);
                cout += "  " + qj.ToString(CultureInfo.InvariantCulture).PadLeft(10);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  The X and Y derivatives should sum to 0.");
        Console.WriteLine("");
        Console.WriteLine("  dPhidX sum, dPhidY sum:");
        Console.WriteLine("");

        for (i = 0; i < NODE_NUM; i++)
        {
            p[0] = t[0 + i * 2];
            p[1] = t[1 + i * 2];

            double sum_x = 0.0;
            double sum_y = 0.0;
            for (j = 0; j < NODE_NUM; j++)
            {
                basis_11_t3(t, j + 1, p, ref qj, ref dqjdx, ref dqjdy);
                sum_x += dqjdx;
                sum_y += dqjdy;
            }

            Console.WriteLine("  "
                              + sum_x.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + sum_y.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

    }

    public static void basis_11_t4(double[] t, int i, double[] p, ref double phi,
            ref double dphidx, ref double dphidy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_MN_T4: one basis at one point for a T4 element.
        //
        //  Discussion:
        //
        //    The T4 element is the cubic bubble triangle.
        //
        //    The routine is given the coordinates of the vertices of a triangle.
        //    It works directly with these coordinates, and does not refer to a 
        //    reference element.
        //
        //    The sides of the triangle DO NOT have to lie along a coordinate
        //    axis.
        //
        //    The routine evaluates the basis functions associated with each vertex,
        //    and their derivatives with respect to X and Y.
        //
        //  Physical Element T4: 
        //       
        //            3
        //           . .
        //          .   .
        //         .  4  .
        //        .       .
        //       1---------2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*4], the coordinates of the vertices
        //    of the triangle, and the coordinates of the centroid.  
        //    It is common to list the first three points in counter clockwise
        //    order.
        //
        //    Input, int I, the index of the basis function.
        //
        //    Input, double[] p, the point where the basis function
        //    is to be evaluated.
        //
        //    Output, ref double PHI, the value of the basis function
        //    at the evaluation point.
        //
        //    Output, ref double DPHIDX, *DPHIDY, the value of the 
        //    derivatives at the evaluation point.
        //
        //  Local parameters:
        //
        //    Local, double AREA, is (twice) the area of the triangle.
        //
    {
        double[] dpsidx = new double[4];
        double[] dpsidy = new double[4];
        int j;
        double[] psi = new double[4];

        double area = t[0 + 0 * 2] * (t[1 + 1 * 2] - t[1 + 2 * 2])
                      + t[0 + 1 * 2] * (t[1 + 2 * 2] - t[1 + 0 * 2])
                      + t[0 + 2 * 2] * (t[1 + 0 * 2] - t[1 + 1 * 2]);

        psi[0] = (t[0 + 2 * 2] - t[0 + 1 * 2]) * (p[1] - t[1 + 1 * 2])
                 - (t[1 + 2 * 2] - t[1 + 1 * 2]) * (p[0] - t[0 + 1 * 2]);
        dpsidx[0] = -(t[1 + 2 * 2] - t[1 + 1 * 2]);
        dpsidy[0] = t[0 + 2 * 2] - t[0 + 1 * 2];

        psi[1] = (t[0 + 0 * 2] - t[0 + 2 * 2]) * (p[1] - t[1 + 2 * 2])
                 - (t[1 + 0 * 2] - t[1 + 2 * 2]) * (p[0] - t[0 + 2 * 2]);
        dpsidx[1] = -(t[1 + 0 * 2] - t[1 + 2 * 2]);
        dpsidy[1] = t[0 + 0 * 2] - t[0 + 2 * 2];

        psi[2] = (t[0 + 1 * 2] - t[0 + 0 * 2]) * (p[1] - t[1 + 0 * 2])
                 - (t[1 + 1 * 2] - t[1 + 0 * 2]) * (p[0] - t[0 + 0 * 2]);
        dpsidx[2] = -(t[1 + 1 * 2] - t[1 + 0 * 2]);
        dpsidy[2] = t[0 + 1 * 2] - t[0 + 0 * 2];
        //
        //  Normalize the first three functions.
        //
        for (j = 0; j < 3; j++)
        {
            psi[j] /= area;
            dpsidx[j] /= area;
            dpsidy[j] /= area;
        }

        //
        //  Compute the cubic bubble function.
        //
        psi[3] = 27.0 * psi[0] * psi[1] * psi[2];

        dpsidx[3] = 27.0 * (
            dpsidx[0] * psi[1] * psi[2]
            + psi[0] * dpsidx[1] * psi[2]
            + psi[0] * psi[1] * dpsidx[2]);

        dpsidy[3] = 27.0 * (
            dpsidy[0] * psi[1] * psi[2]
            + psi[0] * dpsidy[1] * psi[2]
            + psi[0] * psi[1] * dpsidy[2]);
        //
        //  Subtract 1/3 of the cubic bubble function from each of the three linears.
        //
        for (j = 0; j < 3; j++)
        {
            psi[j] -= psi[3] / 3.0;
            dpsidx[j] -= dpsidx[3] / 3.0;
            dpsidy[j] -= dpsidy[3] / 3.0;
        }

        phi = psi[i - 1];
        dphidx = dpsidx[i - 1];
        dphidy = dpsidy[i - 1];
    }

    public static void basis_11_t4_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_11_T4_TEST verifies BASIS_11_T4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    None
        //
    {
        const int NODE_NUM = 4;

        double dqjdx = 0;
        double dqjdy = 0;
        int i;
        int j;
        double qj = 0;
        double[] p = new double[2];
        double[] t =
        {
            2.0, 0.0,
            4.0, 3.0,
            0.0, 4.0,
            0.0, 0.0
        };
        //
        //  The node associated with the fourth basis function is the centroid.
        //
        t[0 + 3 * 2] = (t[0 + 0 * 2] + t[0 + 1 * 2] + t[0 + 2 * 2]) / 3.0;
        t[1 + 3 * 2] = (t[1 + 0 * 2] + t[1 + 1 * 2] + t[1 + 2 * 2]) / 3.0;

        Console.WriteLine("");
        Console.WriteLine("BASIS_11_T4_TEST:");
        Console.WriteLine("  Verify basis functions for element T4.");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + NODE_NUM + "");

        Console.WriteLine("");
        Console.WriteLine("  Physical Nodes:");
        Console.WriteLine("");
        for (j = 0; j < NODE_NUM; j++)
        {
            Console.WriteLine("  "
                              + t[0 + j * 2].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + t[1 + j * 2].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  The basis function values at basis nodes");
        Console.WriteLine("  should form the identity matrix.");
        Console.WriteLine("");

        for (i = 0; i < NODE_NUM; i++)
        {
            p[0] = t[0 + i * 2];
            p[1] = t[1 + i * 2];

            string cout = "";
            for (j = 0; j < NODE_NUM; j++)
            {
                basis_11_t4(t, j + 1, p, ref qj, ref dqjdx, ref dqjdy);
                cout += "  " + qj.ToString(CultureInfo.InvariantCulture).PadLeft(10);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  The X and Y derivatives should sum to 0.");
        Console.WriteLine("");
        Console.WriteLine("  dPhidX sum, dPhidY sum:");
        Console.WriteLine("");

        for (i = 0; i < NODE_NUM; i++)
        {
            p[0] = t[0 + i * 2];
            p[1] = t[1 + i * 2];

            double sum_x = 0.0;
            double sum_y = 0.0;
            for (j = 0; j < NODE_NUM; j++)
            {
                basis_11_t4(t, j + 1, p, ref qj, ref dqjdx, ref dqjdy);
                sum_x += dqjdx;
                sum_y += dqjdy;
            }

            Console.WriteLine("  "
                              + sum_x.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + sum_y.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

    }

    public static void basis_11_t6(double[] t, int i, double[] p, ref double bi,
            ref double dbidx, ref double dbidy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_11_T6: one basis at one point for the T6 element.
        //
        //  Discussion:
        //
        //    The routine is given the coordinates of the nodes of a triangle. 
        //       
        //           3
        //          . .
        //         6   5
        //        .     .
        //       1---4---2
        //
        //    It evaluates the quadratic basis function B(I)(X,Y) associated with
        //    node I, which has the property that it is a quadratic function
        //    which is 1 at node I and zero at the other five nodes.
        //
        //    This routine assumes that the sides of the triangle are straight,
        //    so that the midside nodes fall on the line between two vertices.
        //
        //    This routine relies on the fact that each basis function can be
        //    written as the product of two linear factors, which are easily
        //    computed and normalized.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 February 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*6], the coordinates of the nodes.
        //
        //    Input, int I, the index of the desired basis function.
        //    I should be between 1 and 6.
        //
        //    Input, double[] p, the coordinates of a point at which the basis
        //    function is to be evaluated.
        //
        //    Output, ref double BI, *DBIDX, *DBIDY, the values of the basis function
        //    and its X and Y derivatives.
        //
    {
        int j1;
        int j2;
        int k1;
        int k2;

        switch (i)
        {
            case < 1:
            case > 6:
                Console.WriteLine("");
                Console.WriteLine("BASIS_11_T6 - Fatal error!");
                Console.WriteLine("  Basis index I is not between 1 and 6.");
                Console.WriteLine("  I = " + i + "");
                return;
            //
            //  Determine the pairs of nodes.
            //
            case <= 3:
                j1 = typeMethods.i4_wrap(i + 1, 1, 3);
                j2 = typeMethods.i4_wrap(i + 2, 1, 3);
                k1 = i + 3;
                k2 = typeMethods.i4_wrap(i + 5, 4, 6);
                break;
            default:
                j1 = i - 3;
                j2 = typeMethods.i4_wrap(i - 3 + 2, 1, 3);
                k1 = typeMethods.i4_wrap(i - 3 + 1, 1, 3);
                k2 = typeMethods.i4_wrap(i - 3 + 2, 1, 3);
                break;
        }

        //
        //  For C++ indexing, it is helpful to knock the indices down by one.
        //
        i -= 1;
        j1 -= 1;
        j2 -= 1;
        k1 -= 1;
        k2 -= 1;
        //
        //  Evaluate the two linear factors GF and HF, 
        //  and their normalizers GN and HN.
        //
        double gf = (p[0] - t[0 + j1 * 2]) * (t[1 + j2 * 2] - t[1 + j1 * 2])
                    - (t[0 + j2 * 2] - t[0 + j1 * 2]) * (p[1] - t[1 + j1 * 2]);

        double gn = (t[0 + i * 2] - t[0 + j1 * 2]) * (t[1 + j2 * 2] - t[1 + j1 * 2])
                    - (t[0 + j2 * 2] - t[0 + j1 * 2]) * (t[1 + i * 2] - t[1 + j1 * 2]);

        double hf = (p[0] - t[0 + k1 * 2]) * (t[1 + k2 * 2] - t[1 + k1 * 2])
                    - (t[0 + k2 * 2] - t[0 + k1 * 2]) * (p[1] - t[1 + k1 * 2]);

        double hn = (t[0 + i * 2] - t[0 + k1 * 2]) * (t[1 + k2 * 2] - t[1 + k1 * 2])
                    - (t[0 + k2 * 2] - t[0 + k1 * 2]) * (t[1 + i * 2] - t[1 + k1 * 2]);
        //
        //  Construct the basis function and its derivatives.
        //
        bi = gf / gn
             * (hf / hn);

        dbidx = (t[1 + j2 * 2] - t[1 + j1 * 2]) / gn
                * (hf / hn)
                + gf / gn
                * ((t[1 + k2 * 2] - t[1 + k1 * 2]) / hn);

        dbidy = -((t[0 + j2 * 2] - t[0 + j1 * 2]) / gn)
                * (hf / hn)
                - gf / gn
                * ((t[0 + k2 * 2] - t[0 + k1 * 2]) / hn);
    }

    public static void basis_11_t6_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BASIS_11_T6_TEST verifies BASIS_11_T6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 February 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    None
        //
    {
        const int NODE_NUM = 6;

        double[] dphidx = new double[NODE_NUM * NODE_NUM];
        double[] dphidy = new double[NODE_NUM * NODE_NUM];
        int i;
        int j;
        double[] p = new double[2];
        double[] phi = new double[NODE_NUM * NODE_NUM];
        double[] t =
        {
            2.0, 0.0,
            4.0, 3.0,
            0.0, 4.0,
            3.0, 1.5,
            2.0, 3.5,
            1.0, 2.0
        };
        double v1 = 0;
        double v2 = 0;
        double v3 = 0;

        Console.WriteLine("");
        Console.WriteLine("BASIS_11_T6_TEST:");
        Console.WriteLine("  Verify basis functions for element T6.");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + NODE_NUM + "");

        Console.WriteLine("");
        Console.WriteLine("  Physical Nodes:");
        Console.WriteLine("");
        for (j = 0; j < NODE_NUM; j++)
        {
            Console.WriteLine("  "
                              + j.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + t[0 + j * 2].ToString(CultureInfo.InvariantCulture).PadLeft(7) + "  "
                              + t[1 + j * 2].ToString(CultureInfo.InvariantCulture).PadLeft(7) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  The basis function values at basis nodes");
        Console.WriteLine("  should form the identity matrix.");
        Console.WriteLine("");
        for (i = 1; i <= NODE_NUM; i++)
        {
            for (j = 0; j < NODE_NUM; j++)
            {
                p[0] = t[0 + j * 2];
                p[1] = t[1 + j * 2];

                basis_11_t6(t, i, p, ref v1, ref v2, ref v3);

                phi[i - 1 + j * NODE_NUM] = v1;
                dphidx[i - 1 + j * NODE_NUM] = v2;
                dphidy[i - 1 + j * NODE_NUM] = v3;
            }
        }

        for (i = 0; i < NODE_NUM; i++)
        {
            string cout = "";
            for (j = 0; j < NODE_NUM; j++)
            {
                cout += "  " + phi[i + j * NODE_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(7);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  The X and Y derivatives should sum to 0.");
        Console.WriteLine("");
        Console.WriteLine("  dPhidX sum, dPhidY sum:");
        Console.WriteLine("");
        for (j = 0; j < NODE_NUM; j++)
        {
            double sum_x = 0.0;
            for (i = 0; i < NODE_NUM; i++)
            {
                sum_x += dphidx[i + j * NODE_NUM];
            }

            double sum_y = 0.0;
            for (i = 0; i < NODE_NUM; i++)
            {
                sum_y += dphidy[i + j * NODE_NUM];
            }

            Console.WriteLine("  "
                              + sum_x.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                              + sum_y.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

    }
}