namespace Burkardt.TriangulationNS;

public static partial class Triangulation
{
    public static void triangle_order6_reference_to_physical(double[] t, int n,
            double[] ref_, ref double[] phy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_ORDER6_REFERENCE_TO_PHYSICAL maps reference points to physical points.
        //
        //  Discussion:
        //
        //    Given the vertices of an order 6 physical triangle and a point
        //    (XSI,ETA) in the reference triangle, the routine computes the value
        //    of the corresponding image point (X,Y) in physical space.
        //
        //    The mapping from (XSI,ETA) to (X,Y) has the form:
        //
        //      X(ETA,XSI) = A1 * XSI**2 + B1 * XSI*ETA + C1 * ETA**2
        //                 + D1 * XSI    + E1 * ETA     + F1
        //
        //      Y(ETA,XSI) = A2 * XSI**2 + B2 * XSI*ETA + C2 * ETA**2
        //                 + D2 * XSI    + E2 * ETA     + F2
        //
        //  Reference Element T6:
        //
        //    |
        //    1  3
        //    |  |.
        //    |  | .
        //    S  6  5
        //    |  |   .
        //    |  |    .
        //    0  1--4--2
        //    |
        //    +--0--R--1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*6], the coordinates of the vertices.
        //    The vertices are assumed to be the images of (0,0), (1,0),
        //    (0,1),(1/2,0), (1/2,1/2) and (0,1/2) respectively.
        //
        //    Input, integer N, the number of points to transform.
        //
        //    Input, double REF[2*N], points in the reference triangle.
        //
        //    Output, double PHY[2*N], corresponding points in the
        //    physical triangle.
        //
    {
        double[] a = new double[2];
        double[] b = new double[2];
        double[] c = new double[2];
        double[] d = new double[2];
        double[] e = new double[2];
        double[] f = new double[2];
        int i;
        int j;

        for (i = 0; i < 2; i++)
        {
            a[i] = 2.0 * t[i + 0 * 2] + 2.0 * t[i + 1 * 2]
                   - 4.0 * t[i + 3 * 2];

            b[i] = 4.0 * t[i + 0 * 2]
                - 4.0 * t[i + 3 * 2] + 4.0 * t[i + 4 * 2] - 4.0 * t[i + 5 * 2];

            c[i] = 2.0 * t[i + 0 * 2] + 2.0 * t[i + 2 * 2]
                   - 4.0 * t[i + 5 * 2];

            d[i] = -3.0 * t[i + 0 * 2] - t[i + 1 * 2]
                   + 4.0 * t[i + 3 * 2];

            e[i] = -3.0 * t[i + 0 * 2] - t[i + 2 * 2]
                   + 4.0 * t[i + 5 * 2];
            f[i] = t[i + 0 * 2];

        }

        for (i = 0; i < 2; i++)
        {
            for (j = 0; j < n; j++)
            {
                phy[i + j * 2] = a[i] * ref_[
                                     0 + j * 2] * ref_[
                                     0 + j * 2]
                                 +b[i] * ref_[
                                     0 + j * 2] * ref_[
                                     1 + j * 2]
                                 +c[i] * ref_[
                                     1 + j * 2] * ref_[
                                     1 + j * 2]
                                 +d[i] * ref_[
                                     0 + j * 2]
                                 +e[i] * ref_[
                                     1 + j * 2]
                                 +f[i];
            }
        }
    }
}