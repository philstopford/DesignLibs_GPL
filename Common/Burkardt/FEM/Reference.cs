using System;
using Burkardt.Types;
using Burkardt.Uniform;
using Burkardt.Table;

namespace Burkardt.FEM
{
    public static class Reference
    {

        public static void reference_sample(string code, ref int seed, ref double r, ref double s)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    REFERENCE_SAMPLE samples a reference element.
            //
            //  Discussion:
            //
            //    The routine either samples the unit triangle or the unit square.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 September 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string CODE, identifies the element desired.
            //    Legal values include "Q4", "Q8", "Q9", "Q12", "Q16", "QL", "T3", 
            //    "T4", "T6" and "T10".
            //
            //    Input/output, int *SEED, a seed for the random number generator.
            //
            //    Output, ref double R, *S, a random point in the reference element.
            //
        {
            if (code == "Q4" ||
                code == "Q8" ||
                code == "Q9" ||
                code == "Q12" ||
                code == "Q16" ||
                code == "QL")
            {
                r = UniformRNG.r8_uniform_01(ref seed);
                s = UniformRNG.r8_uniform_01(ref seed);
            }
            else if (code == "T3" || code == "T4" || code == "T6" || code == "T10")
            {
                r = UniformRNG.r8_uniform_01(ref seed);
                s = UniformRNG.r8_uniform_01(ref seed);

                if (1.0 < r + s)
                {
                    r = 1.0 - r;
                    s = 1.0 - s;
                }
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("REFERENCE_SAMPLE - Fatal error!");
                Console.WriteLine("  Illegal code = \"" + code + "\".");
            }
        }

        public static void reference_to_physical_q4(double[] q4, int n, double[] rs,
        ref double[] xy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    REFERENCE_TO_PHYSICAL_Q4 maps Q4 reference points to physical points.
        //
        //  Discussion:
        //
        //    XY(R,S) = XY(0,0) * (1-R) * (1-S)
        //            + XY(1,0) *    R  * (1-S)
        //            + XY(1,1) *    R  *    S
        //            + XY(0,1) * (1-R) *    S
        //
        //  Reference Element Q4:
        //
        //    |
        //    1  4-----3
        //    |  |     |
        //    |  |     |
        //    S  |     |
        //    |  |     |
        //    |  |     |
        //    0  1-----2
        //    |
        //    +--0--R--1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double Q4[2*4], the coordinates of the vertices.
        //    The vertices are assumed to be the images of the reference vertices
        //    (0,0), (1,0), (1,1) and (0,1) respectively.
        //
        //    Input, int N, the number of points to transform.
        //
        //    Input, double RS[2*N], (R,S) points in the reference element.
        //
        //    Output, double XY[2*N], (X,Y) points in the physical element.
        //
        {
            int j;
            double[] psi;

            psi = new double[4 * n];

            for (j = 0; j < n; j++)
            {
                psi[0 + j * 2] = (1.0 - rs[0 + j * 2]) * (1.0 - rs[1 + j * 2]);
                psi[1 + j * 2] = rs[0 + j * 2] * (1.0 - rs[1 + j * 2]);
                psi[2 + j * 2] = rs[0 + j * 2] * rs[1 + j * 2];
                psi[3 + j * 2] = (1.0 - rs[0 + j * 2]) * rs[1 + j * 2];
            }

            typeMethods.r8mat_mm(2, 4, n, q4, psi, ref xy);
        }

        public static void reference_to_physical_t3(double[] t, int n, double[] ref_, ref double[] phy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    REFERENCE_TO_PHYSICAL_T3 maps T3 reference points to physical points.
        //
        //  Discussion:
        //
        //    Given the vertices of an order 3 physical triangle and a point
        //    (XSI,ETA) in the reference triangle, the routine computes the value
        //    of the corresponding image point (X,Y) in physical space.
        //
        //    Note that this routine may also be appropriate for an order 6
        //    triangle, if the mapping between reference and physical space
        //    is linear.  This implies, in particular, that the sides of the
        //    image triangle are straight and that the "midside" nodes in the
        //    physical triangle are halfway along the sides of
        //    the physical triangle.
        //
        //  Reference Element T3:
        //
        //    |
        //    1  3
        //    |  ..
        //    |  . .
        //    S  .  .
        //    |  .   .
        //    |  .    .
        //    0  1-----2
        //    |
        //    +--0--R--1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double[] t, the coordinates of the vertices.
        //    The vertices are assumed to be the images of (0,0), (1,0) and
        //    (0,1) respectively.
        //
        //    Input, int N, the number of points to transform.
        //
        //    Input, double REF[2*N], points in the reference triangle.
        //
        //    Output, double PHY[2*N], corresponding points in the
        //    physical triangle.
        //
        {
            int i;
            int j;

            for (i = 0; i < 2; i++)
            {
                for (j = 0; j < n; j++)
                {
                    phy[i + j * 2] = t[i + 0 * 2] * (1.0 - ref_[
                    0 + j * 2] - ref_[
                    1 + j * 2] )
                    +t[i + 1 * 2] * + ref_[
                    0 + j * 2]
                    +t[i + 2 * 2] * + ref_[
                    1 + j * 2];
                }
            }

            return;
        }

        public static void reference_to_physical_t6(double[] t, int n, double[] ref_, double[] phy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    REFERENCE_TO_PHYSICAL_T6 maps T6 reference points to physical points.
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
        //    |  ..
        //    |  . .
        //    S  6  5
        //    |  .   .
        //    |  .    .
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
        //    Input, int N, the number of points to transform.
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
}