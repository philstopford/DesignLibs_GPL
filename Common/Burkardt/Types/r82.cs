using System;
using Burkardt.Uniform;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double r82_dist_l2(double[] a1, double[] a2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R82_DIST_L2 returns the L2 distance between a pair of R82's.
            //
            //  Discussion:
            //
            //    An R82 is a vector of type R8, with two entries.
            //
            //    The vector L2 norm is defined as:
            //
            //      sqrt ( sum ( 1 <= I <= N ) A(I) * A(I) ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A1[2], A2[2], the vectors.
            //
            //    Output, double R82_DIST_L2, the L2 norm of A1 - A2.
            //
        {
            double value;

            value = Math.Sqrt(Math.Pow(a1[0] - a2[0], 2)
                         + Math.Pow(a1[1] - a2[1], 2));

            return value;
        }

        public static void r82_print(double[] a, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R82_PRINT prints an R82.
            //
            //  Discussion:
            //
            //    An R82 is an R8VEC with two entries.
            //
            //    A format is used which suggests a coordinate pair:
            //
            //  Example:
            //
            //    Center : ( 1.23, 7.45 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 July 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A[2], the coordinates of the vector.
            //
            //    Input, string TITLE, a title.
            //
        {
            Console.WriteLine("  " + title + " : " + ": ( " + a[0].ToString().PadLeft(12)
                + ", " + a[1].ToString().PadLeft(12) + " )");
        }

        public static void r82_uniform_ab(double b, double c, ref int seed, ref double[] r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R82_UNIFORM_AB returns a random R82 value in a given range.
            //
            //  Discussion:
            //
            //    An R82 is an R8VEC with two entries.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double B, C, the minimum and maximum values.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double R[2], the randomly chosen value.
            //
        {
            int i;

            for (i = 0; i < 2; i++)
            {
                r[i] = UniformRNG.r8_uniform_ab(b, c, ref seed);
            }
        }

    }
}