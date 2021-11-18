using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void i4_to_pascal(int k, ref int i, ref int j )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_TO_PASCAL converts a linear index to Pascal triangle coordinates.
        //
        //  Discussion:
        //
        //    We describe the grid points in Pascal's triangle in two ways:
        //
        //    As a linear index K:
        //
        //                     1
        //                   2   3
        //                 4   5   6
        //               7   8   9   10
        //
        //    As elements (I,J) of Pascal's triangle:
        //
        //                     0,0
        //                  1,0   0,1
        //               2,0   1,1    0,2
        //            3,0   2,1   1,2    0,3
        //
        //  Example:
        //
        //     K  I  J
        //
        //     1  0  0
        //     2  1  0
        //     3  0  1
        //     4  2  0
        //     5  1  1
        //     6  0  2
        //     7  3  0
        //     8  2  1
        //     9  1  2
        //    10  0  3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int K, the linear index of the (I,J) element.
        //    1 <= K.
        //
        //    Output, int &I, &J, the Pascal indices.
        //
    {
        switch (k)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("I4_TO_PASCAL - Fatal error!");
                Console.WriteLine("  K must be positive.");
                return;
        }

        int d = i4_to_pascal_degree(k);

        j = k - d * (d + 1) / 2 - 1;
        i = d - j;
    }

    public static int i4_to_pascal_degree(int k)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_TO_PASCAL_DEGREE converts a linear index to a Pascal triangle degree.
        //
        //  Discussion:
        //
        //    We describe the grid points in Pascal's triangle in two ways:
        //
        //    As a linear index K:
        //
        //                     1
        //                   2   3
        //                 4   5   6
        //               7   8   9   10
        //
        //    As elements (I,J) of Pascal's triangle:
        //
        //                     0,0
        //                  1,0   0,1
        //               2,0   1,1    0,2
        //            3,0   2,1   1,2    0,3
        //
        //    The quantity D represents the "degree" of the corresponding monomial,
        //    that is, D = I + J.
        //
        //    We can compute D directly from K using the quadratic formula.
        //
        //  Example:
        //
        //     K  I  J  D
        //
        //     1  0  0  0
        //
        //     2  1  0  1
        //     3  0  1  1
        //
        //     4  2  0  2
        //     5  1  1  2
        //     6  0  2  2
        //
        //     7  3  0  3
        //     8  2  1  3
        //     9  1  2  3
        //    10  0  3  3
        //
        //    11  4  0  4
        //    12  3  1  4
        //    13  2  2  4
        //    14  1  3  4
        //    15  0  4  4
        //
        //    16  5  0  5
        //    17  4  1  5
        //    18  3  2  5
        //    19  2  3  5
        //    20  1  4  5
        //    21  0  5  5
        //
        //    22  6  0  6
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int K, the linear index of the (I,J) element.
        //    1 <= K.
        //
        //     Output, int I4_TO_PASCAL_DEGREE, the degree (sum) of the corresponding 
        //     Pascal indices.
        //
    {
        switch (k)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("I4_TO_PASCAL_DEGREE - Fatal error!");
                Console.WriteLine("  K must be positive.");
                return 1;
        }

        double arg = 1 + 8 * (k - 1);

        int d = (int) (0.5 * (-1.0 + Math.Sqrt(arg)));

        return d;
    }

    public static int pascal_to_i4(int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PASCAL_TO_I4 converts Pacal triangle coordinates to a linear index.
        //
        //  Discussion:
        //
        //    We describe the grid points in a Pascal triangle in two ways:
        //
        //    As a linear index K:
        //
        //                     1
        //                   2   3
        //                 4   5   6
        //               7   8   9   10
        //
        //    As elements (I,J) of Pascal's triangle:
        //
        //                     0,0
        //                  1,0   0,1
        //               2,0   1,1    0,2
        //            3,0   2,1   1,2    0,3
        //
        //  Example:
        //
        //     K  I  J
        //
        //     1  0  0
        //     2  1  0
        //     3  0  1
        //     4  2  0
        //     5  1  1
        //     6  0  2
        //     7  3  0
        //     8  2  1
        //     9  1  2
        //    10  0  3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, J, the row and column indices.  I and J 
        //    must be nonnegative.
        //
        //    Output, int PASCAL_TO_I4, the linear index of the (I,J) element.
        //
    {
        switch (i)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("PASCAL_TO_I4 - Fatal error!");
                Console.WriteLine("  I < 0.");
                Console.WriteLine("  I = " + i + "");
                return 1;
        }

        switch (j)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("PASCAL_TO_I4 - Fatal error!");
                Console.WriteLine("  J < 0.");
                Console.WriteLine("  J = " + j + "");
                return 1;
        }

        int d = i + j;

        int k = d * (d + 1) / 2 + j + 1;

        return k;
    }
}