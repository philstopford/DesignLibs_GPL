namespace Burkardt.Function
{
    public static class Pythagorus
    {
        public static void pythag_triple_next ( ref int i, ref int j, ref int a, ref int b, ref int c )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYTHAG_TRIPLE_NEXT computes the next Pythagorean triple.
        //
        //  Example:
        //
        //     I       J       A       B       C    A^2+B^2 = C^2
        //
        //     2       1       3       4       5      25
        //     3       2       5      12      13     169
        //     4       1      15       8      17     289
        //     4       3       7      24      25     625
        //     5       2      21      20      29     841
        //     5       4       9      40      41    1681
        //     6       1      35      12      37    1369
        //     6       3      27      36      45    2025
        //     6       5      11      60      61    3721
        //     7       2      45      28      53    2809
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, ref int I, &J, the generators.
        //    On first call, set I = J = 0.  On repeated calls, leave I and J
        //    at their output values from the previous call.
        //
        //    Output, ref int A, &B, &C, the next Pythagorean triple.
        //    A, B, and C are positive integers which have no common factors,
        //    and A^2 + B^2 = C^2.
        //
        {
            //
            //  I starts at 2 and increases;
            //
            //  J starts out at 2 if I is odd, or 1 if I is even, increases by 2,
            //    but is always less than I.
            //
  
            if ( i == 0 && j == 0 )
            {
                i = 2;
                j = 1;
            }
            else if ( j + 2 < i )
            {
                j = j + 2;
            }
            else
            {
                i = i + 1;
                j = ( i % 2 ) + 1;
            }

            a = i * i - j * j;
            b = 2 * i * j;
            c = i * i + j * j;
        }
    }
}