namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void l4vec_next ( int n, ref bool[] l4vec )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    L4VEC_NEXT generates the next logical vector.
        //
        //  Discussion:
        //
        //    In the following discussion, we will let '0' stand for FALSE and
        //    '1' for TRUE.
        //
        //    The vectors have the order
        //
        //      (0,0,...,0),
        //      (0,0,...,1), 
        //      ...
        //      (1,1,...,1)
        //
        //    and the "next" vector after (1,1,...,1) is (0,0,...,0).  That is,
        //    we allow wrap around.
        //
        //  Example:
        //
        //    N = 3
        //
        //    Input      Output
        //    -----      ------
        //    0 0 0  =>  0 0 1
        //    0 0 1  =>  0 1 0
        //    0 1 0  =>  0 1 1
        //    0 1 1  =>  1 0 0
        //    1 0 0  =>  1 0 1
        //    1 0 1  =>  1 1 0
        //    1 1 0  =>  1 1 1
        //    1 1 1  =>  0 0 0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 May 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the vectors.
        //
        //    Input/output, bool L4VEC[N], on output, the successor to the
        //    input vector.  
        //
    {
        int i;

        for ( i = n - 1; 0 <= i; i-- )
        {
            switch (l4vec[i])
            {
                case false:
                    l4vec[i] = true;
                    return;
                default:
                    l4vec[i] = false;
                    break;
            }
        }
    }
}