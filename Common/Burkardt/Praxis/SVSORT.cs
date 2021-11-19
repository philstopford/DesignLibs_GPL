namespace Burkardt.Praxis;

public static class SVSORT
{
    public static void svsort ( int n, ref double[] d, ref double[] v ) 

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SVSORT descending sorts D and adjusts the corresponding columns of V.
        //
        //  Discussion:
        //
        //    A simple bubble sort is used on D.
        //
        //    In our application, D contains singular values, and the columns of V are
        //    the corresponding right singular vectors.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Richard Brent.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Richard Brent,
        //    Algorithms for Minimization with Derivatives,
        //    Prentice Hall, 1973,
        //    Reprinted by Dover, 2002.
        //
        //  Parameters:
        //
        //    Input, int N, the length of D, and the order of V.
        //
        //    Input/output, double D[N], the vector to be sorted.  
        //    On output, the entries of D are in descending order.
        //
        //    Input/output, double V[N,N], an N by N array to be adjusted 
        //    as D is sorted.  In particular, if the value that was in D(I) on input is
        //    moved to D(J) on output, then the input column V(*,I) is moved to
        //    the output column V(*,J).
        //
    {
        int j1;

        for ( j1 = 0; j1 < n - 1; j1++ )
        {
            //
            //  Find J3, the index of the largest entry in D[J1:N-1].
            //  MAXLOC apparently requires its output to be an array.
            //
            int j3 = j1;
            int j2;
            for ( j2 = j1 + 1; j2 < n; j2++ )
            {
                if ( d[j3] < d[j2] )
                {
                    j3 = j2;
                }
            }
            //
            //  If J1 != J3, swap D[J1] and D[J3], and columns J1 and J3 of V.
            //
            if (j1 == j3)
            {
                continue;
            }

            double t = d[j1];
            d[j1] = d[j3];
            d[j3] = t;
            int i;
            for ( i = 0; i < n; i++ )
            {
                t         = v[i+j1*n];
                v[i+j1*n] = v[i+j3*n];
                v[i+j3*n] = t;
            }
        }
    }
}