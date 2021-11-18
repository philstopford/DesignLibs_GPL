using System;
using Burkardt.Types;

namespace Burkardt.Transform;

public static class Walsh
{
    public static void ffwt(int n, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FFWT performs an in-place fast Walsh transform.
        //
        //  Discussion:
        //
        //    This routine performs a fast Walsh transform on an input series X
        //    leaving the transformed results in X. 
        //    X is dimensioned N, which must be a power of 2.
        //    The results of this Walsh transform are in sequency order.
        //
        //    The output sequence could be normalized by dividing by N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 March 2011
        //
        //  Author:
        //
        //    Ken Beauchamp
        //
        //  Reference:
        //
        //    Ken Beauchamp,
        //    Walsh functions and their applications,
        //    Academic Press, 1975,
        //    ISBN: 0-12-084050-2,
        //    LC: QA404.5.B33.
        //
        //  Parameters:
        //
        //    Input, int N, the number of items in X.
        //    N must be a power of 2.
        //
        //    Input/output, double X[N], the data to be transformed.
        //
    {
        double hold;
        int i;
        int ii = 0;
        int k;
        int l;
        int nz2 = 0;

        int m = (int) Math.Log2(n);

        int[] two_power = new int[m];

        for (i = 0; i < m; i++)
        {
            two_power[i] = (int) Math.Pow(2, m - 1 - i);
        }

        for (l = 0; l < m; l++)
        {
            int nz = (int) Math.Pow(2, l);
            int nzi = 2 * nz;
            int nzn = n / nzi;
            nz2 = nz2 switch
            {
                0 => 1,
                _ => nz / 2
            };

            for (i = 0; i < nzn; i++)
            {
                int js = i * nzi;
                double z = 1.0;
                for (ii = 0; ii < 2; ii++)
                {
                    int j;
                    for (j = 0; j < nz2; j++)
                    {
                        js += 1;
                        int j2 = js + nz;
                        hold = x[js - 1] + z * x[j2 - 1];
                        z = -z;
                        x[j2 - 1] = x[js - 1] + z * x[j2 - 1];
                        x[js - 1] = hold;
                        z = -z;
                    }

                    if (l == 0)
                    {
                        break;
                    }

                    z = -1.0;
                }
            }
        }

        //
        //  Bit reversal section.
        //
        int nw = 0;
        for (k = 0; k < n; k++)
        {
            //
            //  Choose correct index and switch elements if not already switched.
            //
            if (k < nw)
            {
                hold = x[nw];
                x[nw] = x[k];
                x[k] = hold;
            }

            //
            //  Bump up series by 1.
            //
            for (i = 0; i < m; i++)
            {
                ii = i;
                if (nw < two_power[i])
                {
                    break;
                }

                int mw = nw / two_power[i];
                int mw1 = mw / 2;
                if (mw <= 2 * mw1)
                {
                    break;
                }

                nw -= two_power[i];
            }

            nw += two_power[ii];
        }

    }

    public static void fwt(int n, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FWT performs a fast Walsh transform.
        //
        //  Discussion:
        //
        //    This routine performs a fast Walsh transform on an input series X
        //    leaving the transformed results in X. 
        //    X is dimensioned N, which must be a power of 2.
        //    The results of this Walsh transform are in sequency order.
        //
        //    The output sequence could be normalized by dividing by N.
        //
        //    Note that the program text in the reference included the line
        //      y(jd) = abs ( x(j) - x(j2) )
        //    which has been corrected to:
        //      y(jd) = x(j) - x(j2)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 March 2011
        //
        //  Author:
        //
        //    Ken Beauchamp
        //
        //  Reference:
        //
        //    Ken Beauchamp,
        //    Walsh functions and their applications,
        //    Academic Press, 1975,
        //    ISBN: 0-12-084050-2,
        //    LC: QA404.5.B33.
        //
        //  Parameters:
        //
        //    Input, int N, the number of items in X.
        //    N must be a power of 2.
        //
        //    Input/output, double X[N], the data to be transformed.
        //
    {
        int l;

        double[] y = new double[n];

        int n2 = n / 2;
        int m = (int) Math.Log2(n);

        for (l = 1; l <= m; l++)
        {
            int ny = 0;
            int nz = (int) Math.Pow(2, l - 1);
            int nzi = 2 * nz;
            int nzn = n / nzi;
            int i;
            for (i = 1; i <= nzn; i++)
            {
                int nx = ny + 1;
                ny += nz;
                int js = (i - 1) * nzi;
                int jd = js + nzi + 1;
                int j;
                for (j = nx; j <= ny; j++)
                {
                    js += 1;
                    int j2 = j + n2;
                    y[js - 1] = x[j - 1] + x[j2 - 1];
                    jd -= 1;
                    y[jd - 1] = x[j - 1] - x[j2 - 1];
                }
            }

            typeMethods.r8vec_copy(n, y, ref x);
        }

    }
        
    public static void walsh ( int n, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WALSH performs a fast Walsh transform.
        //
        //  Discussion:
        //
        //    This routine performs a fast Wash transform on an input series X
        //    leaving the transformed results in X.  The array Y is working space.
        //    X and Y are dimensioned N, which must be a power of 2.
        //    The results of this Walsh transform are in sequency order.
        //
        //    The output sequence could be normalized by dividing by N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 March 2011
        //
        //  Author:
        //
        //    Ken Beauchamp
        //
        //  Reference:
        //
        //    Ken Beauchamp,
        //    Walsh functions and their applications,
        //    Academic Press, 1975,
        //    ISBN: 0-12-084050-2,
        //    LC: QA404.5.B33.
        //
        //  Parameters:
        //
        //    Input, int N, the number of items in X.
        //    N must be a power of 2.
        //
        //    Input/output, double X[N], the data to be transformed.
        //
    {
        int j;

        int n2 = n / 2;
        double[] y = new double[n2];
        int m = (int)Math.Log2 ( n );
        const double z = - 1.0;

        for ( j = 1; j <= m; j++ )
        {
            int n1 = (int)Math.Pow ( 2, m - j + 1 );
            int j1 = (int)Math.Pow ( 2, j - 1 );
            int l;
            for ( l = 1; l <= j1; l++ )
            {
                int is_ = ( l - 1 ) * n1 + 1;
                int i1 = 0;
                double w = z;
                int i;
                for ( i = is_; i <= is_ + n1 - 1; i += 2 )
                {
                    double a = x[i-1];
                    x[is_+i1-1] = a + x[i];
                    i1 += 1;
                    y[i1-1] = ( x[i] - a ) * w;
                    w *= z;
                }
                for ( i = 1; i <= n1 / 2; i++ )
                {
                    x[n1/2+is_+i-2] = y[i-1];
                }
            }
        }
    }
}