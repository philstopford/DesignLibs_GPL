﻿using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] r8vec_chebyshev_new ( int n, double a_first, double a_last )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_CHEBYSHEV_NEW creates a vector of Chebyshev spaced values.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double A_FIRST, A_LAST, the first and last entries.
        //
        //    Output, double R8VEC_CHEBYSHEV_NEW[N], a vector of Chebyshev spaced data.
        //
    {
        double[] a = new double[n];

        switch (n)
        {
            case 1:
                a[0] = ( a_first + a_last ) / 2.0;
                break;
            default:
            {
                int i;
                for ( i = 0; i < n; i++ )
                {
                    double theta = (n - i - 1) * Math.PI / (n - 1);

                    double c = Math.Cos ( theta );

                    switch (n % 2)
                    {
                        case 1:
                        {
                            if ( 2 * i + 1 == n )
                            {
                                c = 0.0;
                            }

                            break;
                        }
                    }

                    a[i] = ( ( 1.0 - c ) * a_first  
                             + ( 1.0 + c ) * a_last ) 
                           /   2.0;

                }

                break;
            }
        }

        return a;
    }
    public static double[] r8vec_cheby_extreme_new(int n, double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_CHEBY_EXTERME_NEW creates Chebyshev Extreme values in [A,B].
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double A, B, the interval.
        //
        //    Output, double R8VEC_CHEBY_EXTREME_NEW[N], a vector of Chebyshev spaced data.
        //
    {
        double[] x = new double[n];

        switch (n)
        {
            case 1:
                x[0] = (a + b) / 2.0;
                break;
            default:
            {
                int i;
                for (i = 0; i < n; i++)
                {
                    double theta = (n - i - 1) * Math.PI / (n - 1);

                    double c = Math.Cos(theta);

                    switch (n % 2)
                    {
                        case 1:
                        {
                            if (2 * i + 1 == n)
                            {
                                c = 0.0;
                            }

                            break;
                        }
                    }

                    x[i] = ((1.0 - c) * a
                            + (1.0 + c) * b)
                           / 2.0;
                }

                break;
            }
        }

        return x;
    }

    public static double[] r8vec_cheby_zero_new(int n, double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_CHEBY_ZERO_NEW creates Chebyshev Zero values in [A,B].
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double A, B, the interval.
        //
        //    Output, double R8VEC_CHEBY_ZERO_NEW[N], a vector of Chebyshev spaced data.
        //
    {
        double[] x = new double[n];

        switch (n)
        {
            case 1:
                x[0] = (a + b) / 2.0;
                break;
            default:
            {
                int i;
                for (i = 0; i < n; i++)
                {
                    double theta = (2 * (n - i) - 1) * Math.PI / (2 * n);

                    double c = Math.Cos(theta);

                    switch (n % 2)
                    {
                        case 1:
                        {
                            if (2 * i + 1 == n)
                            {
                                c = 0.0;
                            }

                            break;
                        }
                    }

                    x[i] = ((1.0 - c) * a
                            + (1.0 + c) * b)
                           / 2.0;
                }

                break;
            }
        }

        return x;
    }

    public static double[] r8vec_cheby1space_new(int n, double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_CHEBY1SPACE_NEW creates Type 1 Chebyshev values in [A,B].
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double A, B, the interval.
        //
        //    Output, double R8VEC_CHEBY1SPACE_NEW[N], a vector of Type 1
        //    Chebyshev spaced data.
        //
    {
        double[] x = new double[n];

        switch (n)
        {
            case 1:
                x[0] = (a + b) / 2.0;
                break;
            default:
            {
                for (int i = 0; i < n; i++)
                {
                    double theta = (n - i - 1) * Math.PI / (n - 1);

                    double c = Math.Cos(theta);

                    switch (n % 2)
                    {
                        case 1:
                        {
                            if (2 * i + 1 == n)
                            {
                                c = 0.0;
                            }

                            break;
                        }
                    }

                    x[i] = ((1.0 - c) * a
                            + (1.0 + c) * b)
                           / 2.0;
                }

                break;
            }
        }

        return x;
    }

    public static double[] r8vec_cheby2space_new(int n, double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_CHEBY2SPACE_NEW creates Type 2 Chebyshev values in [A,B].
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 July 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double A, B, the interval.
        //
        //    Output, double R8VEC_CHEBY2SPACE_NEW[N], a vector of Type 2
        //    Chebyshev spaced data.
        //
    {
        int i;

        double[] x = new double[n];

        for (i = 0; i < n; i++)
        {
            double theta = (n - i) * Math.PI / (n + 1);

            double c = Math.Cos(theta);

            x[i] = ((1.0 - c) * a
                    + (1.0 + c) * b)
                   / 2.0;
        }

        return x;
    }
        
}