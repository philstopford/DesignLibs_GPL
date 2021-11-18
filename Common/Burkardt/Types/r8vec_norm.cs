﻿using System;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8vec_norm_l0(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_NORM_L0 returns the l0 "norm" of an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The l0 "norm" simply counts the number of nonzero entries in the vector.
        //    It is not a true norm, but has some similarities to one.  It is useful
        //    in the study of compressive sensing.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double A(N), the vector.
        //
        //    Output, double R8VEC_NORM_L0, the value of the norm.
        //
    {
        int i;

        double value = 0.0;
        for (i = 0; i < n; i++)
        {
            if (a[i] != 0.0)
            {
                value += 1.0;
            }
        }

        return value;
    }

    public static double r8vec_norm_l1(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_NORM_L1 returns the L1 norm of an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The vector L1 norm is defined as:
        //
        //      R8VEC_NORM_L1 = sum ( 1 <= I <= N ) abs ( A(I) ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 March 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, double A[N], the vector whose L1 norm is desired.
        //
        //    Output, double R8VEC_NORM_L1, the L1 norm of A.
        //
    {
        int i;

        double v = 0.0;

        for (i = 0; i < n; i++)
        {
            v += Math.Abs(a[i]);
        }

        return v;
    }

    public static double r8vec_norm_l2(int n, double[] a, int index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_NORM_L2 returns the L2 norm of an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The vector L2 norm is defined as:
        //
        //      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 March 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, double A[N], the vector whose L2 norm is desired.
        //
        //    Output, double R8VEC_NORM_L2, the L2 norm of A.
        //
    {
        int i;

        double v = 0.0;

        for (i = 0; i < n; i++)
        {
            v += a[i + index] * a[i + index];
        }

        v = Math.Sqrt(v);

        return v;
    }

    public static double r8vec_norm_li(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_NORM_LI returns the L-oo norm of an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The vector L-oo norm is defined as:
        //
        //      R8VEC_NORM_LI = max ( 1 <= I <= N ) abs ( A(I) ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 March 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, double A[N], the vector whose L-oo norm is desired.
        //
        //    Output, double R8VEC_NORM_LI, the L-oo norm of A.
        //
    {
        int i;

        double v1 = 0.0;

        for (i = 0; i < n; i++)
        {
            double v2 = Math.Abs(a[i]);

            if (v1 < v2)
            {
                v1 = v2;
            }
        }

        return v1;
    }

    public static double r8vec_norm_lp(int n, double[] a, double p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_NORM_LP returns the LP norm of an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The vector LP norm is defined as:
        //
        //      R8VEC_NORM_LP = ( sum ( 1 <= I <= N ) ( abs ( A(I) ) )^P )^(1/P).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 March 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, double A[N], the vector whose LP norm is desired.
        //
        //    Input, double P, the index of the norm.
        //
        //    Output, double R8VEC_NORML_LP, the LP norm of A.
        //
    {
        int i;

        double v = 0.0;

        switch (p)
        {
            case 1.0:
            {
                for (i = 0; i < n; i++)
                {
                    v += Math.Abs(a[i]);
                }

                break;
            }
            case 2.0:
            {
                for (i = 0; i < n; i++)
                {
                    v += a[i] * a[i];
                }

                v = Math.Sqrt(v);
                break;
            }
            default:
            {
                for (i = 0; i < n; i++)
                {
                    v += Math.Pow(Math.Abs(a[i]), p);
                }

                v = Math.Pow(v, 1.0 / p);
                break;
            }
        }

        return v;
    }

    public static double r8vec_norm_rms(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_NORM_RMS returns the RMS norm of an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The vector RMS norm is defined as:
        //
        //      R8VEC_RMS_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)^2 / N ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, double A[N], the vector.
        //
        //    Output, double R8VEC_NORM_RMS, the RMS norm of A.
        //
    {
        double v = 0.0;

        switch (n)
        {
            case > 0:
            {
                int i;
                for (i = 0; i < n; i++)
                {
                    v += a[i] * a[i];
                }

                v = Math.Sqrt(v / n);
                break;
            }
        }

        return v;
    }

    public static double r8vec_norm ( int n, double[] a, int index = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_NORM returns the L2 norm of an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The vector L2 norm is defined as:
        //
        //      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 March 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, double A[N], the vector whose L2 norm is desired.
        //
        //    Output, double R8VEC_NORM, the L2 norm of A.
        //
    {
        int i;

        double v = 0.0;

        for ( i = 0; i < n; i++ )
        {
            v += a[index + i] * a[index + i];
        }
        v = Math.Sqrt ( v );

        return v;
    }
        
    public static double r8vec_norm_affine(int n, double[] v0, double[] v1, int v0Index = 0, int v1Index = 0)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_NORM_AFFINE returns the affine L2 norm of an R8VEC.
        //
        //  Discussion:
        //
        //    The affine vector L2 norm is defined as:
        //
        //      R8VEC_NORM_AFFINE(V0,V1)
        //        = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the vectors.
        //
        //    Input, double V0[N], the base vector.
        //
        //    Input, double V1[N], the vector.
        //
        //    Output, double R8VEC_NORM_AFFINE, the affine L2 norm.
        //
    {
        double value = 0.0;

        for (int i = 0; i < n; i++)
        {
            value += (v1[(i + v1Index ) % v1.Length] - v0[(i + v0Index ) % v0.Length]) * (v1[(i + v1Index ) % v1.Length] - v0[(i + v0Index ) % v0.Length]);
        }

        value = Math.Sqrt(value);

        return value;
    }

    public static void r8vec_normal_01(int n, ref int seed, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The standard normal probability distribution function (PDF) has
        //    mean 0 and standard deviation 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of values desired.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double X[N], a sample of the standard normal PDF.
        //
        //  Local parameters:
        //
        //    Local, double R[N+1], is used to store some uniform random values.
        //    Its dimension is N+1, but really it is only needed to be the
        //    smallest even number greater than or equal to N.
        //
        //    Local, int X_LO, X_HI, records the range of entries of
        //    X that we need to compute.
        //
    {
        double[] r;
        //
        //  Record the range of X we need to fill in.
        //
        const int x_lo = 1;
        int x_hi = n;
        switch (x_hi - x_lo + 1)
        {
            //
            //  If we need just one new value, do that here to avoid null arrays.
            //
            case 1:
                r = UniformRNG.r8vec_uniform_01_new(2, ref seed);

                x[x_hi - 1] = Math.Sqrt(-2.0 * Math.Log(r[0])) * Math.Cos(2.0 * Math.PI * r[1]);
                break;
            //
            default:
            {
                int i;
                int m;
                switch ((x_hi - x_lo + 1) % 2)
                {
                    case 0:
                    {
                        m = (x_hi - x_lo + 1) / 2;

                        r = UniformRNG.r8vec_uniform_01_new(2 * m, ref seed);

                        for (i = 0; i <= 2 * m - 2; i += 2)
                        {
                            x[x_lo + i - 1] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Cos(2.0 * Math.PI * r[i + 1]);
                            x[x_lo + i] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Sin(2.0 * Math.PI * r[i + 1]);
                        }

                        break;
                    }
                    //
                    default:
                    {
                        x_hi -= 1;

                        m = (x_hi - x_lo + 1) / 2 + 1;

                        r = UniformRNG.r8vec_uniform_01_new(2 * m, ref seed);

                        for (i = 0; i <= 2 * m - 4; i += 2)
                        {
                            x[x_lo + i - 1] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Cos(2.0 * Math.PI * r[i + 1]);
                            x[x_lo + i] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Sin(2.0 * Math.PI * r[i + 1]);
                        }

                        i = 2 * m - 2;

                        x[x_lo + i - 1] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Cos(2.0 * Math.PI * r[i + 1]);
                        break;
                    }
                }

                break;
            }
        }
    }

    public class r8vecNormalData
    {
        public int made;
        public int saved;
        public double y;
    }
        
    public static double[] r8vec_normal_01_new(int n, ref r8vecNormalData data, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_NORMAL_01_NEW returns a unit pseudonormal R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The standard normal probability distribution function (PDF) has
        //    mean 0 and standard deviation 1.
        //
        //    This routine can generate a vector of values on one call.  It
        //    has the feature that it should provide the same results
        //    in the same order no matter how we break up the task.
        //
        //    Before calling this routine, the user may call RANDOM_SEED
        //    in order to set the seed of the random number generator.
        //
        //    The Box-Muller method is used, which is efficient, but
        //    generates an even number of values each time.  On any call
        //    to this routine, an even number of new values are generated.
        //    Depending on the situation, one value may be left over.
        //    In that case, it is saved for the next call.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 February 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of values desired.  If N is negative,
        //    then the code will flush its internal memory; in particular,
        //    if there is a saved value to be used on the next call, it is
        //    instead discarded.  This is useful if the user has reset the
        //    random number seed, for instance.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8VEC_NORMAL_01_NEW[N], a sample of the standard normal PDF.
        //
        //  Local parameters:
        //
        //    Local, int MADE, records the number of values that have
        //    been computed.  On input with negative N, this value overwrites
        //    the return value of N, so the user can get an accounting of
        //    how much work has been done.
        //
        //    Local, double R[N+1], is used to store some uniform random values.
        //    Its dimension is N+1, but really it is only needed to be the
        //    smallest even number greater than or equal to N.
        //
        //    Local, int SAVED, is 0 or 1 depending on whether there is a
        //    single saved value left over from the previous call.
        //
        //    Local, int X_LO, X_HI, records the range of entries of
        //    X that we need to compute.  This starts off as 1:N, but is adjusted
        //    if we have a saved value that can be immediately stored in X(1),
        //    and so on.
        //
        //    Local, double Y, the value saved from the previous call, if
        //    SAVED is 1.
        //
    {
        double[] r;
        switch (n)
        {
            //
            //  I'd like to allow the user to reset the internal data.
            //  But this won't work properly if we have a saved value Y.
            //  I'm making a crock option that allows the user to signal
            //  explicitly that any internal memory should be flushed,
            //  by passing in a negative value for N.
            //
            case < 0:
                data.made = 0;
                data.saved = 0;
                data.y = 0.0;
                return null;
            case 0:
                return null;
        }

        double[] x = new double[n];
        //
        //  Record the range of X we need to fill in.
        //
        int x_lo = 1;
        int x_hi = n;
        switch (data.saved)
        {
            //
            //  Use up the old value, if we have it.
            //
            case 1:
                x[0] = data.y;
                data.saved = 0;
                x_lo = 2;
                break;
        }

        switch (x_hi - x_lo + 1)
        {
            //
            //  Maybe we don't need any more values.
            //
            case 0:
                break;
            //
            //  If we need just one new value, do that here to avoid null arrays.
            //
            case 1:
                r = UniformRNG.r8vec_uniform_01_new(2, ref seed);

                x[x_hi - 1] = Math.Sqrt(-2.0 * Math.Log(r[0])) * Math.Cos(2.0 * Math.PI * r[1]);
                data.y = Math.Sqrt(-2.0 * Math.Log(r[0])) * Math.Sin(2.0 * Math.PI * r[1]);

                data.saved = 1;

                data.made += 2;
                break;
            //
            default:
            {
                int i;
                int m;
                switch ((x_hi - x_lo + 1) % 2)
                {
                    case 0:
                    {
                        m = (x_hi - x_lo + 1) / 2;

                        r = UniformRNG.r8vec_uniform_01_new(2 * m, ref seed);

                        for (i = 0; i <= 2 * m - 2; i += 2)
                        {
                            x[x_lo + i - 1] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Cos(2.0 * Math.PI * r[i + 1]);
                            x[x_lo + i] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Sin(2.0 * Math.PI * r[i + 1]);
                        }

                        data.made = data.made + x_hi - x_lo + 1;
                        break;
                    }
                    //
                    default:
                    {
                        x_hi -= 1;

                        m = (x_hi - x_lo + 1) / 2 + 1;

                        r = UniformRNG.r8vec_uniform_01_new(2 * m, ref seed);

                        for (i = 0; i <= 2 * m - 4; i += 2)
                        {
                            x[x_lo + i - 1] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Cos(2.0 * Math.PI * r[i + 1]);
                            x[x_lo + i] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Sin(2.0 * Math.PI * r[i + 1]);
                        }

                        i = 2 * m - 2;

                        x[x_lo + i - 1] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Cos(2.0 * Math.PI * r[i + 1]);
                        data.y = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Sin(2.0 * Math.PI * r[i + 1]);

                        data.saved = 1;

                        data.made = data.made + x_hi - x_lo + 2;
                        break;
                    }
                }

                break;
            }
        }

        return x;
    }

    public static double[] r8vec_any_normal(int dim_num, double[] v1)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_ANY_NORMAL returns some normal vector to V1.
        //
        //  Discussion:
        //
        //    If DIM_NUM < 2, then no normal vector can be returned.
        //
        //    If V1 is the zero vector, then any unit vector will do.
        //
        //    No doubt, there are better, more robust algorithms.  But I will take
        //    just about ANY reasonable unit vector that is normal to V1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double V1[DIM_NUM], the vector.
        //
        //    Output, double R8VEC_ANY_NORMAL[DIM_NUM], a vector that is
        //    normal to V2, and has unit Euclidean length.
        //
    {
        int i;

        switch (dim_num)
        {
            case < 2:
                Console.WriteLine("");
                Console.WriteLine("R8VEC_ANY_NORMAL - Fatal error!");
                Console.WriteLine("  Called with DIM_NUM < 2.");
                return null;
        }

        double[] v2 = new double[dim_num];

        if (r8vec_norm(dim_num, v1) == 0.0)
        {
            r8vec_zeros(dim_num, ref v2);
            v2[0] = 1.0;
            return v2;
        }

        //
        //  Seek the largest entry in V1, VJ = V1(J), and the
        //  second largest, VK = V1(K).
        //
        //  Since V1 does not have zero norm, we are guaranteed that
        //  VJ, at least, is not zero.
        //
        int j = -1;
        double vj = 0.0;

        int k = -1;
        double vk = 0.0;

        for (i = 0; i < dim_num; i++)
        {
            if (!(Math.Abs(vk) < Math.Abs(v1[i])) && k != -1)
            {
                continue;
            }

            if (Math.Abs(vj) < Math.Abs(v1[i]) || j == -1)
            {
                k = j;
                vk = vj;
                j = i;
                vj = v1[i];
            }
            else
            {
                k = i;
                vk = v1[i];
            }
        }

        //
        //  Setting V2 to zero, except that V2(J) = -VK, and V2(K) = VJ,
        //  will just about do the trick.
        //
        r8vec_zeros(dim_num, ref v2);

        v2[j] = -vk / Math.Sqrt(vk * vk + vj * vj);
        v2[k] = vj / Math.Sqrt(vk * vk + vj * vj);

        return v2;
    }

    public static double[] r8vec_normal_ab_new(int n, double b, double c, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_NORMAL_AB_NEW returns a scaled pseudonormal R8VEC.
        //
        //  Discussion:
        //
        //    The scaled normal probability distribution function (PDF) has
        //    mean A and standard deviation B.
        //
        //    This routine can generate a vector of values on one call.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of values desired.
        //
        //    Input, double B, C, the mean and standard deviation.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8VEC_NORMAL_AB_NEW[N], a sample of the standard normal PDF.
        //
        //  Local parameters:
        //
        //    Local, double R(N+1), is used to store some uniform random values.
        //    Its dimension is N+1, but really it is only needed to be the
        //    smallest even number greater than or equal to N.
        //
        //    Local, int X_LO, X_HI, records the range of entries of
        //    X that we need to compute.
        //
    {
        int i;
        double[] r;

        double[] x = new double[n];
        //
        //  Record the range of X we need to fill in.
        //
        const int x_lo = 1;
        int x_hi = n;
        switch (x_hi - x_lo + 1)
        {
            //
            //  If we need just one new value, do that here to avoid null arrays.
            //
            case 1:
                r = UniformRNG.r8vec_uniform_01_new(2, ref seed);

                x[x_hi - 1] = Math.Sqrt(-2.0 * Math.Log(r[0])) * Math.Cos(2.0 * Math.PI * r[1]);
                break;
            //
            default:
            {
                int m;
                switch ((x_hi - x_lo + 1) % 2)
                {
                    case 0:
                    {
                        m = (x_hi - x_lo + 1) / 2;

                        r = UniformRNG.r8vec_uniform_01_new(2 * m, ref seed);

                        for (i = 0; i <= 2 * m - 2; i += 2)
                        {
                            x[x_lo + i - 1] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Cos(2.0 * Math.PI * r[i + 1]);
                            x[x_lo + i] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Sin(2.0 * Math.PI * r[i + 1]);
                        }

                        break;
                    }
                    //
                    default:
                    {
                        x_hi -= 1;

                        m = (x_hi - x_lo + 1) / 2 + 1;

                        r = UniformRNG.r8vec_uniform_01_new(2 * m, ref seed);

                        for (i = 0; i <= 2 * m - 4; i += 2)
                        {
                            x[x_lo + i - 1] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Cos(2.0 * Math.PI * r[i + 1]);
                            x[x_lo + i] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Sin(2.0 * Math.PI * r[i + 1]);
                        }

                        i = 2 * m - 2;

                        x[x_lo + i - 1] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Cos(2.0 * Math.PI * r[i + 1]);
                        break;
                    }
                }

                break;
            }
        }

        for (i = 0; i < n; i++)
        {
            x[i] = b + c * x[i];
        }

        return x;
    }

    public static double r8vec_normsq ( int n, double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_NORMSQ returns the squared L2 norm of an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The squared vector L2 norm is defined as:
        //
        //      R8VEC_NORMSQ = sum ( 1 <= I <= N ) A(I)^2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the vector dimension.
        //
        //    Input, double A[N], the vector.
        //
        //    Output, double R8VEC_NORMSQ, the squared L2 norm.
        //
    {
        int i;

        double v = 0.0;

        for ( i = 0; i < n; i++ )
        {
            v += a[i] * a[i];
        }
        return v;
    }

    public static double r8vec_normsq_affine ( int n, double[] v0, double[] v1, int v0Index = 0, int v1Index = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_NORMSQ_AFFINE returns the squared affine L2 norm of an R8VEC.
        //
        //  Discussion:
        //
        //    The squared affine vector L2 norm is defined as:
        //
        //      R8VEC_NORMSQ_AFFINE(V0,V1)
        //        = sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the vectors.
        //
        //    Input, double V0[N], the base vector.
        //
        //    Input, double V1[N], the vector whose squared affine L2 norm is desired.
        //
        //    Output, double R8VEC_NORMSQ_AFFINE, the squared affine L2 norm.
        //
    {
        int i;

        double value = 0.0;

        for ( i = 0; i < n; i++ )
        {
            value += ( v1[(i + v1Index) % v1.Length] - v0[(i + v0Index) % v0.Length] ) * ( v1[(i + v1Index) % v1.Length] - v0[(i + v0Index) % v0.Length] );
        }
        return value;
    }
        
}