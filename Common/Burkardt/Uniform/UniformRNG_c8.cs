using System;
using System.Numerics;

namespace Burkardt.Uniform;

public static partial class UniformRNG
{
    public static Complex c8_uniform_01(ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8_UNIFORM_01 returns a unit double complex pseudorandom number.
        //
        //  Discussion:
        //
        //    The angle should be uniformly distributed between 0 and 2 * PI,
        //    the square root of the radius uniformly distributed between 0 and 1.
        //
        //    This results in a uniform distribution of values in the unit circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 April 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, the "seed" value, which should NOT be 0.
        //    On output, SEED has been updated.
        //
        //    Output, complex <double> C8_UNIFORM_01, a pseudorandom complex value.
        //
    {
        int k = seed / 127773;

        seed = 16807 * (seed - k * 127773) - k * 2836;

        switch (seed)
        {
            case < 0:
                seed += 2147483647;
                break;
        }

        double r = Math.Sqrt(seed * 4.656612875E-10);

        k = seed / 127773;

        seed = 16807 * (seed - k * 127773) - k * 2836;

        switch (seed)
        {
            case < 0:
                seed += 2147483647;
                break;
        }

        double theta = 2.0 * Math.PI * (seed * 4.656612875E-10);

        Complex value = new(r * Math.Cos(theta), r * Math.Sin(theta));

        return value;
    }
        
    public static Complex[] c8vec_uniform_01_new ( int n, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8VEC_UNIFORM_01_NEW returns a unit pseudorandom C8VEC.
        //
        //  Discussion:
        //
        //    A C8VEC is a vector of complex <double> values.
        //
        //    The angles should be uniformly distributed between 0 and 2 * PI,
        //    the square roots of the radius uniformly distributed between 0 and 1.
        //
        //    This results in a uniform distribution of values in the unit circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of values to compute.
        //
        //    Input/output, int &SEED, the "seed" value, which should NOT be 0.
        //    On output, SEED has been updated.
        //
        //    Output, complex <double> C8VEC_UNIFORM_01_NEW[N], the pseudorandom 
        //    complex vector.
        //
    {
        int i;

        Complex[] c = new Complex [n];

        for ( i = 0; i < n; i++ )
        {
            int k = seed / 127773;

            seed = 16807 * ( seed - k * 127773 ) - k * 2836;

            switch (seed)
            {
                case < 0:
                    seed += 2147483647;
                    break;
            }

            double r = Math.Sqrt ( seed * 4.656612875E-10 );

            k = seed / 127773;

            seed = 16807 * ( seed - k * 127773 ) - k * 2836;

            switch (seed)
            {
                case < 0:
                    seed += 2147483647;
                    break;
            }

            double theta = 2.0 * Math.PI * ( seed * 4.656612875E-10 );

            c[i] = r * new Complex ( Math.Cos ( theta ), Math.Sin ( theta ) );
        }

        return c;
    }

    public static Complex[] c8mat_uniform_01_new(int m, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8MAT_UNIFORM_01_NEW returns a unit pseudorandom C8MAT.
        //
        //  Discussion:
        //
        //    A C8MAT is an array of complex <double> values.
        //
        //    The angles should be uniformly distributed between 0 and 2 * PI,
        //    the square roots of the radius uniformly distributed between 0 and 1.
        //
        //    This results in a uniform distribution of values in the unit circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the matrix.
        //
        //    Input/output, int &SEED, the "seed" value, which should NOT be 0.
        //    On output, SEED has been updated.
        //
        //    Output, complex <double> C8MAT_UNIFORM_01_NEW[M*N], the pseudorandom 
        //    complex matrix.
        //
    {
        int j;

        Complex[] c = new Complex[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                int k = seed / 127773;

                seed = 16807 * (seed - k * 127773) - k * 2836;

                switch (seed)
                {
                    case < 0:
                        seed += 2147483647;
                        break;
                }

                double r = Math.Sqrt(seed * 4.656612875E-10);

                k = seed / 127773;

                seed = 16807 * (seed - k * 127773) - k * 2836;

                switch (seed)
                {
                    case < 0:
                        seed += 2147483647;
                        break;
                }

                double theta = 2.0 * Math.PI * (seed * 4.656612875E-10);

                c[i + j * m] = r * new Complex(Math.Cos(theta), Math.Sin(theta));
            }
        }

        return c;
    }
}