using System;

namespace Burkardt.Cube
{
    public class Integrals
    {
        public static double cube_monomial(double[] a, double[] b, int[] expon )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE_MONOMIAL integrates a monomial over a cube in 3D.
        //
        //  Discussion:
        //
        //    This routine integrates a monomial of the form
        //
        //      product ( 1 <= dim <= 3 ) x(dim)^expon(dim)
        //
        //    The combination 0^0 should be treated as 1.
        //
        //    The integration region is:
        //      A(1) <= X <= B(1)
        //      A(2) <= Y <= B(2)
        //      A(3) <= Z <= B(3)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[3], B[3], the lower and upper limits.
        //
        //    Input, int EXPON[3], the exponents.
        //
        //    Output, double CUBE_MONOMIAL, the integral of the monomial.
        //
        {
            int i;
            double value;

            for (i = 0; i < 3; i++)
            {
                if (expon[i] == -1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("CUBE_MONOMIAL - Fatal error!");
                    Console.WriteLine("  Exponent of -1 encountered.");
                    return 1;
                }
            }

            value = 1.0;

            for (i = 0; i < 3; i++)
            {
                value = value
                        * (Math.Pow(b[i], expon[i] + 1) - Math.Pow(a[i], expon[i] + 1))
                        / (double) (expon[i] + 1);
            }

            return value;
        }

        public static void cube_monomial_test(int degree_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CUBE_MONOMIAL_TEST tests CUBE_MONOMIAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 September 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DEGREE_MAX, the maximum total degree of the
            //    monomials to check.
            //
        {
            double[] a =  {
                -1.0, -1.0, -1.0
            }
            ;
            int alpha;
            double[] b =  {
                +1.0, +1.0, +1.0
            }
            ;
            int beta;
            int[] expon = new int[3];
            int gamma;
            double value;

            Console.WriteLine("");
            Console.WriteLine("CUBE_MONOMIAL_TEST");
            Console.WriteLine("  For the unit hexahedron,");
            Console.WriteLine("  CUBE_MONOMIAL returns the exact value of the");
            Console.WriteLine("  integral of X^ALPHA Y^BETA Z^GAMMA");
            Console.WriteLine("");
            Console.WriteLine("  Volume = " + cube_volume(a, b) + "");
            Console.WriteLine("");
            Console.WriteLine("     ALPHA      BETA     GAMMA      INTEGRAL");
            Console.WriteLine("");

            for (alpha = 0; alpha <= degree_max; alpha++)
            {
                expon[0] = alpha;
                for (beta = 0; beta <= degree_max - alpha; beta++)
                {
                    expon[1] = beta;
                    for (gamma = 0; gamma <= degree_max - alpha - beta; gamma++)
                    {
                        expon[2] = gamma;

                        value = cube_monomial(a, b, expon);

                        Console.WriteLine("  " + expon[0].ToString().PadLeft(8)
                            + "  " + expon[1].ToString().PadLeft(8)
                            + "  " + expon[2].ToString().PadLeft(8)
                            + "  " + value.ToString().PadLeft(14) + "");
                    }
                }
            }
        }
        
        public static double cube_volume ( double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE_VOLUME: volume of a cube in 3D.
        //
        //  Discussion:
        //
        //    The integration region is:
        //      A(1) <= X <= B(1)
        //      A(2) <= Y <= B(2)
        //      A(3) <= Z <= B(3)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[3], B[3], the lower and upper limits.
        //
        //    Output, double CUBE_VOLUME, the volume.
        //
        {
            int i;
            double value;

            value = 1.0;

            for ( i = 0; i < 3; i++ )
            {
                value = value * ( b[i] - a[i] );
            }

            return value;
        }
    }
}