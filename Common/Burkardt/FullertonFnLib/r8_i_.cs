using System;

namespace Burkardt.FullertonFnLib
{
    public static partial class FullertonLib
    {
        public static int r8_inits(double[] dos, int nos, double eta )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_INITS initializes a Chebyshev series.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 April 2016
        //
        //  Author:
        //
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Roger Broucke,
        //    Algorithm 446:
        //    Ten Subroutines for the Manipulation of Chebyshev Series,
        //    Communications of the ACM,
        //    Volume 16, Number 4, April 1973, pages 254-256.
        //
        //  Parameters:
        //
        //    Input, double DOS[NOS], the Chebyshev coefficients.
        //
        //    Input, int NOS, the number of coefficients.
        //
        //    Input, double ETA, the desired accuracy.
        //
        //    Output, int R8_INITS, the number of terms of the series needed
        //    to ensure the requested accuracy.
        //
        {
            double err;
            int i;
            int value = 0;

            if (nos < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("R8_INITS - Fatal error!");
                Console.WriteLine("  Number of coefficients < 1.");
                return (1);
            }

            if (eta < dos[nos - 1])
            {
                Console.WriteLine("");
                Console.WriteLine("R8_INITS - Warning!");
                Console.WriteLine("  ETA may be too small.");
                Console.WriteLine("  The requested accuracy cannot be guaranteed");
                Console.WriteLine("  even if all available coefficients are used.");
                value = nos;
            }
            else
            {
                err = 0.0;

                for (i = nos - 1; 0 <= i; i--)
                {
                    value = i + 1;
                    err = err + Math.Abs(dos[i]);
                    if (eta < err)
                    {
                        break;
                    }
                }
            }

            return value;
        }

        public static double r8_int(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_INT returns the integer part of an R8 argument.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 September 2011
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Wayne Fullerton.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Wayne Fullerton,
            //    Portable Special Function Routines,
            //    in Portability of Numerical Software,
            //    edited by Wayne Cowell,
            //    Lecture Notes in Computer Science, Volume 57,
            //    Springer 1977,
            //    ISBN: 978-3-540-08446-4,
            //    LC: QA297.W65.
            //
            //  Parameters:
            //
            //    Input, double X, the argument.
            //
            //    Output, double R8_INT, the integer part of X.
            //
        {
            int i;
            int ibase;
            int ipart;
            int npart = 0;
            double part;
            double scale = 0.0;
            double value;
            double xbig = 0.0;
            double xmax = 0.0;
            double xscl;

            if (npart == 0)
            {
                ibase = i4_mach(10);
                xmax = 1.0 / r8_mach(4);
                xbig = r8_min((double) (i4_mach(9)), 1.0 / r8_mach(4));
                scale = (double) i4_pow(ibase,
                    (int) (Math.Log(xbig) / Math.Log((double) (ibase)) - 0.5));
                npart = (int)(Math.Log(xmax) / Math.Log(scale) + 1.0);
            }

            //
            //  X may be too small.
            //
            if (x < -xmax)
            {
                value = x;
            }
            else if (x < -xbig)
            {
                xscl = -x;

                for (i = 1; i <= npart; i++)
                {
                    xscl = xscl / scale;
                }

                value = 0.0;
                for (i = 1; i <= npart; i++)
                {
                    xscl = xscl * scale;
                    ipart = (int) (xscl);
                    part = (double) (ipart);
                    xscl = xscl - part;
                    value = value * scale + part;
                }

                value = -value;
            }
            else if (x <= xbig)
            {
                value = (int) (x);
            }
            else if (x <= xmax)
            {
                xscl = x;

                for (i = 1; i <= npart; i++)
                {
                    xscl = xscl / scale;
                }

                value = 0.0;
                for (i = 1; i <= npart; i++)
                {
                    xscl = xscl * scale;
                    ipart = (int) (xscl);
                    part = (double) (ipart);
                    xscl = xscl - part;
                    value = value * scale + part;
                }
            }
            //
            //  X may be too large.
            //
            else
            {
                value = x;
            }

            return value;
        }
    }
}