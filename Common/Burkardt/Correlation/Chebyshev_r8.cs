using System;

namespace Burkardt.CorrelationNS;

public static partial class Correlation
{
    public static int r8_inits ( double[] dos, int nos, double eta )

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
        //    15 September 2011
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
        int value;

        switch (nos)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("R8_INITS - Fatal error!");
                Console.WriteLine("  Number of coefficients < 1.");
                return 1;
        }

        err = 0.0;

        for ( i = nos - 1; 0 <= i; i-- )
        {
            err += Math.Abs ( dos[i] );
            if ( eta < err )
            {
                value = i + 1;
                return value;
            }
        }

        value = i;
        Console.WriteLine("");
        Console.WriteLine("R8_INITS - Warning!");
        Console.WriteLine("  ETA may be too small.");

        return value;
    }
        
    public static double r8_csevl ( double x, double[] a, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CSEVL evaluates a Chebyshev series.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 September 2011
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
        //    Input, double X, the evaluation point.
        //
        //    Input, double CS[N], the Chebyshev coefficients.
        //
        //    Input, int N, the number of Chebyshev coefficients.
        //
        //    Output, double R8_CSEVL, the Chebyshev series evaluated at X.
        //
    {
        double b0;
        double b1;
        double b2 = 0;
        int i;
        double twox;
        double value = 0;

        switch (n)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("R8_CSEVL - Fatal error!");
                Console.WriteLine("  Number of terms <= 0.");
                return 1;
            case > 1000:
                Console.WriteLine("");
                Console.WriteLine("R8_CSEVL - Fatal error!");
                Console.WriteLine("  Number of terms greater than 1000.");
                return 1;
        }

        switch (x)
        {
            case < -1.1:
            case > 1.1:
                Console.WriteLine("");
                Console.WriteLine("R8_CSEVL - Fatal error!");
                Console.WriteLine("  X outside (-1,+1).");
                return 1;
        }

        twox = 2.0 * x;
        b1 = 0.0;
        b0 = 0.0;

        for ( i = n - 1; 0 <= i; i-- )
        {
            b2 = b1;
            b1 = b0;
            b0 = twox * b1 - b2 + a[i];
        }

        value = 0.5 * ( b0 - b2 );

        return value;
    }
}