using System;
using Burkardt.Uniform;
using entropyRNG;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static int get_seed()
        {
            return RNG.nextint(1, Int32.MaxValue);
        }

        public static int[] perm_uniform ( int n, int base_, ref int seed )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM_UNIFORM selects a random permutation of N objects.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 October 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms,
            //    Academic Press, 1978, second edition,
            //    ISBN 0-12-519260-6.
            //
            //  Parameters:
            //
            //    Input, int N, the number of objects to be permuted.
            //
            //    Input, int BASE, is 0 for a 0-based permutation and 1 for 
            //    a 1-based permutation.
            //
            //    Input/output, int *SEED, a seed for the random number generator.
            //
            //    Output, int PERM_UNIFORM[N], a permutation of (BASE, BASE+1, ..., BASE+N-1).
            //
        {
            int[] p = new int[n];
 
            for (int i = 0; i < n; i++ )
            {
                p[i] = i + base_;
            }

            for (int i = 0; i < n; i++ )
            {
                int j = UniformRNG.i4_uniform( i, n - 1, ref seed );
                int k    = p[i];
                p[i] = p[j];
                p[j] = k;
            }
 
            return p;
        }
        
        public static int[] perm_uniform_new ( int n, ref int seed )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_UNIFORM_NEW selects a random permutation of N objects.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Albert Nijenhuis, Herbert Wilf,
        //    Combinatorial Algorithms for Computers and Calculators,
        //    Second Edition,
        //    Academic Press, 1978,
        //    ISBN: 0-12-519260-6,
        //    LC: QA164.N54.
        //
        //  Parameters:
        //
        //    Input, int N, the number of objects to be permuted.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int PERM_UNIFORM_NEW[N], a permutation of (1, 2, ..., N).
        //
        {
            int[] p = new int[n];

            for (int i = 0; i < n; i++ )
            {
                p[i] = i + 1;
            }

            for (int i = 0; i < n; i++ )
            {
                int j = UniformRNG.i4_uniform_ab ( i, n - 1, ref seed );
                int k = p[i];
                p[i] = p[j];
                p[j] = k;
            }
 
            return p;
        }
        
        public static void perm_check ( int n, int[] p )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_CHECK checks that a vector represents a permutation.
        //
        //  Discussion:
        //
        //    The routine verifies that each of the integers from 1
        //    to N occurs among the N entries of the permutation.
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
        //    Input, int N, the number of entries.
        //
        //    Input, int P[N], the permutation, in standard index form.
        //
        {
            bool error;
            int ifind;
            int iseek;

            for ( iseek = 1; iseek <= n; iseek++ )
            {
                error = true;

                for ( ifind = 1; ifind <= n; ifind++ )
                {
                    if ( p[ifind-1] == iseek )
                    {
                        error = false;
                        break;
                    }
                }

                if ( error )
                {
                    Console.WriteLine();
                    Console.WriteLine("PERM_CHECK - Fatal error!");
                    Console.WriteLine("  The input permutation is not legal.");
                }
            }
        }
        
        public static void perm_print ( int n, int[] p, string title )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_PRINT prints a permutation.
        //
        //  Discussion:
        //
        //    The permutation is assumed to be zero-based.
        //
        //  Example:
        //
        //    Input:
        //
        //      P = 6 1 2 0 4 2 5
        //
        //    Printed output:
        //
        //      "This is the permutation:"
        //
        //      0 1 2 3 4 5 6
        //      6 1 2 0 4 2 5
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 May 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of objects permuted.
        //
        //    Input, int P[N], the permutation, in standard index form.
        //
        //    Input, string TITLE, a title.
        //    If no title is supplied, then only the permutation is printed.
        //
        {
            int inc = 20;

            if ( s_len_trim ( title ) != 0 )
            {
                Console.WriteLine();
                Console.WriteLine(title);

                for (int ilo = 0; ilo < n; ilo = ilo + inc )
                {
                    int ihi = ilo + inc;
                    if ( n < ihi ) 
                    {
                        ihi = n;
                    }

                    Console.WriteLine();
                    string cout = "  ";
                    for (int i = ilo; i < ihi; i++ )
                    {
                        cout += i.ToString().PadLeft(4);
                    }
                    Console.WriteLine(cout);

                    Console.WriteLine();
                    cout = "  ";
                    for (int i = ilo; i < ihi; i++ )
                    {
                        cout += p[i].ToString().PadLeft(4);
                    }
                    Console.WriteLine(cout);
                    Console.WriteLine();
                }
            }
            else
            {
                for (int ilo = 0; ilo < n; ilo = ilo + inc )
                {
                    int ihi = ilo + inc;
                    if ( n < ihi ) 
                    {
                        ihi = n;
                    }
                    string cout = "  ";
                    for (int i = ilo; i < ihi; i++ )
                    {
                        cout += p[i].ToString().PadLeft(4);
                    }
                    Console.WriteLine(cout);
                    Console.WriteLine();
                }
            }
        }
        
        
        public static double[] monomial_value ( int m, int n, int[] e, double[] x )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONOMIAL_VALUE evaluates a monomial.
        //
        //  Discussion:
        //
        //    This routine evaluates a monomial of the form
        //
        //      product ( 1 <= i <= m ) x(i)^e(i)
        //
        //    where the exponents are nonnegative integers.  Note that
        //    if the combination 0^0 is encountered, it should be treated
        //    as 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N, the number of points at which the
        //    monomial is to be evaluated.
        //
        //    Input, int E[M], the exponents.
        //
        //    Input, double X[M*N], the point coordinates.
        //
        //    Output, double MONOMIAL_VALUE[N], the value of the monomial.
        //
        {
            double[] v = new double[n];

            for (int j = 0; j < n; j++ )
            {
                v[j] = 1.0;
            }

            for (int i = 0; i < m; i++ )
            {
                if ( 0 != e[i] )
                {
                    for (int j = 0; j < n; j++ )
                    {
                        v[j] = v[j] * Math.Pow ( x[i+j*m], e[i] );
                    }
                }
            }

            return v;
        }

        public static int inits ( double[] dos, int nos, double eta )
//****************************************************************************80
//
//  Purpose:
//
//    INITS initializes a Chebyshev series.
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
//    Output, int INITS, the number of terms of the series needed
//    to ensure the requested accuracy.
//
        {
            double err;
            int i;
            int value;

            if ( nos < 1 )
            {
                Console.WriteLine("");
                Console.WriteLine("INITS - Fatal error!");
                Console.WriteLine("  Number of coefficients < 1.");
                return ( 1 );
            }

            err = 0.0;

            for ( i = nos - 1; 0 <= i; i-- )
            {
                err = err + Math.Abs ( dos[i] );
                if ( eta < err )
                {
                    value = i + 1;
                    return value;
                }
            }

            value = i;
            Console.WriteLine("");
            Console.WriteLine("INITS - Warning!");
            Console.WriteLine("  ETA may be too small.");

            return value;
        }
        
    public static double csevl ( double x, double[] a, int n )
//****************************************************************************80
//
//  Purpose:
//
//    CSEVL evaluates a Chebyshev series.
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
//    Output, double CSEVL, the Chebyshev series evaluated at X.
//
        {
            double b0;
            double b1;
            double b2 = 0;
            int i;
            double twox;
            double value;

            if ( n < 1 )
            {
                Console.WriteLine("");
                Console.WriteLine("CSEVL - Fatal error!");
                Console.WriteLine("  Number of terms <= 0.");
                return ( 1 );
            }

            if ( 1000 < n )
            {
                Console.WriteLine("");
                Console.WriteLine("CSEVL - Fatal error!");
                Console.WriteLine("  Number of terms greater than 1000.");
                return ( 1 );
            }

            if ( x < -1.1 || 1.1 < x )
            {
                Console.WriteLine("");
                Console.WriteLine("CSEVL - Fatal error!");
                Console.WriteLine("  X outside (-1,+1).");
                return ( 1 );
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
}