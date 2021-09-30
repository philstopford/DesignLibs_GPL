using System;
using Burkardt;
using Burkardt.Quadrature;

namespace FastGaussLegendreTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for FASTGL_TEST.
            //
            //  Licensing:
            //
            //    This code is distributed under the BSD license. 
            //
            //  Modified:
            //
            //    22 December 2015
            //
            //  Author:
            //
            //    Ignace Bogaert
            //
            //  Reference:
            //
            //    Ignace Bogaert,
            //    Iteration-free computation of Gauss-Legendre quadrature nodes and weights,
            //    SIAM Journal on Scientific Computing,
            //    Volume 36, Number 3, 2014, pages A1008-1026.
            //
        {
            Console.WriteLine("");
            Console.WriteLine("FASTGL_TEST");
            Console.WriteLine("  Test the FASTGL library.");
            //
            // Some information on the origin of this code.
            //
            Console.WriteLine("");
            Console.WriteLine("This program shows usage examples for the Gauss-Legendre quadrature rules, computed with fastgl::GLPair(l, k)");
            Console.WriteLine("\t--> l is the number of nodes in the rule, k is the index of the node that will be computed.");
            Console.WriteLine();
            Console.WriteLine("The computation of the nodes and weights is based on the following paper:");
            Console.WriteLine("\tIgnace Bogaert, 'Iteration-Free Computation of Gauss-Legendre Quadrature Nodes and Weights',");
            Console.WriteLine("\tto appear in the SIAM Journal of Scientific Computing.");
            Console.WriteLine();
            Console.WriteLine("The main features of this software are:");
            Console.WriteLine("\t- Speed: due to the simple formulas and the O(1) complexity computation of individual");
            Console.WriteLine("\t  Gauss-Legendre quadrature nodes and weights. This also makes it perfectly compatible");
            Console.WriteLine("\t  with parallel computing paradigms such as multithreading and MPI.");
            Console.WriteLine("\t- Accuracy: the error on the nodes and weights is within a few ulps (see the paper for details).");
            //
            //  Test the numerical integration of Math.Exp(x) over the range [-1,1]
            //  for varying number of Gauss-Legendre quadrature nodes l.
            //
            Console.WriteLine("");
            Console.WriteLine("First test-case: int(Math.Exp(x), x = -1..1):");
            Console.WriteLine("");

            int l;
            
            for (l = 5; l <= 9; ++l)
            {
                double Res = 0.0;
                for (int k = 1; k <= l; ++k)
                {
                    FastGaussLegendre.QuadPair p = FastGaussLegendre.GLPair(l, k);
                    Res += p.weight * Math.Exp(p.x());
                }

                Console.WriteLine("Gauss-Legendre " + l + "-node result = " + Res);
            }

            Console.WriteLine("Exact result                 = " + (Math.Exp(1.0) - Math.Exp(-1.0)).ToString("0.################"));
            //	
            //  Test the numerical integration of cos(1000 x) over the range [-1,1]
            //  for varying number of Gauss-Legendre quadrature nodes l.
            //  The fact that only twelve digits of accuracy are obtained is due to the 
            //  condition number of the summation.
            //
            Console.WriteLine("");
            Console.WriteLine("Second test-case: int(cos(1000x), x = -1..1):");
            Console.WriteLine("");
            for (l = 500; l <= 600; l += 20)
            {
                double Res = 0.0;
                for (int k = 1; k <= l; ++k)
                {
                    FastGaussLegendre.QuadPair p = FastGaussLegendre.GLPair(l, k);
                    Res += p.weight * Math.Cos(1000.0 * p.x());
                }

                Console.WriteLine("Gauss-Legendre " + l + "-node result = " + Res);
            }

            Console.WriteLine("Exact result                   = " + 0.002 * Math.Sin(1000.0));
            //	
            //  Test the numerical integration of ln(x) over the range [0,1]
            //  Normally, one would not use Gauss-Legendre quadrature for this,
            //  but for the sake of having an example with l > 100, this is included.
            //
            Console.WriteLine("");
            Console.WriteLine("Third test-case: int(ln(x), x = 0..1):");
            Console.WriteLine("");
            l = 1;
            for (int p = 0; p <= 6; ++p)
            {
                double Res = 0.0;
                for (int k = 1; k <= l; ++k)
                {
                    FastGaussLegendre.QuadPair p_ = FastGaussLegendre.GLPair(l, k);
                    Res += 0.5 * p_.weight * Math.Log(0.5 * (p_.x() + 1.0));
                }

                Console.WriteLine("Gauss-Legendre " + l + "-node result = " + Res);
                l *= 10;
            }

            Console.WriteLine("Exact result                       = " + -1.0);
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("FASTGL_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}