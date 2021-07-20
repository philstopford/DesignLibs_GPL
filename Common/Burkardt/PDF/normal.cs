using System;

namespace Burkardt.PDFLib
{
    public static partial class PDF
    {
        public static bool normal_check ( double mu, double sigma )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_CHECK checks the parameters of the Normal PDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double MU, SIGMA, the parameters of the PDF.
            //    0.0 < SIGMA.
            //
            //    Output, bool NORMAL_CHECK, is true if the parameters are legal.
            //
        {
            if ( sigma <= 0.0 )
            {
                Console.WriteLine("");
                Console.WriteLine("NORMAL_CHECK - Fatal error!");
                Console.WriteLine("  SIGMA <= 0.");
                return false;
            }

            return true;
        }
    }
}