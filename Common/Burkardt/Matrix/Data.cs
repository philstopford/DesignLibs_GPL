using System;
using Burkardt.PolynomialNS;
using Burkardt.Types;

namespace Burkardt.MatrixNS;

public static class Data
{
    public static void data_to_dif(int ntab, double[] xtab, double[] ytab, ref double[] diftab )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DATA_TO_DIF sets up a divided difference table from raw data.
        //
        //  Discussion:
        //
        //    Space can be saved by using a single array for both the DIFTAB and
        //    YTAB dummy parameters.  In that case, the difference table will
        //    overwrite the Y data without interfering with the computation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Carl deBoor,
        //    A Practical Guide to Splines,
        //    Springer, 2001,
        //    ISBN: 0387953663,
        //    LC: QA1.A647.v27.
        //
        //  Parameters:
        //
        //    Input, int NTAB, the number of pairs of points
        //    (XTAB[I],YTAB[I]) which are to be used as data.
        //
        //    Input, double XTAB[NTAB], the X values at which data was taken.
        //    These values must be distinct.
        //
        //    Input, double YTAB[NTAB], the corresponding Y values.
        //
        //    Output, double DIFTAB[NTAB], the divided difference coefficients
        //    corresponding to the input (XTAB,YTAB).
        //
    {
        int i;
        int j;
        //
        //  Copy the data values into DIFTAB.
        //
        for (i = 0; i < ntab; i++)
        {
            diftab[i] = ytab[i];
        }

        //
        //  Make sure the abscissas are distinct.
        //
        for (i = 0; i < ntab; i++)
        {
            for (j = i + 1; j < ntab; j++)
            {
                switch (xtab[i] - xtab[j])
                {
                    case 0.0:
                        Console.WriteLine("");
                        Console.WriteLine("DATA_TO_DIF - Fatal error!");
                        Console.WriteLine("  Two entries of XTAB are equal!");
                        Console.WriteLine("  XTAB[%d] = " + xtab[i] + "");
                        Console.WriteLine("  XTAB[%d] = " + xtab[j] + "");
                        return;
                }
            }
        }

        //
        //  Compute the divided differences.
        //
        for (i = 1; i <= ntab - 1; i++)
        {
            for (j = ntab - 1; i <= j; j--)
            {
                diftab[j] = (diftab[j] - diftab[j - 1]) / (xtab[j] - xtab[j - i]);
            }
        }
    }

    public static double[] data_to_dif_new(int ntab, double[] xtab, double[] ytab )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DATA_TO_DIF_NEW sets up a divided difference table from raw data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Carl deBoor,
        //    A Practical Guide to Splines,
        //    Springer, 2001,
        //    ISBN: 0387953663,
        //    LC: QA1.A647.v27.
        //
        //  Parameters:
        //
        //    Input, int NTAB, the number of pairs of points
        //    (XTAB[I],YTAB[I]) which are to be used as data.
        //
        //    Input, double XTAB[NTAB], the X values at which data was taken.
        //    These values must be distinct.
        //
        //    Input, double YTAB[NTAB], the corresponding Y values.
        //
        //    Output, double DATA_TO_DIF_NEW[NTAB], the divided difference coefficients
        //    corresponding to the input.
        //
    {
        double[] diftab;
        int i;
        int j;
        //
        //  Make sure the abscissas are distinct.
        //
        for (i = 0; i < ntab; i++)
        {
            for (j = i + 1; j < ntab; j++)
            {
                switch (xtab[i] - xtab[j])
                {
                    case 0.0:
                        Console.WriteLine("");
                        Console.WriteLine("DATA_TO_DIF_NEW - Fatal error!");
                        Console.WriteLine("  Two entries of XTAB are equal!");
                        Console.WriteLine("  XTAB[%d] = " + xtab[i] + "");
                        Console.WriteLine("  XTAB[%d] = " + xtab[j] + "");
                        return null;
                }
            }
        }

        //
        //  Copy the Y data into DIFTAB.
        //
        diftab = new double[ntab];

        for (i = 0; i < ntab; i++)
        {
            diftab[i] = ytab[i];
        }

        //
        //  Compute the divided differences.
        //
        for (i = 1; i <= ntab - 1; i++)
        {
            for (j = ntab - 1; i <= j; j--)
            {
                diftab[j] = (diftab[j] - diftab[j - 1])
                            / (xtab[j] - xtab[j - i]);
            }
        }

        return diftab;
    }

    public static void data_to_dif_display(int ntab, double[] xtab, double[] ytab,
            ref double[] diftab )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DATA_TO_DIF_DISPLAY sets up a divided difference table and prints intermediate data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NTAB, the number of pairs of points
        //    (XTAB[I],YTAB[I]) which are to be used as data.
        //
        //    Input, double XTAB[NTAB], the X values at which data was taken.
        //    These values must be distinct.
        //
        //    Input, double YTAB[NTAB], the corresponding Y values.
        //
        //    Output, double DIFTAB[NTAB], the divided difference coefficients
        //    corresponding to the input (XTAB,YTAB).
        //
    {
        int i;
        int j;

        if (!typeMethods.r8vec_distinct(ntab, xtab))
        {
            Console.WriteLine("");
            Console.WriteLine("DATA_TO_DIF_DISPLAY - Fatal error!");
            Console.WriteLine("  Two entries of XTAB are equal!");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("  The divided difference table:");
        Console.WriteLine("");
        string cout = "        ";
        for (i = 0; i < ntab; i++)
        {
            cout += xtab[i].ToString().PadLeft(10) + "  ";
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        cout = 0.ToString().PadLeft(6) + "  ";
        for (i = 0; i < ntab; i++)
        {
            cout += ytab[i].ToString().PadLeft(10) + "  ";
        }

        Console.WriteLine(cout);
        //
        //  Copy the data values into DIFTAB.
        //
        for (i = 0; i < ntab; i++)
        {
            diftab[i] = ytab[i];
        }

        //
        //  Compute the divided differences.
        //
        for (i = 1; i <= ntab - 1; i++)
        {
            cout = i.ToString().PadLeft(6) + "  ";
            for (j = ntab - 1; i <= j; j--)
            {
                diftab[j] = (diftab[j] - diftab[j - 1]) / (xtab[j] - xtab[j - i]);
            }

            for (j = i; j < ntab; j++)
            {
                cout += diftab[j].ToString().PadLeft(10) + "  ";
            }

            Console.WriteLine(cout);
        }
    }

    public static void data_to_r8poly(int ntab, double[] xtab, double[] ytab, ref double[] c )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DATA_TO_R8POLY computes the coefficients of a polynomial interpolating data.
        //
        //  Discussion:
        //
        //    Space can be saved by using a single array for both the C and
        //    YTAB parameters.  In that case, the coefficients will
        //    overwrite the Y data without interfering with the computation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Carl deBoor,
        //    A Practical Guide to Splines,
        //    Springer, 2001,
        //    ISBN: 0387953663,
        //    LC: QA1.A647.v27.
        //
        //  Parameters:
        //
        //    Input, int NTAB, the number of data points.
        //
        //    Input, double XTAB[NTAB], YTAB[NTAB], the data values.
        //
        //    Output, double C[NTAB], the coefficients of the polynomial that passes
        //    through the data (XTAB,YTAB).  C(0) is the constant term.
        //
    {
        if (!typeMethods.r8vec_distinct(ntab, xtab))
        {
            Console.WriteLine("");
            Console.WriteLine("DATA_TO_R8POLY - Fatal error!");
            Console.WriteLine("  Two entries of XTAB are equal.");
            return;
        }

        data_to_dif(ntab, xtab, ytab, ref c);

        Dif.dif_to_r8poly(ntab, xtab, c, ref c);
    }
}