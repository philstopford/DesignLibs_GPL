/////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2003 CenterSpace Software, LLC                            //
//                                                                         //
// This code is free software under the Artistic license.                  //
//                                                                         //
// CenterSpace Software                                                    //
// 301 SW 4th Street - Suite #240                                          //
// Corvallis, Oregon, 97333                                                //
// USA                                                                     //
// http://www.centerspace.net                                              //
/////////////////////////////////////////////////////////////////////////////
using System;
using utility;

namespace CenterSpace.Free;

/// <summary>
/// Summary description for Class1.
/// </summary>
internal class HistogramExample
{
    /// <summary>
    /// The main entry point for the application.
    /// </summary>
    [STAThread]
    private static void Main(string[] args)
    {
        int length = 500;
        double[] data = new double[length];
        Random randGen = new();
        for (int i = 0; i < length; ++i)
        {
            data[i] = randGen.NextDouble();
        }

        Histo h = new(20, 0, 1);
        h.AddData(data);

        double numLessThanOneHalf = 0, numGreaterThanOneHalf = 0;
        for (int i = 0; i < h.NumBins; ++i)
        {
            switch (h.BinBoundaries[i + 1])
            {
                case <= 0.5:
                    numLessThanOneHalf += h.Counts[i];
                    break;
                default:
                    numGreaterThanOneHalf += h.Counts[i];
                    break;
            }
        }

        Histo.FloatFormat = "F4";
        Console.WriteLine("{0} uniform random numbers in [0,1]", length);
        Console.WriteLine("  {0} <= 1/2, and {1} > 1/2", numLessThanOneHalf, numGreaterThanOneHalf);
        Console.WriteLine(h.StemLeaf());
    }
}