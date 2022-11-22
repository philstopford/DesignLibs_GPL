using Clipper2Lib;
using geoWrangler;
using utility;

namespace shapeEngine;

public static class distortShape
{
    public static PathsD distortion(PathsD input, bool[] drawn, decimal lDC1, decimal lDC2, double resolution, int scaleFactorForOperation)
    {
        Fragmenter fragment = new(resolution, scaleFactorForOperation);
        PathsD ret = new ();
        for (int poly = 0; poly < input.Count; poly++)
        {
            ret.Add(new (input[poly]));
            int poly1 = poly;
            switch (drawn[poly])
            {
                // Now let's get some barrel distortion sorted out. Only for non-drawn polygons, and skip if both coefficients are zero to avoid overhead.
                case false when lDC1 != 0 || lDC2 != 0:
                {
                    int pCount = input[poly].Count;
#if !SHAPEENGINESINGLETHREADED
                    Parallel.For(0, pCount, point =>
#else
                    for (Int32 point = 0; point < pCount; point++)
#endif
                        {
                            double px = ret[poly1][point].x;
                            double py = ret[poly1][point].y;

                            // Need to calculate a new 'radius' from the origin for each point in the polygon, then scale the X/Y values accordingly in the polygon.
                            // Use scale factor to try and guarantee a -1 to +1 value range
                            px /= scaleFactorForOperation;
                            py /= scaleFactorForOperation;

                            double oRadius = Math.Sqrt(Utils.myPow(px, 2) + Utils.myPow(py, 2));
                            // Polynomial radial distortion.
                            // rd = r(1 + (k1 * r^2) + (k2 * r^4)) from Zhang, 1999 (https://www.microsoft.com/en-us/research/wp-content/uploads/2016/11/zhan99.pdf)
                            // we only want a scaling factor for our X, Y coordinates.
                            // '1 -' or '1 +' drive the pincushion/barrel tone. Coefficients being negative will have the same effect, so just pick a direction and stick with it.
                            int amplifier = 1000; // scales up end-user values to work within this approach.
                            double t1 = Convert.ToDouble(lDC1) * amplifier * Utils.myPow(Math.Abs(oRadius), 2);
                            double t2 = Convert.ToDouble(lDC1) * Utils.myPow(amplifier, 2) * Utils.myPow(Math.Abs(oRadius), 4);
                            double sFactor = 1 - (t1 + t2);

                            px *= sFactor * scaleFactorForOperation;
                            py *= sFactor * scaleFactorForOperation;

                            ret[poly1][point] = new (px, py);

                        }
#if !SHAPEENGINESINGLETHREADED
                    );
#endif

                    // Re-fragment
                    ret[poly] = fragment.fragmentPath(ret[poly]);
                    break;
                }
            }
        }

        return ret;
    }
    
}