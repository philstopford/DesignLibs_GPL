using System;
using geoLib;
using geoWrangler;

namespace GeoWranglerTest
{
    class Program
    {
        static void Main(string[] args)
        {
            GeoLibPointF[] source = new GeoLibPointF[]
            {
                new GeoLibPointF(0.02985, 0.18999),
                new GeoLibPointF(0.00864, 0.21120),
                new GeoLibPointF(0.01217, 0.21474),
                new GeoLibPointF(0.01571, 0.21827),
                new GeoLibPointF(0.02278, 0.21120),
                new GeoLibPointF(0.02631, 0.21474),
                new GeoLibPointF(0.02985, 0.21827),
                new GeoLibPointF(0.00156, 0.24656),
                new GeoLibPointF(-0.00196, 0.24302),
                new GeoLibPointF(-0.00550, 0.23949),
                new GeoLibPointF(-0.01257, 0.24656),
                new GeoLibPointF(-0.04121, 0.21791),
                new GeoLibPointF(-0.04829, 0.21084),
                new GeoLibPointF(-0.05500, 0.20413),
                new GeoLibPointF(-0.06207, 0.21120),
                new GeoLibPointF(-0.05500, 0.21827),
                new GeoLibPointF(-0.06207, 0.22534),
                new GeoLibPointF(-0.06560, 0.22181),
                new GeoLibPointF(-0.06914, 0.21828),
                new GeoLibPointF(-0.07621, 0.22535),
                new GeoLibPointF(-0.07267, 0.22888),
                new GeoLibPointF(-0.06914, 0.23242),
                new GeoLibPointF(-0.06560, 0.23595),
                new GeoLibPointF(-0.06207, 0.23949),
                new GeoLibPointF(-0.05853, 0.24302),
                new GeoLibPointF(-0.05500, 0.24656),
                new GeoLibPointF(-0.04793, 0.23949),
                new GeoLibPointF(-0.05499, 0.23242),
                new GeoLibPointF(-0.05146, 0.22888),
                new GeoLibPointF(-0.04792, 0.22535),
                new GeoLibPointF(-0.01964, 0.25363),
                new GeoLibPointF(-0.02671, 0.26070),
                new GeoLibPointF(-0.02318, 0.26424),
                new GeoLibPointF(-0.01964, 0.26777),
                new GeoLibPointF(-0.06277, 0.31091),
                new GeoLibPointF(-0.05570, 0.31798),
                new GeoLibPointF(-0.04792, 0.31020),
                new GeoLibPointF(-0.04085, 0.31727),
                new GeoLibPointF(-0.04792, 0.32434),
                new GeoLibPointF(-0.04085, 0.33141),
                new GeoLibPointF(-0.01964, 0.31020),
                new GeoLibPointF(-0.02671, 0.30313),
                new GeoLibPointF(-0.03378, 0.31020),
                new GeoLibPointF(-0.04085, 0.30313),
                new GeoLibPointF(-0.01257, 0.27484),
                new GeoLibPointF(-0.00550, 0.28191),
                new GeoLibPointF(0.00156, 0.27484),
                new GeoLibPointF(0.04399, 0.31727),
                new GeoLibPointF(0.05106, 0.31020),
                new GeoLibPointF(0.04753, 0.30666),
                new GeoLibPointF(0.04399, 0.30313),
                new GeoLibPointF(0.05106, 0.29606),
                new GeoLibPointF(0.05813, 0.30313),
                new GeoLibPointF(0.06520, 0.29606),
                new GeoLibPointF(0.06167, 0.29252),
                new GeoLibPointF(0.05813, 0.28899),
                new GeoLibPointF(0.05460, 0.28545),
                new GeoLibPointF(0.04399, 0.27484),
                new GeoLibPointF(0.03692, 0.28191),
                new GeoLibPointF(0.04399, 0.28899),
                new GeoLibPointF(0.04045, 0.29252),
                new GeoLibPointF(0.03692, 0.29606),
                new GeoLibPointF(0.00864, 0.26777),
                new GeoLibPointF(0.01571, 0.26070),
                new GeoLibPointF(0.00864, 0.25363),
                new GeoLibPointF(0.05177, 0.21050),
                new GeoLibPointF(0.04470, 0.20343),
                new GeoLibPointF(0.03692, 0.21120),
                new GeoLibPointF(0.03338, 0.20767),
                new GeoLibPointF(0.02985, 0.20413),
                new GeoLibPointF(0.03692, 0.19706),
                new GeoLibPointF(0.03692, 0.19706)
            };

            GeoLibPointF[] cleaned = GeoWrangler.stripColinear(source, tolerance:0.2);

            int xx = 2;
        }
    }
}
