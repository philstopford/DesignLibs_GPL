using System;
using System.Collections.Generic;
using System.Text;

namespace geoLib
{
    public class GeoLibArray
    {
        public GeoLibPoint point { get; set; }
        public GeoLibPoint pitch { get; set; }
        public GeoLibPoint count { get; set; }
    }

    public class GeoLibArrayF
    {
        public GeoLibPointF point { get; set; }
        public GeoLibPointF pitch { get; set; }
        public GeoLibPoint count { get; set; }
    }
}
