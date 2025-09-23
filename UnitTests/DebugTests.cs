using System;
using System.Collections.Generic;
using NUnit.Framework;
using geoCoreLib;
using utility;

namespace DebugTests
{
    [TestFixture]
    public class DebugConsistencyTests
    {
        [Test]
        public static void debug_gds_units()
        {
            string baseDir = "../../../../../geocore_test/";
            
            GeoCoreHandler gH_OAS = new();
            gH_OAS.updateGeoCoreHandler(baseDir + "/consistency/triage.oas", GeoCore.fileType.oasis);
            GeoCore gcOAS = gH_OAS.getGeo();
            Assert.That(gcOAS.isValid(), Is.True);

            Console.WriteLine($"Original OASIS drawing: databaseunits={gcOAS.getDrawing().databaseunits}, userunits={gcOAS.getDrawing().userunits}");
            Console.WriteLine($"Original drawing scale: {gcOAS.getDrawing().getDrawingScale()}");

            string outFile = "../../../../../geocore_out/debug_gds_units.gds";
            if (System.IO.File.Exists(outFile))
            {
                System.IO.File.Delete(outFile);
            }
            gds.gdsWriter gw = new(gcOAS, outFile);
            gw.save();

            GeoCoreHandler gH_GDS = new();
            gH_GDS.updateGeoCoreHandler(outFile, GeoCore.fileType.gds);
            GeoCore gcGDS = gH_GDS.getGeo();
            Assert.That(gcGDS.isValid(), Is.True);

            Console.WriteLine($"GDS drawing after load: databaseunits={gcGDS.getDrawing().databaseunits}, userunits={gcGDS.getDrawing().userunits}");
            Console.WriteLine($"GDS drawing scale after load: {gcGDS.getDrawing().getDrawingScale()}");
        }
    }
}