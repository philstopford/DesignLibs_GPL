using System;
using geoCoreLib;

class Program 
{
    static void Main(string[] args)
    {
        try
        {
            string testFile = "/home/runner/work/DesignLibs_GPL/DesignLibs_GPL/geocore_test/cellref_array_irregular.oas";
            Console.WriteLine($"Testing OASIS file: {testFile}");
            
            GeoCoreHandler gH = new GeoCoreHandler();
            gH.updateGeoCoreHandler(testFile, GeoCore.fileType.oasis);
            GeoCore gc = gH.getGeo();
            
            if (gc.isValid())
            {
                Console.WriteLine("✅ OASIS file loaded successfully!");
                GCDrawingfield drawing = gc.getDrawing();
                Console.WriteLine($"Drawing has {drawing.cellList.Count} cells");
                foreach (var cell in drawing.cellList)
                {
                    Console.WriteLine($"  Cell: {cell.cellName} with {cell.elementList.Count} elements");
                }
            }
            else
            {
                Console.WriteLine("❌ OASIS file failed to load");
                foreach (string error in gc.error_msgs)
                {
                    Console.WriteLine($"Error: {error}");
                }
            }
        }
        catch (Exception ex)
        {
            Console.WriteLine($"❌ Exception: {ex.Message}");
        }
    }
}