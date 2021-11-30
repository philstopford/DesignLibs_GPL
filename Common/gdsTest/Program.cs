using geoCoreLib;

namespace gdsTest;

internal static class Program
{
    private static void Main()
    {
        // string filename = "C:\\Users\\phil\\OneDrive\\Documents\\Visual Studio 2015\\Projects\\MonteCarlo_Csharp\\oasTest\\test.gds";
        const string filename = "C:\\Users\\phil\\OneDrive\\Documents\\Visual Studio 2015\\Projects\\MonteCarlo_Csharp\\gdsFiles\\clip1.gds";
        // string filename = "C:\\Users\\phil\\OneDrive\\Documents\\Visual Studio 2015\\Projects\\MonteCarlo_Csharp\\gdsFiles\\doetest.gds";

        GeoCore g = new();
        g.updateGeoCore(filename, GeoCore.fileType.gds);
    }
}