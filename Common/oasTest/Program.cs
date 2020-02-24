using geoCoreLib;

namespace oasisTest
{
    class Program
    {
        static void Main(string[] args)
        {
            // string filename = "C:\\Users\\phil\\OneDrive\\Documents\\Visual Studio 2015\\Projects\\MonteCarlo_Csharp\\oasTest\\test.oas";
            string filename = "C:\\Users\\phil\\OneDrive\\Documents\\Visual Studio 2015\\Projects\\MonteCarlo_Csharp\\geoCore\\oasTest\\clip2.oas";
            // string filename = "C:\\Users\\phil\\OneDrive\\Documents\\Visual Studio 2015\\Projects\\MonteCarlo_Csharp\\gdsFiles\\doetest.oas";

            GeoCore g = new GeoCore();
            g.updateGeoCore(filename, GeoCore.fileType.oasis);

            int HoldMe = 2;
        }
    }
}
