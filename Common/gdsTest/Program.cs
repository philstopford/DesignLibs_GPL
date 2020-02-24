using geoCoreLib;

namespace gdsTest
{
    class Program
    {
        static void Main(string[] args)
        {
            // string filename = "C:\\Users\\phil\\OneDrive\\Documents\\Visual Studio 2015\\Projects\\MonteCarlo_Csharp\\oasTest\\test.gds";
            string filename = "C:\\Users\\phil\\OneDrive\\Documents\\Visual Studio 2015\\Projects\\MonteCarlo_Csharp\\gdsFiles\\clip1.gds";
            // string filename = "C:\\Users\\phil\\OneDrive\\Documents\\Visual Studio 2015\\Projects\\MonteCarlo_Csharp\\gdsFiles\\doetest.gds";

            GeoCore g = new GeoCore();
            g.updateGeoCore(filename, GeoCore.fileType.gds);

            int HoldMe = 2;
        }
    }
}
