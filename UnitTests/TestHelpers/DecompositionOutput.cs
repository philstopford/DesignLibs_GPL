using System;
using System.IO;
using System.Linq;
using System.Reflection;

namespace TestHelpers
{
    public static class DecompositionOutput
    {
        private const string EnvVar = "DECOMPOSITION_OUTPUT";
        private static readonly string[] RepoMarkers = new[] { ".git", "*.sln", "UnitTests.csproj" };

        public static string GetPath()
        {
            var env = Environment.GetEnvironmentVariable(EnvVar);
            if (!string.IsNullOrWhiteSpace(env))
            {
                var p = Path.GetFullPath(env);
                Directory.CreateDirectory(p);
                return p;
            }

            var gh = Environment.GetEnvironmentVariable("GITHUB_WORKSPACE");
            if (!string.IsNullOrWhiteSpace(gh))
            {
                var p = Path.Combine(gh, "artifacts", "decomposition_out");
                Directory.CreateDirectory(p);
                return Path.GetFullPath(p);
            }

            var start = GetTestAssemblyDirectory();
            if (!string.IsNullOrEmpty(start))
            {
                var dir = new DirectoryInfo(start);
                while (dir != null)
                {
                    if (dir.EnumerateFileSystemInfos().Any(fi => RepoMarkers.Any(marker =>
                        marker.StartsWith("*.") ? fi.Extension.Equals(marker.Substring(1), StringComparison.OrdinalIgnoreCase) : fi.Name.Equals(marker, StringComparison.OrdinalIgnoreCase))))
                    {
                        var p = Path.Combine(dir.FullName, "artifacts", "decomposition_out");
                        Directory.CreateDirectory(p);
                        return Path.GetFullPath(p);
                    }

                    if (Directory.Exists(Path.Combine(dir.FullName, ".git")))
                    {
                        var p = Path.Combine(dir.FullName, "artifacts", "decomposition_out");
                        Directory.CreateDirectory(p);
                        return Path.GetFullPath(p);
                    }

                    dir = dir.Parent;
                }
            }

            var tmp = Path.Combine(Path.GetTempPath(), "decomposition_out");
            Directory.CreateDirectory(tmp);
            return Path.GetFullPath(tmp);
        }

        private static string GetTestAssemblyDirectory()
        {
            try
            {
                var nunitType = Type.GetType("NUnit.Framework.TestContext, nunit.framework");
                if (nunitType != null)
                {
                    var currentContextProp = nunitType.GetProperty("CurrentContext", BindingFlags.Public | BindingFlags.Static);
                    var ctx = currentContextProp?.GetValue(null);
                    var testDirProp = ctx?.GetType().GetProperty("TestDirectory", BindingFlags.Public | BindingFlags.Instance);
                    var testDir = testDirProp?.GetValue(ctx) as string;
                    if (!string.IsNullOrEmpty(testDir))
                        return Path.GetFullPath(testDir);
                }
            }
            catch { }

            var baseDir = AppContext.BaseDirectory;
            if (!string.IsNullOrEmpty(baseDir))
                return Path.GetFullPath(baseDir);

            var asm = Assembly.GetExecutingAssembly();
            var asmLocation = asm?.Location;
            if (!string.IsNullOrEmpty(asmLocation))
                return Path.GetDirectoryName(Path.GetFullPath(asmLocation));

            return null;
        }
    }
}
