using System;
using System.IO;
using System.Linq;
using System.Reflection;

namespace TestHelpers
{
    public static class TestOutput
    {
        private const string EnvVarDefault = "DECOMPOSITION_OUTPUT";
        private static readonly string[] RepoMarkers = new[] { ".git", "*.sln", "UnitTests.csproj" };

        // Convenience default: decomposition_out
        public static string GetPath()
        {
            return GetPath("decomposition_out");
        }

        // Get a writable path for a named test output directory.
        // Resolution order:
        // 1) Environment variables (NAME, NAME_OUTPUT, NAME_OUT, DECOMPOSITION_OUTPUT)
        // 2) GITHUB_WORKSPACE/artifacts/<name>
        // 3) Walk up from test assembly directory to repo root and use <repo>/artifacts/<name>
        // 4) Fallback to system temp
        public static string GetPath(string name)
        {
            if (string.IsNullOrWhiteSpace(name))
                name = "decomposition_out";

            string normalized = name.Trim();
            string up = normalized.ToUpperInvariant();
            string[] candidates = new[] { up, up + "_OUTPUT", up + "_OUT", EnvVarDefault };

            foreach (var envName in candidates)
            {
                try
                {
                    var env = Environment.GetEnvironmentVariable(envName);
                    if (!string.IsNullOrWhiteSpace(env))
                    {
                        var p = Path.GetFullPath(env);
                        Directory.CreateDirectory(p);
                        return p;
                    }
                }
                catch { }
            }

            var gh = Environment.GetEnvironmentVariable("GITHUB_WORKSPACE");
            if (!string.IsNullOrWhiteSpace(gh))
            {
                var p = Path.Combine(gh, "artifacts", normalized);
                Directory.CreateDirectory(p);
                return Path.GetFullPath(p);
            }

            var start = GetTestAssemblyDirectory();
            if (!string.IsNullOrEmpty(start))
            {
                var dir = new DirectoryInfo(start);
                while (dir != null)
                {
                    try
                    {
                        if (dir.EnumerateFileSystemInfos().Any(fi => RepoMarkers.Any(marker =>
                            marker.StartsWith("*.") ? fi.Extension.Equals(marker.Substring(1), StringComparison.OrdinalIgnoreCase) : fi.Name.Equals(marker, StringComparison.OrdinalIgnoreCase))))
                        {
                            var p = Path.Combine(dir.FullName, "artifacts", normalized);
                            Directory.CreateDirectory(p);
                            return Path.GetFullPath(p);
                        }

                        if (Directory.Exists(Path.Combine(dir.FullName, ".git")))
                        {
                            var p = Path.Combine(dir.FullName, "artifacts", normalized);
                            Directory.CreateDirectory(p);
                            return Path.GetFullPath(p);
                        }
                    }
                    catch { }

                    dir = dir.Parent;
                }
            }

            var tmp = Path.Combine(Path.GetTempPath(), normalized);
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