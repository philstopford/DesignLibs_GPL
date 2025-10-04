using System;
using System.IO;
using System.Linq;
using System.Reflection;

namespace TestHelpers
{
    /// <summary>
    /// Helper class to get test output paths for various test suites.
    /// Provides consistent path resolution across different environments (local dev, CI, etc.)
    /// </summary>
    public static class TestPath
    {
        private static readonly string[] RepoMarkers = new[] { ".git", "*.sln", "UnitTests.csproj" };

        /// <summary>
        /// Gets the output path for a specific test subdirectory (e.g., "decomposition_out", "geocore_out", etc.)
        /// </summary>
        /// <param name="subdirectory">The subdirectory name (e.g., "decomposition_out", "geocore_out")</param>
        /// <returns>Full path to the test output directory</returns>
        public static string Get(string subdirectory)
        {
            // Check for specific environment variable (e.g., DECOMPOSITION_OUTPUT, GEOCORE_OUTPUT)
            var envVarName = subdirectory.ToUpperInvariant().Replace("_OUT", "_OUTPUT");
            var env = Environment.GetEnvironmentVariable(envVarName);
            if (!string.IsNullOrWhiteSpace(env))
            {
                var p = Path.GetFullPath(env);
                Directory.CreateDirectory(p);
                return p;
            }

            // Check for generic TEST_OUTPUT environment variable
            var testOutput = Environment.GetEnvironmentVariable("TEST_OUTPUT");
            if (!string.IsNullOrWhiteSpace(testOutput))
            {
                var p = Path.Combine(Path.GetFullPath(testOutput), subdirectory);
                Directory.CreateDirectory(p);
                return p;
            }

            // Check for GitHub Actions workspace
            var gh = Environment.GetEnvironmentVariable("GITHUB_WORKSPACE");
            if (!string.IsNullOrWhiteSpace(gh))
            {
                var p = Path.Combine(gh, "artifacts", subdirectory);
                Directory.CreateDirectory(p);
                return Path.GetFullPath(p);
            }

            // Try to find repository root and use artifacts folder
            var start = GetTestAssemblyDirectory();
            if (!string.IsNullOrEmpty(start))
            {
                var dir = new DirectoryInfo(start);
                while (dir != null)
                {
                    if (dir.EnumerateFileSystemInfos().Any(fi => RepoMarkers.Any(marker =>
                        marker.StartsWith("*.") ? fi.Extension.Equals(marker.Substring(1), StringComparison.OrdinalIgnoreCase) : fi.Name.Equals(marker, StringComparison.OrdinalIgnoreCase))))
                    {
                        var p = Path.Combine(dir.FullName, "artifacts", subdirectory);
                        Directory.CreateDirectory(p);
                        return Path.GetFullPath(p);
                    }

                    if (Directory.Exists(Path.Combine(dir.FullName, ".git")))
                    {
                        var p = Path.Combine(dir.FullName, "artifacts", subdirectory);
                        Directory.CreateDirectory(p);
                        return Path.GetFullPath(p);
                    }

                    dir = dir.Parent;
                }
            }

            // Fallback to temp directory
            var tmp = Path.Combine(Path.GetTempPath(), subdirectory);
            Directory.CreateDirectory(tmp);
            return Path.GetFullPath(tmp);
        }

        /// <summary>
        /// Backward compatibility method for decomposition tests
        /// </summary>
        public static string GetDecompositionPath() => Get("decomposition_out");

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

    /// <summary>
    /// Backward compatibility - use TestPath instead
    /// </summary>
    [Obsolete("Use TestPath.GetDecompositionPath() or TestPath.Get(\"decomposition_out\") instead")]
    public static class DecompositionOutput
    {
        public static string GetPath() => TestPath.GetDecompositionPath();
    }
}
