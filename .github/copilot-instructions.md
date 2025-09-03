# DesignLibs_GPL

DesignLibs_GPL is a .NET 8.0 C# library collection providing geometric, mathematical, and graphics libraries primarily built for use by the Variance and Quilt tools. The solution contains 30+ projects including core libraries, Eto.Forms graphics components, unit tests, and prototyping applications.

Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.

## Working Effectively

- Bootstrap, build, and test the repository:
  - `dotnet restore` -- takes ~18 seconds. Expect warnings about WPF projects on Linux and Eto.Forms CI package versions. NEVER CANCEL.
  - `dotnet build UnitTests/UnitTests.csproj` -- takes ~3 seconds for core libraries build. NEVER CANCEL.
  - `dotnet test UnitTests/UnitTests.csproj` -- takes ~2 minutes, runs 430 tests. NEVER CANCEL. Set timeout to 5+ minutes.

- Platform-specific limitations:
  - **Windows-only projects**: `Eto.Veldrid.Wpf` and `TestEtoVeldrid.Wpf` cannot build on Linux
  - **Graphics dependencies**: Veldrid SPIRV native libraries may be missing on some platforms
  - **Build workaround**: Use `dotnet build UnitTests/UnitTests.csproj` instead of full solution build to avoid platform issues

- Build individual library projects:
  - `dotnet build Common/geoCore/geoCore.csproj` -- builds core geometry library (~2 seconds)
  - `dotnet build Common/geoLib/geoLib.csproj` -- builds basic geometric primitives
  - `dotnet build Common/clipper/clipper.csproj` -- builds polygon clipping library
  - All individual library builds complete in 1-3 seconds

## Validation

- Always test the core functionality after making changes:
  - `dotnet run --project HistoTest/HistoTest.csproj` -- runs successfully, demonstrates statistical computations
  - Tests should show histogram output of random number distribution
  - If this fails, core mathematical libraries have issues

- Always run the unit test suite to validate changes:
  - `dotnet test UnitTests/UnitTests.csproj` -- NEVER CANCEL, takes ~2 minutes
  - Expect 327 passing tests, 103 failing tests (known failures in clipper and shape engine tests)
  - New failures beyond the known 103 indicate regressions

- Always check code formatting before committing:
  - `dotnet format --verify-no-changes UnitTests/UnitTests.csproj` -- takes ~30 seconds
  - Expect many whitespace formatting violations in existing code
  - Run `dotnet format UnitTests/UnitTests.csproj` to fix formatting issues

## Common Tasks

The following are validated commands and their expected behavior:

### Repository Structure
```
DesignLibs_GPL/
├── Common/                          # Core library projects
│   ├── clipper/                     # Polygon clipping (Clipper2 C# port)
│   ├── geoCore/                     # GDSII/Oasis file parsing
│   ├── geoLib/                      # Basic geometric primitives
│   ├── geoWrangler/                 # Geometry transformations and operations
│   ├── Color/                       # Color representations for viewports
│   ├── LibTessDotNet/              # Tessellation library
│   ├── MersenneTwister/            # Random number generation
│   ├── Noise/                      # Perlin/Simplex noise functions
│   └── [other libraries]/
├── Eto/                            # Eto.Forms UI components (some Windows-only)
│   ├── errorReporter/              # Cross-platform error handling
│   ├── Eto.VeldridSurface/        # 2D graphics viewport
│   └── [other UI components]/
├── UnitTests/                      # NUnit test project
├── HistoTest/                      # Working example application
└── Prototyping/                    # Various test applications
```

### Build Commands
```bash
# Full restore (expect warnings)
dotnet restore                      # ~18 seconds, NEVER CANCEL

# Build core libraries successfully  
dotnet build UnitTests/UnitTests.csproj    # ~3 seconds

# Build individual libraries
dotnet build Common/geoCore/geoCore.csproj # ~2 seconds
dotnet build Common/clipper/clipper.csproj # ~2 seconds

# AVOID: Full solution build fails due to platform-specific projects
# dotnet build  # FAILS on Linux due to WPF dependencies
```

### Test Commands
```bash
# Run full test suite
dotnet test UnitTests/UnitTests.csproj     # ~2 minutes, NEVER CANCEL, timeout 5+ minutes
# Expected: 327 passed, 103 failed (known failures)

# Run working example application
dotnet run --project HistoTest/HistoTest.csproj    # ~2 seconds
# Expected output: Histogram of random number distribution
```

### Validation Commands
```bash
# Check code formatting
dotnet format --verify-no-changes UnitTests/UnitTests.csproj   # ~30 seconds
# Expected: Many whitespace violations in existing code

# Fix formatting issues
dotnet format UnitTests/UnitTests.csproj   # ~30 seconds
```

### Known Issues and Workarounds
```bash
# Issue: Windows Desktop SDK missing on Linux
# Workaround: Build individual projects instead of full solution

# Issue: Veldrid SPIRV native dependencies missing
# Workaround: Use -p:NoPackageVeldrid=true for builds that don't need graphics

# Issue: Eto.Forms CI package version mismatches
# Status: Warnings only, builds still succeed

# Issue: 103 failing unit tests in clipper and shape engine
# Status: Known issues, not regressions unless count increases
```

## Important Library Information

- **geoCore**: GDSII and Oasis file format parsing and writing
- **geoLib**: Basic geometric primitives (points, rectangles, matrices) 
- **geoWrangler**: Geometry transformations, raycasting, fragmentation, decimation
- **clipper**: Polygon clipping operations using Clipper2 library
- **shapeEngine**: Shape processing and analysis
- **LibTessDotNet**: Triangle tessellation for complex polygons
- **Eto.VeldridSurface**: Cross-platform 2D graphics viewport for end-user applications

## Development Workflow

1. Always run `dotnet restore` first (NEVER CANCEL, ~18 seconds)
2. Build with `dotnet build UnitTests/UnitTests.csproj` (NEVER CANCEL, ~3 seconds)  
3. Test with `dotnet test UnitTests/UnitTests.csproj` (NEVER CANCEL, ~2 minutes)
4. Validate functionality with `dotnet run --project HistoTest/HistoTest.csproj`
5. Check formatting with `dotnet format --verify-no-changes UnitTests/UnitTests.csproj`
6. Fix formatting if needed with `dotnet format UnitTests/UnitTests.csproj`

CRITICAL: Always use timeouts of 5+ minutes for test commands and never cancel long-running operations. The test suite genuinely takes ~2 minutes to complete all 430 tests.