# DesignLibs_GPL
GPLv3 release of DesignLibs

DesignLibs comprises a multitude of libraries. These are primarily built for use by the Variance (https://github.com/philstopford/Variance_GPL) and Quilt (https://github.com/philstopford/Quilt_GPL) tools. Gradually, documentation will be added to the wiki.

## Common

### Burkardt:

This is a WIP port of the huge resource available at https://people.sc.fsu.edu/~jburkardt/cpp_src/cpp_src.html. Validation is being tracked at https://github.com/philstopford/DesignLibs_GPL/issues/2


### Color:

This provides color representations and common variables that are used by the Eto.Forms-based viewports, and end-user tools.

### KDTree:

This is a copy of the KDTree formerly hosted at https://code.google.com/p/kd-sharp/

### LibTessDotNet:

This is a copy of https://github.com/speps/LibTessDotNet

### MersenneTwister:

This is a copy of http://www.centerspace.net/resources/free-stuff/mersenne-twister

### MiscUtil:

This is a modified copy of http://yoda.arachsys.com/csharp/miscutil/. The modification was the removal of some code that could not be supported under .NET Core 3.1, but which was not needed by the projects DesignLibs serves.

### Noise:

This provides Perlin, Simplex and OpenSimplex noise functions.

### Clipper:

This is a copy of the C# version of https://sourceforge.net/projects/polyclipping/. This is heavily used for polygon clipping, including raycasting.

### info.lundin.math

This is a copy of http://lundin.info/mathparser

### EntropyRNG:

This provides a threadsafe RNG system (interacting also optionally with the MersenneTwister).

### GeoLib:

This provides the basic geometric primitives (integer and floating point both). Points, arrays, matrices, and rectangles are available. These were written due to the multi-platform operation, where System.Drawing was not necessarily available, and relying on Eto.Forms introduced awkward dependencies on headless platforms.

### GeoWrangler:

This is a general purpose library to wrangle geometry of different types. Integer and floating points from GeoLib, as well as the Path(s) constructs from ClipperLib are all supported. This library allows for various transformations and conversions of geometry.

It also provides:
 - raycaster
 - fragmentation and decimation features
 - keyholding
 - arraying
 - inversion (including n-polygon forms)
 - sanitization (point sequence orientation, re-ordering, etc.)
 
### GeoCore:

This is a layout parsing library for GDSII and Oasis files. The design is a mix of inspiration from various references. Example usage is available in the root level geoCoreTest and also the gdsTest and oasTest projects.

## Eto

### Eto.VeldridSurface:

This uses the Eto.Veldrid project to deliver a common 2D viewport library for use inside end-user applications. The choice of API depends on platform. This viewport is supported on .NET Core 3.1 and modern .NET, where etoViewport cannot run due to OpenTK issues. The shaders folder is relevant to this project.

### etoViewport:

This is legacy code. It uses Eto.OpenTK to provide common 2D viewports, but has been replaced by Eto.Veldrid due to general incompatibility with .NET Core 3.1 and later.

### errorReporter:

Simple error handler, intended to allow end-user code to call errorReporter in headless or Eto without worrying about how the error is being shown to the user.
