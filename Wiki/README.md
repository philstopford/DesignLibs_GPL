# DesignLibs_GPL API Documentation

This directory contains comprehensive API documentation for all major DesignLibs_GPL libraries.

## 📚 Documentation Files

| File | Library | Description |
|------|---------|-------------|
| `Home.md` | Overview | Main navigation hub and library overview |
| `geoLib-API.md` | geoLib | Basic geometric primitives (points, arrays, matrices, rectangles) |
| `geoCore-API.md` | geoCore | GDSII and Oasis file format parsing and writing |
| `geoWrangler-API.md` | geoWrangler | Geometry transformations, raycasting, fragmentation, and decimation |
| `shapeEngine-API.md` | shapeEngine | Shape generation and manipulation (rectangles, L-shapes, T-shapes, etc.) |
| `clipper-API.md` | clipper | Advanced polygon clipping operations (Clipper2 C# port) |
| `utility-API.md` | utility | General purpose utilities (compression, hashing, mathematical operations) |
| `MersenneTwister-API.md` | MersenneTwister | High-quality pseudorandom number generation |
| `Noise-API.md` | Noise | Perlin, Simplex, and OpenSimplex noise functions |
| `Eto.VeldridSurface-API.md` | Eto.VeldridSurface | Cross-platform 2D viewport for end-user applications |

## 🚀 Deploying to GitHub Wiki

The documentation in this directory is **NOT automatically published** to the GitHub wiki. You need to manually deploy it.

### Quick Deploy (Recommended)
Use the provided automated script:
```bash
./deploy-wiki.sh
```

### Manual Deploy
See `WIKI_DEPLOYMENT_GUIDE.md` for detailed manual deployment instructions.

## 📊 Documentation Statistics

- **10 comprehensive API documentation files**
- **5,529+ total lines** of detailed documentation
- **Complete coverage** of all major DesignLibs_GPL libraries
- **Professional-quality** API reference material
- **Ready-to-deploy** content formatted for GitHub wiki

## 🔄 Updating Documentation

To update the wiki after making changes:
1. Edit the `.md` files in this directory
2. Run `./deploy-wiki.sh` to deploy changes
3. The GitHub wiki will be updated automatically

## 📖 Content Features

Each library documentation includes:
- **Complete API Reference** - All classes, methods, properties, and enumerations
- **Comprehensive Examples** - Real-world usage patterns and code samples  
- **Integration Guidance** - How libraries work together in the DesignLibs ecosystem
- **Performance Tips** - Optimization guidelines and best practices
- **Error Handling** - Proper validation and error management patterns
- **Cross-References** - Links between related libraries showing their interactions

This documentation provides everything needed to effectively use the DesignLibs_GPL library ecosystem.