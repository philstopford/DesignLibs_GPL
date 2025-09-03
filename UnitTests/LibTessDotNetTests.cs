using NUnit.Framework;
using LibTessDotNet.Double;
using System;
using System.Collections.Generic;
using System.Linq;

namespace UnitTests;

[TestFixture]
public class LibTessDotNetTests
{
    [Test]
    public void Tess_Constructor_ShouldCreateValidInstance()
    {
        var tess = new Tess();
        Assert.That(tess, Is.Not.Null);
    }

    [Test]
    public void ContourVertex_Constructor_ShouldCreateValidVertex()
    {
        var vertex = new ContourVertex(new Vec3(1.0, 2.0, 0.0));

        Assert.That(vertex.Position.X, Is.EqualTo(1.0));
        Assert.That(vertex.Position.Y, Is.EqualTo(2.0));
        Assert.That(vertex.Position.Z, Is.EqualTo(0.0));
        Assert.That(vertex.Data, Is.Null);
    }

    [Test]
    public void ContourVertex_ConstructorWithData_ShouldStoreData()
    {
        var data = "test_data";
        var vertex = new ContourVertex(new Vec3(1.0, 2.0, 0.0), data);

        Assert.That(vertex.Data, Is.EqualTo(data));
    }

    [Test]
    public void AddContour_SimpleTriangle_ShouldAcceptContour()
    {
        var tess = new Tess();
        var vertices = new ContourVertex[]
        {
            new ContourVertex(new Vec3(0.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(1.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(0.5, 1.0, 0.0))
        };

        Assert.DoesNotThrow(() => tess.AddContour(vertices));
    }

    [Test]
    public void AddContour_WithList_ShouldAcceptContour()
    {
        var tess = new Tess();
        var vertices = new List<ContourVertex>
        {
            new ContourVertex(new Vec3(0.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(1.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(0.5, 1.0, 0.0))
        };

        Assert.DoesNotThrow(() => tess.AddContour(vertices));
    }

    [Test]
    public void Tessellate_SimpleTriangle_ShouldProduceValidResult()
    {
        var tess = new Tess();
        var vertices = new ContourVertex[]
        {
            new ContourVertex(new Vec3(0.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(1.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(0.5, 1.0, 0.0))
        };

        tess.AddContour(vertices);
        tess.Tessellate(WindingRule.EvenOdd, ElementType.Polygons, 3);

        Assert.That(tess.Vertices, Is.Not.Null);
        Assert.That(tess.Vertices.Length, Is.GreaterThan(0));
        Assert.That(tess.Elements, Is.Not.Null);
        Assert.That(tess.Elements.Length, Is.GreaterThan(0));
    }

    [Test]
    public void Tessellate_Square_ShouldProduceTwoTriangles()
    {
        var tess = new Tess();
        var vertices = new ContourVertex[]
        {
            new ContourVertex(new Vec3(0.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(1.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(1.0, 1.0, 0.0)),
            new ContourVertex(new Vec3(0.0, 1.0, 0.0))
        };

        tess.AddContour(vertices);
        tess.Tessellate(WindingRule.EvenOdd, ElementType.Polygons, 3);

        Assert.That(tess.Vertices.Length, Is.EqualTo(4));
        Assert.That(tess.Elements.Length, Is.EqualTo(6)); // 2 triangles = 6 indices
    }

    [Test]
    public void Tessellate_WithDifferentWindingRules_ShouldWork()
    {
        var tess = new Tess();
        var vertices = new ContourVertex[]
        {
            new ContourVertex(new Vec3(0.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(1.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(1.0, 1.0, 0.0)),
            new ContourVertex(new Vec3(0.0, 1.0, 0.0))
        };

        tess.AddContour(vertices);

        // Test different winding rules
        Assert.DoesNotThrow(() => tess.Tessellate(WindingRule.EvenOdd, ElementType.Polygons, 3));
        Assert.DoesNotThrow(() => tess.Tessellate(WindingRule.NonZero, ElementType.Polygons, 3));
        Assert.DoesNotThrow(() => tess.Tessellate(WindingRule.Positive, ElementType.Polygons, 3));
        Assert.DoesNotThrow(() => tess.Tessellate(WindingRule.Negative, ElementType.Polygons, 3));
        Assert.DoesNotThrow(() => tess.Tessellate(WindingRule.AbsGeqTwo, ElementType.Polygons, 3));
    }

    [Test]
    public void Tessellate_WithDifferentElementTypes_ShouldWork()
    {
        var tess = new Tess();
        var vertices = new ContourVertex[]
        {
            new ContourVertex(new Vec3(0.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(1.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(1.0, 1.0, 0.0)),
            new ContourVertex(new Vec3(0.0, 1.0, 0.0))
        };

        tess.AddContour(vertices);

        // Test different element types
        Assert.DoesNotThrow(() => tess.Tessellate(WindingRule.EvenOdd, ElementType.Polygons, 3));
        Assert.DoesNotThrow(() => tess.Tessellate(WindingRule.EvenOdd, ElementType.ConnectedPolygons, 3));
        Assert.DoesNotThrow(() => tess.Tessellate(WindingRule.EvenOdd, ElementType.BoundaryContours));
    }

    [Test]
    public void Tessellate_EmptyContour_ShouldProduceEmptyResult()
    {
        var tess = new Tess();
        tess.Tessellate(WindingRule.EvenOdd, ElementType.Polygons, 3);

        // LibTessDotNet may return null arrays when there's no input
        Assert.That(tess.Vertices?.Length ?? 0, Is.EqualTo(0));
        Assert.That(tess.Elements?.Length ?? 0, Is.EqualTo(0));
    }

    [Test]
    public void Tessellate_ComplexPolygon_ShouldProduceValidResult()
    {
        var tess = new Tess();

        // Create a polygon with a hole
        var outerVertices = new ContourVertex[]
        {
            new ContourVertex(new Vec3(0.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(2.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(2.0, 2.0, 0.0)),
            new ContourVertex(new Vec3(0.0, 2.0, 0.0))
        };

        var holeVertices = new ContourVertex[]
        {
            new ContourVertex(new Vec3(0.5, 0.5, 0.0)),
            new ContourVertex(new Vec3(0.5, 1.5, 0.0)),
            new ContourVertex(new Vec3(1.5, 1.5, 0.0)),
            new ContourVertex(new Vec3(1.5, 0.5, 0.0))
        };

        tess.AddContour(outerVertices, ContourOrientation.Original);
        tess.AddContour(holeVertices, ContourOrientation.Clockwise); // Hole should be clockwise
        tess.Tessellate(WindingRule.EvenOdd, ElementType.Polygons, 3);

        Assert.That(tess.Vertices.Length, Is.GreaterThan(0));
        Assert.That(tess.Elements.Length, Is.GreaterThan(0));
        Assert.That(tess.Elements.Length % 3, Is.EqualTo(0)); // Should be triangles
    }

    [Test]
    public void Tessellate_StarShape_ShouldProduceValidResult()
    {
        var tess = new Tess();
        var vertices = new List<ContourVertex>();

        // Create a 5-pointed star
        var numPoints = 5;
        var outerRadius = 1.0;
        var innerRadius = 0.4;

        for (int i = 0; i < numPoints * 2; i++)
        {
            var angle = i * Math.PI / numPoints;
            var radius = (i % 2 == 0) ? outerRadius : innerRadius;
            var x = Math.Cos(angle) * radius;
            var y = Math.Sin(angle) * radius;
            vertices.Add(new ContourVertex(new Vec3(x, y, 0.0)));
        }

        tess.AddContour(vertices);
        tess.Tessellate(WindingRule.EvenOdd, ElementType.Polygons, 3);

        Assert.That(tess.Vertices.Length, Is.GreaterThan(0));
        Assert.That(tess.Elements.Length, Is.GreaterThan(0));
        Assert.That(tess.Elements.Length % 3, Is.EqualTo(0));
    }

    [Test]
    public void Vec3_Operations_ShouldWorkCorrectly()
    {
        var v1 = new Vec3(1.0, 2.0, 3.0);
        var v2 = new Vec3(4.0, 5.0, 6.0);

        Assert.That(v1.X, Is.EqualTo(1.0));
        Assert.That(v1.Y, Is.EqualTo(2.0));
        Assert.That(v1.Z, Is.EqualTo(3.0));

        // Test ToString
        var str = v1.ToString();
        Assert.That(str, Is.Not.Null);
        Assert.That(str, Contains.Substring("1"));
        Assert.That(str, Contains.Substring("2"));
        Assert.That(str, Contains.Substring("3"));
    }

    [Test]
    public void Tessellate_WithCustomPoolFactory_ShouldWork()
    {
        var pool = new DefaultPool();
        var tess = new Tess(pool);

        var vertices = new ContourVertex[]
        {
            new ContourVertex(new Vec3(0.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(1.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(0.5, 1.0, 0.0))
        };

        tess.AddContour(vertices);
        tess.Tessellate(WindingRule.EvenOdd, ElementType.Polygons, 3);

        Assert.That(tess.Vertices.Length, Is.GreaterThan(0));
        Assert.That(tess.Elements.Length, Is.GreaterThan(0));
    }

    [Test]
    public void Tessellate_BoundaryContours_ShouldProduceBoundaries()
    {
        var tess = new Tess();
        var vertices = new ContourVertex[]
        {
            new ContourVertex(new Vec3(0.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(1.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(1.0, 1.0, 0.0)),
            new ContourVertex(new Vec3(0.0, 1.0, 0.0))
        };

        tess.AddContour(vertices);
        tess.Tessellate(WindingRule.EvenOdd, ElementType.BoundaryContours);

        Assert.That(tess.Vertices.Length, Is.GreaterThan(0));
        Assert.That(tess.Elements.Length, Is.GreaterThan(0));
    }

    [Test]
    public void Tessellate_MultipleContours_ShouldHandleCorrectly()
    {
        var tess = new Tess();

        // Add first contour
        var vertices1 = new ContourVertex[]
        {
            new ContourVertex(new Vec3(0.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(1.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(1.0, 1.0, 0.0)),
            new ContourVertex(new Vec3(0.0, 1.0, 0.0))
        };

        // Add second contour (separate square)
        var vertices2 = new ContourVertex[]
        {
            new ContourVertex(new Vec3(2.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(3.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(3.0, 1.0, 0.0)),
            new ContourVertex(new Vec3(2.0, 1.0, 0.0))
        };

        tess.AddContour(vertices1);
        tess.AddContour(vertices2);
        tess.Tessellate(WindingRule.EvenOdd, ElementType.Polygons, 3);

        Assert.That(tess.Vertices.Length, Is.GreaterThan(0));
        Assert.That(tess.Elements.Length, Is.GreaterThan(0));
        Assert.That(tess.Elements.Length % 3, Is.EqualTo(0));
    }

    [Test]
    public void ContourOrientation_ShouldAffectTessellation()
    {
        var tess = new Tess();
        var vertices = new ContourVertex[]
        {
            new ContourVertex(new Vec3(0.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(1.0, 0.0, 0.0)),
            new ContourVertex(new Vec3(1.0, 1.0, 0.0)),
            new ContourVertex(new Vec3(0.0, 1.0, 0.0))
        };

        // Test different orientations
        Assert.DoesNotThrow(() =>
        {
            tess.AddContour(vertices, ContourOrientation.Original);
            tess.Tessellate(WindingRule.EvenOdd, ElementType.Polygons, 3);
        });

        var tess2 = new Tess();
        Assert.DoesNotThrow(() =>
        {
            tess2.AddContour(vertices, ContourOrientation.Clockwise);
            tess2.Tessellate(WindingRule.EvenOdd, ElementType.Polygons, 3);
        });

        var tess3 = new Tess();
        Assert.DoesNotThrow(() =>
        {
            tess3.AddContour(vertices, ContourOrientation.CounterClockwise);
            tess3.Tessellate(WindingRule.EvenOdd, ElementType.Polygons, 3);
        });
    }
}