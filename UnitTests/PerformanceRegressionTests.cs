using Clipper2Lib;
using geoWrangler;
using NUnit.Framework;

namespace UnitTests
{
    [TestFixture]
    public class PerformanceRegressionTests
    {
        private const double FLOATING_TOLERANCE = 1e-9;

        [Test]
        public void raycaster_numerical_accuracy_test()
        {
            // Test case 1: Simple rectangle
            PathD rect_source = new PathD
            {
                new PointD(0, 0),
                new PointD(100, 0),
                new PointD(100, 100),
                new PointD(0, 100),
                new PointD(0, 0)
            };

            PathD rect_coll = new PathD
            {
                new PointD(-50, -50),
                new PointD(150, -50),
                new PointD(150, 150),
                new PointD(-50, 150),
                new PointD(-50, -50)
            };

            // Create raycaster with deterministic parameters
            RayCast rc = new RayCast(rect_source, rect_coll, 1000, projectCorners: true);
            PathsD clipped_rays = rc.getClippedRays();

            // Expected results for a simple rectangle case - these should be consistent
            Assert.That(clipped_rays.Count, Is.EqualTo(5), "Clipped rays count mismatch");

            // Verify specific point accuracy
            Assert.That(Math.Abs(clipped_rays[0][1].x - (-1000.0)), Is.LessThanOrEqualTo(FLOATING_TOLERANCE), "X coordinate precision loss");
            Assert.That(Math.Abs(clipped_rays[0][1].y - 0.0), Is.LessThanOrEqualTo(FLOATING_TOLERANCE), "Y coordinate precision loss");
        }

        [Test]
        public void raycaster_inversion_mode_accuracy_test()
        {
            PathD rect_source = new PathD
            {
                new PointD(10, 10),
                new PointD(20, 10),
                new PointD(20, 20),
                new PointD(10, 20),
                new PointD(10, 10)
            };

            PathD rect_coll = new PathD
            {
                new PointD(0, 0),
                new PointD(30, 0),
                new PointD(30, 30),
                new PointD(0, 30),
                new PointD(0, 0)
            };

            // Test different inversion modes for consistency
            RayCast rc_none = new RayCast(rect_source, rect_coll, 100, invert: RayCast.inversionMode.none);
            RayCast rc_x = new RayCast(rect_source, rect_coll, 100, invert: RayCast.inversionMode.x);
            RayCast rc_y = new RayCast(rect_source, rect_coll, 100, invert: RayCast.inversionMode.y);

            PathsD rays_none = rc_none.getClippedRays();
            PathsD rays_x = rc_x.getClippedRays();
            PathsD rays_y = rc_y.getClippedRays();

            // Verify that all modes produce consistent output structure
            Assert.That(rays_none.Count, Is.EqualTo(rays_x.Count), "Inversion mode x count mismatch");
            Assert.That(rays_none.Count, Is.EqualTo(rays_y.Count), "Inversion mode y count mismatch");

            // Verify numerical stability between runs
            RayCast rc_none_repeat = new RayCast(rect_source, rect_coll, 100, invert: RayCast.inversionMode.none);
            PathsD rays_none_repeat = rc_none_repeat.getClippedRays();

            Assert.That(rays_none.Count, Is.EqualTo(rays_none_repeat.Count), "Repeated execution count mismatch");
            for (int i = 0; i < rays_none.Count; i++)
            {
                for (int j = 0; j < rays_none[i].Count; j++)
                {
                    Assert.That(Math.Abs(rays_none[i][j].x - rays_none_repeat[i][j].x), Is.LessThanOrEqualTo(FLOATING_TOLERANCE), $"X coordinate differs at ray {i}, point {j}");
                    Assert.That(Math.Abs(rays_none[i][j].y - rays_none_repeat[i][j].y), Is.LessThanOrEqualTo(FLOATING_TOLERANCE), $"Y coordinate differs at ray {i}, point {j}");
                }
            }
        }

        [Test]
        public void decomposition_numerical_accuracy_test()
        {
            // Simple L-shaped polygon for decomposition test
            PathD lShape = new PathD
            {
                new PointD(0, 0),
                new PointD(100, 0),
                new PointD(100, 50),
                new PointD(50, 50),
                new PointD(50, 100),
                new PointD(0, 100),
                new PointD(0, 0)
            };

            PathsD source = new PathsD { lShape };

            // Test basic decomposition
            PathsD decomposed = GeoWrangler.decompose(source);
            
            // Verify that decomposition preserves area
            double originalArea = Math.Abs(Clipper.Area(lShape));
            double decomposedArea = 0.0;
            for (int i = 0; i < decomposed.Count; i++)
            {
                decomposedArea += Math.Abs(Clipper.Area(decomposed[i]));
            }

            Assert.That(Math.Abs(originalArea - decomposedArea), Is.LessThanOrEqualTo(FLOATING_TOLERANCE), "Area not preserved in decomposition");

            // Test consistency across multiple runs
            PathsD decomposed_repeat = GeoWrangler.decompose(source);
            Assert.That(decomposed.Count, Is.EqualTo(decomposed_repeat.Count), "Decomposition count not consistent");

            // Verify numerical consistency
            for (int i = 0; i < decomposed.Count; i++)
            {
                double area1 = Clipper.Area(decomposed[i]);
                double area2 = Clipper.Area(decomposed_repeat[i]);
                Assert.That(Math.Abs(area1 - area2), Is.LessThanOrEqualTo(FLOATING_TOLERANCE), $"Area differs for decomposed polygon {i}");
            }
        }

        [Test]
        public void rectangular_decomposition_numerical_accuracy_test()
        {
            // Create a complex polygon that should decompose into rectangles
            PathD complexPoly = new PathD
            {
                new PointD(0, 0),
                new PointD(60, 0),
                new PointD(60, 20),
                new PointD(40, 20),
                new PointD(40, 40),
                new PointD(60, 40),
                new PointD(60, 60),
                new PointD(0, 60),
                new PointD(0, 0)
            };

            bool abort = false;
            PathsD rectangles = GeoWrangler.rectangular_decomposition(ref abort, complexPoly, 1000, 0, true);

            // Verify that decomposition completes without aborting
            Assert.That(abort, Is.False, "Decomposition was aborted");

            // Verify area preservation
            double originalArea = Math.Abs(Clipper.Area(complexPoly));
            double decomposedArea = 0.0;
            for (int i = 0; i < rectangles.Count; i++)
            {
                decomposedArea += Math.Abs(Clipper.Area(rectangles[i]));
            }

            Assert.That(Math.Abs(originalArea - decomposedArea), Is.LessThanOrEqualTo(FLOATING_TOLERANCE), "Area not preserved in rectangular decomposition");

            // Test with horizontal vs vertical decomposition for consistency
            bool abort2 = false;
            PathsD rectangles_horiz = GeoWrangler.rectangular_decomposition(ref abort2, complexPoly, 1000, 0, false);

            Assert.That(abort2, Is.False, "Horizontal decomposition was aborted");

            double horizArea = 0.0;
            for (int i = 0; i < rectangles_horiz.Count; i++)
            {
                horizArea += Math.Abs(Clipper.Area(rectangles_horiz[i]));
            }

            Assert.That(Math.Abs(originalArea - horizArea), Is.LessThanOrEqualTo(FLOATING_TOLERANCE), "Area not preserved in horizontal decomposition");
        }

        [Test]
        public void raycaster_multisample_accuracy_test()
        {
            PathD rect_source = new PathD
            {
                new PointD(0, 0),
                new PointD(50, 0),
                new PointD(50, 50),
                new PointD(0, 50),
                new PointD(0, 0)
            };

            PathD rect_coll = new PathD
            {
                new PointD(-25, -25),
                new PointD(75, -25),
                new PointD(75, 75),
                new PointD(-25, 75),
                new PointD(-25, -25)
            };

            // Test with different multisample ray counts
            RayCast rc_0 = new RayCast(rect_source, rect_coll, 200, multisampleRayCount: 0);
            RayCast rc_2 = new RayCast(rect_source, rect_coll, 200, multisampleRayCount: 2);
            RayCast rc_4 = new RayCast(rect_source, rect_coll, 200, multisampleRayCount: 4);

            PathsD rays_0 = rc_0.getClippedRays();
            PathsD rays_2 = rc_2.getClippedRays();
            PathsD rays_4 = rc_4.getClippedRays();

            // All should have the same number of emission points
            Assert.That(rays_0.Count, Is.EqualTo(rays_2.Count), "Multisample count mismatch for 2 samples");
            Assert.That(rays_0.Count, Is.EqualTo(rays_4.Count), "Multisample count mismatch for 4 samples");

            // Verify that multisampling produces reasonable results (not drastically different)
            for (int i = 0; i < rays_0.Count; i++)
            {
                Assert.That(rays_0[i].Count, Is.EqualTo(2), "Ray should have start and end point");
                Assert.That(rays_2[i].Count, Is.EqualTo(2), "Multisampled ray should have start and end point");
                Assert.That(rays_4[i].Count, Is.EqualTo(2), "Multisampled ray should have start and end point");

                // Check that multisampled results are within reasonable bounds
                double dist_0 = GeoWrangler.distanceBetweenPoints(rays_0[i][0], rays_0[i][1]);
                double dist_2 = GeoWrangler.distanceBetweenPoints(rays_2[i][0], rays_2[i][1]);
                double dist_4 = GeoWrangler.distanceBetweenPoints(rays_4[i][0], rays_4[i][1]);

                Assert.That(dist_0, Is.GreaterThan(0), "Ray length should be positive");
                Assert.That(dist_2, Is.GreaterThan(0), "Ray length should be positive");
                Assert.That(dist_4, Is.GreaterThan(0), "Ray length should be positive");

                // The multisampled results should not be wildly different from single sample
                Assert.That(Math.Abs(dist_0 - dist_2) / dist_0, Is.LessThan(0.5), "Multisample 2 result too different");
                Assert.That(Math.Abs(dist_0 - dist_4) / dist_0, Is.LessThan(0.5), "Multisample 4 result too different");
            }
        }
    }
}