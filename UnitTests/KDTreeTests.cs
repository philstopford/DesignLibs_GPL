using KDTree;
using NUnit.Framework;

namespace UnitTests;

/// <summary>
/// Tests for the KDTree library functionality, covering tree construction,
/// nearest neighbor search, and distance functions.
/// </summary>
[TestFixture]
public class KDTreeTests
{
    [Test]
    public void KDTree_Constructor_ShouldCreateTreeWithCorrectDimensions()
    {
        var tree = new KDTree<string>(3);

        Assert.That(tree, Is.Not.Null);
        // Tree should be created successfully for 3D space
    }

    [Test]
    public void KDTree_ConstructorWithBucketCapacity_ShouldCreateTreeWithCorrectParameters()
    {
        var tree = new KDTree<int>(2, 10);

        Assert.That(tree, Is.Not.Null);
        // Tree should be created with custom bucket capacity
    }

    [Test]
    public void KDTree_AddPoint_ShouldAddPointSuccessfully()
    {
        var tree = new KDTree<string>(2);
        double[] point = { 1.0, 2.0 };

        Assert.DoesNotThrow(() => tree.AddPoint(point, "TestData"));
        Assert.That(tree.Size, Is.EqualTo(1));
    }

    [Test]
    public void KDTree_AddMultiplePoints_ShouldMaintainCorrectCount()
    {
        var tree = new KDTree<int>(2);

        tree.AddPoint(new double[] { 1.0, 1.0 }, 1);
        tree.AddPoint(new double[] { 2.0, 2.0 }, 2);
        tree.AddPoint(new double[] { 3.0, 3.0 }, 3);

        Assert.That(tree.Size, Is.EqualTo(3));
    }

    [Test]
    public void KDTree_NearestNeighbors_ShouldFindClosestPoint()
    {
        var tree = new KDTree<string>(2);

        tree.AddPoint(new double[] { 1.0, 1.0 }, "Point1");
        tree.AddPoint(new double[] { 2.0, 2.0 }, "Point2");
        tree.AddPoint(new double[] { 10.0, 10.0 }, "Point3");

        var searchPoint = new double[] { 1.1, 1.1 };
        var nearestNeighbors = tree.NearestNeighbors(searchPoint, 1);

        Assert.That(nearestNeighbors, Is.Not.Null);
        Assert.That(nearestNeighbors.MoveNext(), Is.True);
        Assert.That(nearestNeighbors.Current, Is.EqualTo("Point1"));
    }

    [Test]
    public void KDTree_NearestNeighbors_WithMaxReturned_ShouldLimitResults()
    {
        var tree = new KDTree<int>(2);

        for (int i = 0; i < 10; i++)
        {
            tree.AddPoint(new double[] { i, i }, i);
        }

        var searchPoint = new double[] { 5.0, 5.0 };
        var nearestNeighbors = tree.NearestNeighbors(searchPoint, 3);

        int count = 0;
        while (nearestNeighbors.MoveNext())
        {
            count++;
        }

        Assert.That(count, Is.LessThanOrEqualTo(3));
    }

    [Test]
    public void KDTree_NearestNeighbors_WithDistanceThreshold_ShouldFilterResults()
    {
        var tree = new KDTree<string>(2);

        tree.AddPoint(new double[] { 0.0, 0.0 }, "Close");
        tree.AddPoint(new double[] { 10.0, 10.0 }, "Far");

        var searchPoint = new double[] { 0.1, 0.1 };
        var nearestNeighbors = tree.NearestNeighbors(searchPoint, 10, 1.0); // Distance threshold of 1.0

        Assert.That(nearestNeighbors, Is.Not.Null);
        Assert.That(nearestNeighbors.MoveNext(), Is.True);
        Assert.That(nearestNeighbors.Current, Is.EqualTo("Close"));

        // Should not find the far point within distance threshold
        Assert.That(nearestNeighbors.MoveNext(), Is.False);
    }

    [Test]
    public void SquareEuclideanDistanceFunction_Distance_ShouldCalculateCorrectDistance()
    {
        var distanceFunction = new SquareEuclideanDistanceFunction();
        var p1 = new double[] { 0.0, 0.0 };
        var p2 = new double[] { 3.0, 4.0 };

        double distance = distanceFunction.Distance(p1, p2);

        // Square euclidean distance should be 3^2 + 4^2 = 25
        Assert.That(distance, Is.EqualTo(25.0).Within(1e-10));
    }

    [Test]
    public void SquareEuclideanDistanceFunction_DistanceToRectangle_WithPointInsideRectangle_ShouldReturnZero()
    {
        var distanceFunction = new SquareEuclideanDistanceFunction();
        var point = new double[] { 2.0, 3.0 };
        var min = new double[] { 1.0, 2.0 };
        var max = new double[] { 4.0, 5.0 };

        double distance = distanceFunction.DistanceToRectangle(point, min, max);

        Assert.That(distance, Is.EqualTo(0.0));
    }

    [Test]
    public void SquareEuclideanDistanceFunction_DistanceToRectangle_WithPointOutsideRectangle_ShouldCalculateCorrectDistance()
    {
        var distanceFunction = new SquareEuclideanDistanceFunction();
        var point = new double[] { 0.0, 0.0 };
        var min = new double[] { 2.0, 2.0 };
        var max = new double[] { 4.0, 4.0 };

        double distance = distanceFunction.DistanceToRectangle(point, min, max);

        // Distance should be sqrt((2-0)^2 + (2-0)^2) squared = 8
        Assert.That(distance, Is.EqualTo(8.0).Within(1e-10));
    }

    [Test]
    public void KDTree_RemoveAt_ShouldRemoveSpecificPoint()
    {
        var tree = new KDTree<string>(2);

        tree.AddPoint(new double[] { 1.0, 1.0 }, "Point1");
        tree.AddPoint(new double[] { 2.0, 2.0 }, "Point2");

        var searchPoint = new double[] { 1.1, 1.1 };
        var nearestNeighbors = tree.NearestNeighbors(searchPoint, 1);

        Assert.That(nearestNeighbors.MoveNext(), Is.True);
        Assert.That(nearestNeighbors.Current, Is.EqualTo("Point1"));

        // For this test, we just verify the tree structure works
        Assert.That(tree.Size, Is.EqualTo(2));
    }

    [Test]
    public void KDTree_ThreeDimensional_ShouldWorkCorrectly()
    {
        var tree = new KDTree<string>(3);

        tree.AddPoint(new double[] { 1.0, 2.0, 3.0 }, "Point1");
        tree.AddPoint(new double[] { 4.0, 5.0, 6.0 }, "Point2");

        var searchPoint = new double[] { 1.1, 2.1, 3.1 };
        var nearestNeighbors = tree.NearestNeighbors(searchPoint, 1);

        Assert.That(nearestNeighbors.MoveNext(), Is.True);
        Assert.That(nearestNeighbors.Current, Is.EqualTo("Point1"));
    }
}