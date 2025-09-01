using NUnit.Framework;
using Clipper2Lib;
using shapeEngine;
using System;

namespace UnitTests
{
    [TestFixture]
    public class ContourGenTests
    {
        [Test]
        public void makeContour_WithNullPath_ReturnsEmptyPath()
        {
            // Act
            PathD result = contourGen.makeContour(null, 5.0, 5.0, 1.0, 5.0, 0.1, 1.0);
            
            // Assert
            Assert.IsNotNull(result);
            Assert.AreEqual(0, result.Count);
        }

        [Test]
        public void makeContour_WithEmptyPath_ReturnsEmptyPath()
        {
            // Arrange
            PathD emptyPath = new PathD();
            
            // Act
            PathD result = contourGen.makeContour(emptyPath, 5.0, 5.0, 1.0, 5.0, 0.1, 1.0);
            
            // Assert
            Assert.IsNotNull(result);
            Assert.AreEqual(0, result.Count);
        }

        [Test]
        public void makeContour_WithSinglePoint_ReturnsOriginalPoint()
        {
            // Arrange
            PathD singlePointPath = new PathD();
            singlePointPath.Add(new PointD(5, 10));
            
            // Act
            PathD result = contourGen.makeContour(singlePointPath, 5.0, 5.0, 1.0, 5.0, 0.1, 1.0);
            
            // Assert
            Assert.IsNotNull(result);
            Assert.AreEqual(1, result.Count);
            Assert.AreEqual(5, result[0].x);
            Assert.AreEqual(10, result[0].y);
        }

        [Test]
        public void makeContour_WithTwoPoints_ReturnsOriginalPath()
        {
            // Arrange - This is the original issue: 2-point path (0,0) to (0,60)
            PathD twoPointPath = new PathD();
            twoPointPath.Add(new PointD(0, 0));
            twoPointPath.Add(new PointD(0, 60));
            
            // Act
            PathD result = contourGen.makeContour(twoPointPath, 5.0, 5.0, 1.0, 5.0, 0.1, 1.0);
            
            // Assert
            Assert.IsNotNull(result);
            Assert.AreEqual(2, result.Count);
            Assert.AreEqual(0, result[0].x);
            Assert.AreEqual(0, result[0].y);
            Assert.AreEqual(0, result[1].x);
            Assert.AreEqual(60, result[1].y);
        }

        [Test]
        public void makeContour_WithValidTriangle_ProcessesNormally()
        {
            // Arrange
            PathD trianglePath = new PathD();
            trianglePath.Add(new PointD(0, 0));
            trianglePath.Add(new PointD(10, 0));
            trianglePath.Add(new PointD(5, 10));
            trianglePath.Add(new PointD(0, 0)); // Close the path
            
            // Act
            PathD result = contourGen.makeContour(trianglePath, 2.0, 2.0, 1.0, 5.0, 0.1, 1.0);
            
            // Assert
            Assert.IsNotNull(result);
            Assert.Greater(result.Count, 0); // Should return some contoured points
            // For a valid path, the contour generation should produce a result
        }
    }
}