using System.Drawing;
using NetTopologySuite.Index.KdTree;
using Xunit;

namespace HydrologyMaps.Tests;

public class KdTree2DTest
{
    [Fact]
    public void FindExact()
    {
        // Arrange
        KdTree<GraphNode> tree = new KdTree<GraphNode>(node => node.X, node => node.Y);

        tree.Insert(new GraphNode(new Point(1, 2), NodeType.Default));
        tree.Insert(new GraphNode(new Point(3, 5), NodeType.Default));
        tree.Insert(new GraphNode(new Point(6, 1), NodeType.Default));
        tree.Insert(new GraphNode(new Point(9, 3), NodeType.Default));
        tree.Insert(new GraphNode(new Point(4, 7), NodeType.Default));

        // Act
        GraphNode closestNode = tree.FindClosest(6, 1, 1).First();

        // Assert
        Assert.Equal(6, closestNode.X);
        Assert.Equal(1, closestNode.Y);
    }

    [Fact]
    public void FindNearby()
    {
        // Arrange
        KdTree<GraphNode> tree = new KdTree<GraphNode>(node => node.X, node => node.Y);

        tree.Insert(new GraphNode(new Point(1, 2), NodeType.Default));
        tree.Insert(new GraphNode(new Point(3, 5), NodeType.Default));
        tree.Insert(new GraphNode(new Point(6, 1), NodeType.Default));
        tree.Insert(new GraphNode(new Point(9, 3), NodeType.Default));
        tree.Insert(new GraphNode(new Point(4, 7), NodeType.Default));

        // Act
        GraphNode closestNode = tree.FindClosest(5, 4, 1).First();

        // Assert
        Assert.Equal(3, closestNode.X);
        Assert.Equal(5, closestNode.Y);
    }

    [Fact]
    public void FindNearbyFar()
    {
        // Arrange
        KdTree<GraphNode> tree = new KdTree<GraphNode>(node => node.X, node => node.Y);

        tree.Insert(new GraphNode(new Point(1, 2), NodeType.Default));
        tree.Insert(new GraphNode(new Point(3, 5), NodeType.Default));
        tree.Insert(new GraphNode(new Point(6, 1), NodeType.Default));
        tree.Insert(new GraphNode(new Point(9, 3), NodeType.Default));
        tree.Insert(new GraphNode(new Point(4, 7), NodeType.Default));
        tree.Insert(new GraphNode(new Point(6, 7), NodeType.Default));

        // Act
        GraphNode closestNode = tree.FindClosest(5000, 4000, 1).First();

        // Assert
        Assert.Equal(6, closestNode.X);
        Assert.Equal(7, closestNode.Y);
    }

    [Fact]
    public void Find4Nearby()
    {
        // Arrange
        KdTree<GraphNode> tree = new KdTree<GraphNode>(node => node.X, node => node.Y);

        var nonAddedNode1 = new GraphNode(new Point(1, 2), NodeType.Default);
        var nonAddedNode2 = new GraphNode(new Point(-5, 2), NodeType.Default);
        var addedNode1 = new GraphNode(new Point(6, 1), NodeType.Default);

        tree.Insert(nonAddedNode1);
        tree.Insert(nonAddedNode2);
        tree.Insert(addedNode1);
        tree.Insert(new GraphNode(new Point(9, 3), NodeType.Default));
        tree.Insert(new GraphNode(new Point(4, 7), NodeType.Default));
        tree.Insert(new GraphNode(new Point(6, 7), NodeType.Default));

        // Act
        var closestNodes = tree.FindClosest(5000, 4000, 4);
        // Assert
        Assert.False(closestNodes.Contains(nonAddedNode1));
        Assert.False(closestNodes.Contains(nonAddedNode2));
        Assert.True(closestNodes.Contains(addedNode1));
    }
}