using NetTopologySuite.Utilities;

namespace HydrologyMaps;

public class HydrologyKdTree2D
{
    private class KdNode
    {
        public GraphNode Point;
        public KdNode Left;
        public KdNode Right;

        public KdNode(GraphNode point)
        {
            Point = point;
            Left = null;
            Right = null;
        }
    }
    
    private class ComparableGraphNode : IComparable<ComparableGraphNode>
    {
        public GraphNode Node { get; }
        public int Distance { get; set; }

        public ComparableGraphNode(GraphNode node, int distance)
        {
            Node = node;
            Distance = distance;
        }

        public int CompareTo(ComparableGraphNode other)
        {
            return Distance.CompareTo(other.Distance);
        }
    }

    private KdNode root;

    public HydrologyKdTree2D()
    {
        root = null;
    }

    public void Insert(GraphNode point)
    {
        root = InsertRecursively(root, point, 0);
    }

    private KdNode InsertRecursively(KdNode currentNode, GraphNode point, int depth)
    {
        if (currentNode == null)
        {
            return new KdNode(point);
        }

        int currentDimension = depth % 2;
        if (currentDimension == 0)
        {
            if (point.X < currentNode.Point.X)
            {
                currentNode.Left = InsertRecursively(currentNode.Left, point, depth + 1);
            }
            else
            {
                currentNode.Right = InsertRecursively(currentNode.Right, point, depth + 1);
            }
        }
        else
        {
            if (point.Y < currentNode.Point.Y)
            {
                currentNode.Left = InsertRecursively(currentNode.Left, point, depth + 1);
            }
            else
            {
                currentNode.Right = InsertRecursively(currentNode.Right, point, depth + 1);
            }
        }

        return currentNode;
    }

    public List<GraphNode> FindClosest(int x, int y, int n)
    {
        var nearestNodes = new PriorityQueue<ComparableGraphNode>(n, Comparer<ComparableGraphNode>.Create((a, b) => b.Distance - a.Distance));
        FindClosestRecursively(root, x, y, 0, nearestNodes, n);
        return nearestNodes.Select(item => item.Node).ToList();
    }
    
  // public GraphNode FindClosest(int x, int y)
  // {
  //     return FindClosestRecursively(root, x, y, 0).Point;
  // }

  
    private void FindClosestRecursively(KdNode currentNode, int X, int Y, int depth, PriorityQueue<ComparableGraphNode> nearestNodes, int n)
    {
        if (currentNode == null)
        {
            return;
        }

        int currentDimension = depth % 2;
        KdNode nextNode = null;
        KdNode otherNode = null;

        if ((currentDimension == 0 && X < currentNode.Point.X) || (currentDimension == 1 && Y < currentNode.Point.Y))
        {
            nextNode = currentNode.Left;
            otherNode = currentNode.Right;
        }
        else
        {
            nextNode = currentNode.Right;
            otherNode = currentNode.Left;
        }

        FindClosestRecursively(nextNode, X, Y, depth + 1, nearestNodes, n);

        int distToCurrentNode = DistanceSquared(currentNode.Point, X, Y);
        if (nearestNodes.Count() < n || distToCurrentNode < nearestNodes.Peek().Distance)
        {
            if (nearestNodes.Count() == n)
            {
                nearestNodes.Poll();
            }
            nearestNodes.Add(new ComparableGraphNode(currentNode.Point, distToCurrentNode));
        }

        int distToSplittingPlane = 0;
        if (currentDimension == 0)
        {
            distToSplittingPlane = X - currentNode.Point.X;
        }
        else
        {
            distToSplittingPlane = Y - currentNode.Point.Y;
        }

        if (nearestNodes.Count() < n || distToSplittingPlane * distToSplittingPlane < DistanceSquared(nearestNodes.Peek().Node, X, Y))
        {
            FindClosestRecursively(otherNode, X, Y, depth + 1, nearestNodes, n);
        }
    }

    private int DistanceSquared(GraphNode a, int X, int Y)
    {
        int deltaX = a.X - X;
        int deltaY = a.Y - Y;
        return deltaX * deltaX + deltaY * deltaY;
    }
}