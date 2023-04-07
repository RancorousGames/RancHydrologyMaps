using NetTopologySuite.Utilities;

namespace HydrologyMaps;

public class KdTree<T>
{
    private class KdNode
    {
        public T Point;
        public KdNode Left;
        public KdNode Right;

        public KdNode(T point)
        {
            Point = point;
            Left = null;
            Right = null;
        }
    }

    private class ComparableItem : IComparable<ComparableItem>
    {
        public T Item { get; }
        public int Distance { get; set; }

        public ComparableItem(T item, int distance)
        {
            Item = item;
            Distance = distance;
        }

        public int CompareTo(ComparableItem other)
        {
            return Distance.CompareTo(other.Distance);
        }
    }

    private KdNode root;
    private readonly Func<T, int> _getX;
    private readonly Func<T, int> _getY;

    public KdTree(Func<T, int> getX, Func<T, int> getY)
    {
        root = null;
        _getX = getX;
        _getY = getY;
    }

    public void Insert(T point)
    {
        root = InsertRecursively(root, point, 0);
    }

    private KdNode InsertRecursively(KdNode currentNode, T point, int depth)
    {
        if (currentNode == null)
        {
            return new KdNode(point);
        }

        int currentDimension = depth % 2;
        if (currentDimension == 0)
        {
            if (_getX(point) < _getX(currentNode.Point))
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
            if (_getY(point) < _getY(currentNode.Point))
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

    private int DistanceSquared(T a, int X, int Y)
    {
        int deltaX = _getX(a) - X;
        int deltaY = _getY(a) - Y;
        return deltaX * deltaX + deltaY * deltaY;
    }

    public List<T> FindClosest(int x, int y, int n)
    {
        var nearestNodes =
            new PriorityQueue<ComparableItem>(n,
                Comparer<ComparableItem>.Create((a, b) => b.Distance - a.Distance));
        FindClosestRecursively(root, x, y, 0, nearestNodes, n);
        return nearestNodes.Select(item => item.Item).ToList();
    }

// public T FindClosest(int x, int y)
// {
//     return FindClosestRecursively(root, x, y, 0).Point;
// }


    private void FindClosestRecursively(KdNode currentNode, int X, int Y, int depth,
        PriorityQueue<ComparableItem> nearestNodes, int n)
    {
        if (currentNode == null)
        {
            return;
        }

        int currentDimension = depth % 2;
        KdNode nextNode = null;
        KdNode otherNode = null;

        if ((currentDimension == 0 && X < _getX(currentNode.Point)) || (currentDimension == 1 && Y < _getY(currentNode.Point)))
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

            nearestNodes.Add(new ComparableItem(currentNode.Point, distToCurrentNode));
        }

        int distToSplittingPlane = 0;
        if (currentDimension == 0)
        {
            distToSplittingPlane = X -  _getX(currentNode.Point);
        }
        else
        {
            distToSplittingPlane = Y -  _getY(currentNode.Point);
        }

        if (nearestNodes.Count() < n ||
            distToSplittingPlane * distToSplittingPlane < DistanceSquared(nearestNodes.Peek().Item, X, Y))
        {
            FindClosestRecursively(otherNode, X, Y, depth + 1, nearestNodes, n);
        }
    }
}