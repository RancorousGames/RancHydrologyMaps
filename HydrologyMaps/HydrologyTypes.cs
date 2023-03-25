using System.Numerics;
using Voronoi2;

namespace HydrologyMaps;

public class Point
{
    public int X { get; set; }
    public int Y { get; set; }

    public Point(int x, int y)
    {
        X = x;
        Y = y;
    }

    public static Point operator +(Point a, Point b)
    {
        return new Point(a.X + b.X, a.Y + b.Y);
    }

    public static Point operator -(Point a, Point b)
    {
        return new Point(a.X - b.X, a.Y - b.Y);
    }

    public static bool operator ==(Point a, Point b)
    {
        if (ReferenceEquals(a, b))
        {
            return true;
        }

        if (a is null || b is null)
        {
            return false;
        }

        return a.X == b.X && a.Y == b.Y;
    }

    public static Point operator -(Point point)
    {
        return new Point(-point.X, -point.Y);
    }
    
    public static bool operator !=(Point a, Point b)
    {
        return !(a == b);
    }
    
    public override bool Equals(object obj)
    {
        if (obj is Point other)
        {
            return X == other.X && Y == other.Y;
        }

        return false;
    }

    public override int GetHashCode()
    {
        return X ^ Y;
    }
    
    public static Point Zero => new Point(0, 0);

    public Point Abs()
    {
        return new Point(Math.Abs(X), Math.Abs(Y));
    }
    
}

public class Vector2D
{
    public double X { get; set; }
    public double Y { get; set; }

    public Vector2D(double x, double y)
    {
        X = x;
        Y = y;
    }

    public Vector2D(Point point)
    {
        X = point.X;
        Y = point.Y;
    }
    
    public Point AsPoint()
    {
        return new Point((int)Math.Round(X), (int)Math.Round(Y));
    }

    public static Vector2D operator +(Vector2D v1, Vector2D v2)
    {
        return new Vector2D(v1.X + v2.X, v1.Y + v2.Y);
    }

    public static Vector2D operator -(Vector2D v1, Vector2D v2)
    {
        return new Vector2D(v1.X - v2.X, v1.Y - v2.Y);
    }

    public static Vector2D operator *(Vector2D v, double scalar)
    {
        return new Vector2D(v.X * scalar, v.Y * scalar);
    }

    public static Vector2D operator *(double scalar, Vector2D v)
    {
        return new Vector2D(v.X * scalar, v.Y * scalar);
    }

    public static Vector2D operator /(Vector2D v, double scalar)
    {
        if (scalar == 0)
        {
            throw new DivideByZeroException("Cannot divide vector by zero.");
        }
        return new Vector2D(v.X / scalar, v.Y / scalar);
    }

    public double Magnitude()
    {
        return Math.Sqrt(X * X + Y * Y);
    }

    public Vector2D Normalized()
    {
        double magnitude = Magnitude();
        if (magnitude == 0)
        {
            return new Vector2D(0, 0);
        }
        return new Vector2D(X / magnitude, Y / magnitude);
    }

    public double Dot(Vector2D other)
    {
        return X * other.X + Y * other.Y;
    }

    public double Cross(Vector2D other)
    {
        return X * other.Y - Y * other.X;
    }

    public static (double distance, bool isLeft) DistanceBetweenPointAndLineSegment(Point start, Point end, Point point)
    {
        return DistanceBetweenPointAndLineSegment(new Vector2D(start), new Vector2D(end), new Vector2D(point));
    }
    
    public static (double distance, bool isLeft) DistanceBetweenPointAndLineSegment(Vector2D start, Vector2D end, Vector2D point)
    {
        // Find the displacement vector of the line segment
        Vector2D displacement = end - start;

        // Find a vector from start to point
        Vector2D startToPoint = point - start;

        // Calculate the dot product of displacement vector and startToPoint vector
        double dotProduct = displacement.Dot(startToPoint);

        // Calculate the cross product of displacement vector and startToPoint vector
        double crossProduct = displacement.Cross(startToPoint);

        // Determine the position of the point relative to the vector
        bool isLeft = crossProduct > 0;

        // If the dot product is negative, the closest point is start
        if (dotProduct < 0)
        {
            return (startToPoint.Magnitude(), isLeft);
        }
        // Calculate the square of the displacement vector's magnitude
        double displacementSquared = displacement.Dot(displacement);

        // If the dot product is greater than the displacement vector's magnitude squared, the closest point is end
        if (dotProduct > displacementSquared)
        {
            return ((point - end).Magnitude(), isLeft);
        }

        // Calculate the projection of startToPoint onto the displacement vector
        Vector2D projection = start + (displacement * (dotProduct / displacementSquared));

        // Calculate the displacement vector from projection to point
        Vector2D projectionToPoint = point - projection;

        // Calculate and return the magnitude of the projectionToPoint vector and position
        return (projectionToPoint.Magnitude(), isLeft);
    }
}

public struct RiverEdge
{
    public RiverEdge(DirectedNode p1, DirectedNode p2)
    {
        P1 = p1;
        P2 = p2;
    }

    public DirectedNode P2 { get; }

    public DirectedNode P1 { get; }
    
    public static double DistanceBetweenPointAndLineSegment(Vector2D start, Vector2D end, Vector2D point)
    {
        // Find the displacement vector of the line segment
        Vector2D displacement = end - start;

        // Find a unit vector n that is perpendicular to the displacement vector
        Vector2D n = new Vector2D(-displacement.Y, displacement.X).Normalized();

        // Find a vector from start to point
        Vector2D startToPoint = point - start;

        // Calculate the dot product of displacement vector and startToPoint vector
        double dotProduct = displacement.Dot(startToPoint);

        // If the dot product is negative, the closest point is start, otherwise check if it's greater than the displacement vector
        if (dotProduct < 0)
        {
            return startToPoint.Magnitude();
        }
        else if (dotProduct > displacement.Dot(displacement))
        {
            return (point - end).Magnitude();
        }

        // Calculate the projection of startToPoint onto the displacement vector
        Vector2D projection = start + (displacement * (dotProduct / displacement.Dot(displacement)));

        // Calculate the displacement vector from projection to point
        Vector2D projectionToPoint = point - projection;

        // Calculate and return the magnitude of the cross product of the projectionToPoint vector and the unit vector n
        return Math.Abs(projectionToPoint.Cross(n)) / n.Magnitude();
    }
}

public enum NodeType
{
    Default,
    Border,
    River,
    RiverMouth
}

public interface IGraphNode
{
    Point Point { get; set; }
    NodeType Type { get; set; }
}

public class GraphNode : IGraphNode
{
    public int X
    {
        get => Point.X;
        set => Point = new Point(value, Point.Y);
    }

    public int Y
    {
        get => Point.Y;
        set => Point = new Point(Point.X, value);
    }

    public Point Point { get; set; }
    public NodeType Type { get; set; }

    public GraphNode(Point point, NodeType type)
    {
        Point = point;
        Type = type;
        Position = new double[] { point.X, point.Y };
    }

    public double[] Position { get; set; }
}

public class DirectedNode : IGraphNode
{
    public int X
    {
        get => Point.X;
        set => Point = new Point(value, Point.Y);
    }

    public int Y
    {
        get => Point.Y;
        set => Point = new Point(Point.X, value);
    }

    public Point Point { get; set; }
    public NodeType Type { get; set; }
    public Vector2 Direction { get; set; }

    public DirectedNode(Point point, NodeType type, Vector2 direction)
    {
        Point = point;
        Type = type;
        Direction = direction;
        Position = new double[] { point.X, point.Y };
    }

    public double[] Position { get; set; }
}

enum RiverAction
{
    Continue,
    SymmetricBranch,
    AsymmetricBranch
}

public struct HydrologyParameters
{
    public double ContinueProbability = 0.1;
    public double SymmetricBranchProbability = 0.4;
    public double AsymmetricBranchProbability = 0.5;
    public int SpaceBetweenRiverMouthCandidates = 25;
    
    public HydrologyParameters()
    {
    }
}

public class HydrologyMap
{
    public HydrologyMap(float[,] heightMap, List<DirectedNode> allRiverNodes, List<RiverEdge> riverEdges,
        List<GraphEdge> voronoiEdges, List<IGraphNode> voronoiNodes)
    {
        HeightMap = heightMap;
        RiverEdges = riverEdges;
        VoronoiEdges = voronoiEdges;
        VoronoiNodes = voronoiNodes;
        AllRiverNodes = allRiverNodes;
    }

    public float[,] HeightMap { get; }
    public List<RiverEdge> RiverEdges { get; }
    public List<GraphEdge> VoronoiEdges { get; }
    public List<IGraphNode> VoronoiNodes { get; }
    public List<DirectedNode> AllRiverNodes { get; }
}