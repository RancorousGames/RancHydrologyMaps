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

public struct GraphNode : IGraphNode
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

public struct DirectedNode : IGraphNode
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

struct HydrologyParameters
{
    public double ContinueProbability = 0.1;
    public double SymmetricBranchProbability = 0.4;
    public double AsymmetricBranchProbability = 0.5;
    public int SpaceBetwenRiverMouthCandidates = 25;

    public HydrologyParameters()
    {
    }
}

public class HydrologyMap
{
    public HydrologyMap(float[,] heightMap, List<DirectedNode> allRiverNodes, List<RiverEdge> riverEdges,
        List<GraphEdge> voronoiEdges)
    {
        HeightMap = heightMap;
        RiverEdges = riverEdges;
        VoronoiEdges = voronoiEdges;
        AllRiverNodes = allRiverNodes;
    }

    public float[,] HeightMap { get; }
    public List<RiverEdge> RiverEdges { get; }
    public List<GraphEdge> VoronoiEdges { get; }
    public List<DirectedNode> AllRiverNodes { get; }
}