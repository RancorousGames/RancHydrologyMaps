using System.Drawing;
using System.Drawing.Imaging;
using System.Numerics;
using Voronoi2;
using Point = System.Drawing.Point;

namespace HydrologyMaps;

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
    public int SpaceBetwenRiverMouthCandidates = 30;

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

public class HydrologyMapGen
{
    Random Random = new Random();
    HydrologyParameters parameters = new HydrologyParameters();


    public HydrologyMap GenerateIsland(int width, int height)
    {
        int seed = Random.Next();
        Random = new Random(seed);
        Console.WriteLine("Seed is " + seed);

        FastNoiseLite noiseGenerator = new FastNoiseLite();
        noiseGenerator.SetNoiseType(FastNoiseLite.NoiseType.OpenSimplex2);
        noiseGenerator.SetSeed(Random.Next());
        noiseGenerator.SetFrequency(0.003f);
        noiseGenerator.SetFractalGain(0);

        float[,] heightmap = new float[width, height];

        float centerX = width / 2f;
        float centerY = height / 2f;

        for (int y = 0; y < height; y++)
        {
            for (int x = 0; x < width; x++)
            {
                float distanceX = x - centerX;
                float distanceY = y - centerY;
                float distanceToCenter = (float)Math.Sqrt(distanceX * distanceX + distanceY * distanceY);

                float radialGradient = 1f - Math.Min(1f, distanceToCenter / (Math.Min(width, height) / 2.5f));

                float perlinValue =
                    (noiseGenerator.GetNoise(x, y) + 1f) * radialGradient /
                    2f; // Normalize Perlin noise to the range [0, 1]
                if (perlinValue < 0.15) perlinValue = 0;
                else if (perlinValue < 0.2) perlinValue = 0.15f;
                //else perlinValue = 0.5f;
                heightmap[x, y] = perlinValue;
            }
        }

        List<DirectedNode> riverMouthCandidates = GetRiverMouthCandidates(heightmap, parameters.SpaceBetwenRiverMouthCandidates);
        
        if (riverMouthCandidates.Count == 0)
        {
            return new HydrologyMap(heightmap, riverMouthCandidates, new List<RiverEdge>(), new List<GraphEdge>());
        }

        List<DirectedNode> allRiverNodes = riverMouthCandidates.Select(x => x).ToList();// SelectRandomBorderCoordsAsDirectedNodes(borderCoords, 0.2, 0.5);
        Queue<DirectedNode> riverNodes = new Queue<DirectedNode>();
        riverMouthCandidates.ForEach(x => riverNodes.Enqueue(x));

        List<RiverEdge> riverEdges = new List<RiverEdge>();
//


        for (int i = 0; i < 1111; i++)
        {
            var node = riverNodes.Dequeue();
            List<DirectedNode> newRiverNodes =
                FindNewRiverNodes(heightmap, node, 20, 5, 0.2f, 3, allRiverNodes, 20);
//
            foreach (var newDirectedNode in newRiverNodes)
            {
                riverNodes.Enqueue(newDirectedNode);
                allRiverNodes.Add(newDirectedNode);
                riverEdges.Add(new RiverEdge(node, newDirectedNode));
            }

//
            if (riverNodes.Count == 0) break;
        }

        var voronoiCellCenters = allRiverNodes.Concat();

        var voronoiEdges = GenerateVoronoiEdges(allRiverNodes);
     //   voronoiEdges = voronoiEdges.Where(x => Distance((int)x.x1, (int)x.x2, (int)x.y1, (int)x.y2) < 51.0).ToList();
        return new HydrologyMap(heightmap, allRiverNodes, riverEdges, voronoiEdges);
        //return new HydrologyMap(heightmap, borderCoords, new List<RiverEdge>(), new List<GraphEdge>());
    }

    List<DirectedNode> GetRiverMouthCandidates(float[,] heightmap, int skipCount)
    {
        // Initialize variables
        int width = heightmap.GetLength(0);
        int height = heightmap.GetLength(1);
        List<DirectedNode> borderCoordinates = new List<DirectedNode>();

        // Find a border point by doing a quick scan through the middle of the heightmap
        Point firstBorderPoint = new Point(-1, -1);
        for (int x = 0; x < width; x++)
        {
            if (heightmap[x, height / 2] > 0)
            {
                firstBorderPoint = new Point(x, height / 2);
                break;
            }
        }

        if (firstBorderPoint.X == -1)
        {
            return borderCoordinates; // Return an empty list if no border point is found
        }

        // Crawl over boundary of shape to find border coordinates, adding every skipCount to a list as DirectedNode objects
        int skipCounter = 0;
        Point prevPoint = new Point(-1, -1);

        Point nextPoint = Point.Empty;
        Point currentPoint = firstBorderPoint;
        while (true)
        {
            Point[] neighborOffsets =
            {
                new(1, 0), new(0, 1), new(-1, 0), new(0, -1),
                new(1, 1), new(-1, 1), new(1, -1), new(-1, -1)
            };


            foreach (Point offset in neighborOffsets)
            {
                Point candidate = new Point(currentPoint.X + offset.X, currentPoint.Y + offset.Y);

                if (candidate != prevPoint && candidate != currentPoint && heightmap[candidate.X, candidate.Y] > 0)
                {
                    bool isBorderPoint = false;

                    // Check if the candidate point has an ocean neighbor
                    foreach (Point offset2 in neighborOffsets)
                    {
                        Point neighbor = new Point(candidate.X + offset2.X, candidate.Y + offset2.Y);

                        if (neighbor.X >= 0 && neighbor.X < width && neighbor.Y >= 0 && neighbor.Y < height)
                        {
                            if (heightmap[neighbor.X, neighbor.Y] == 0)
                            {
                                isBorderPoint = true;
                                break;
                            }
                        }
                    }

                    if (isBorderPoint)
                    {
                        nextPoint = candidate;
                        break;
                    }
                }
            }

            if (nextPoint == firstBorderPoint) break; // Stop when we are back where we started

            if (skipCounter == 0)
            {
                Vector2 oceanDirection = Vector2.Zero;


                // Calculate ocean direction except for the first node
                if (borderCoordinates.Count > 0)
                {
                    oceanDirection = CalculateOceanDirection(prevPoint, currentPoint, nextPoint, heightmap);
                }

                borderCoordinates.Add(new DirectedNode(currentPoint, NodeType.Border, oceanDirection));
            }


            skipCounter = (skipCounter + 1) % skipCount;

            // Move to the next border point

            prevPoint = currentPoint;
            currentPoint = nextPoint;
        }

        borderCoordinates.RemoveAt(borderCoordinates.Count - 1);

        // Handle the ocean direction for the first and last border coordinates
        if (borderCoordinates.Count > 1)
        {
            // Calculate ocean direction for the first coordinate
            Vector2 firstPointDirection = CalculateOceanDirection(borderCoordinates[borderCoordinates.Count - 1].Point,
                borderCoordinates[0].Point, borderCoordinates[1].Point, heightmap);
            borderCoordinates[0] = new DirectedNode(borderCoordinates[0].Point, NodeType.Border, firstPointDirection);
        }

        // Return result list
        return borderCoordinates;

        Vector2 CalculateOceanDirection(Point prevBorderPoint, Point currentBorderPoint, Point nextBorderPoint,
            float[,] heightmap)
        {
            float midX = (prevBorderPoint.X + nextBorderPoint.X) / 2f;
            float midY = (prevBorderPoint.Y + nextBorderPoint.Y) / 2f;

            float heightAtMidPoint = heightmap[(int)midX, (int)midY];

            Vector2 oceanDirection;

            if (heightAtMidPoint > 0)
            {
                oceanDirection = new Vector2(currentBorderPoint.X - midX, currentBorderPoint.Y - midY);
            }
            else
            {
                oceanDirection = new Vector2(midX - currentBorderPoint.X, midY - currentBorderPoint.Y);
            }

            return oceanDirection / oceanDirection.Length();
        }
    }
/*

static List<Point> GetBorderCoordinates(float[,] heightmap, int skipCount)
{
int width = heightmap.GetLength(0);
int height = heightmap.GetLength(1);
List<Point> borderCoordinates = new List<Point>();

int currentSkip = 0;

for (int y = 0; y < height; y++)
{
for (int x = 0; x < width; x++)
{
if (heightmap[x, y] != 0)
{
   bool hasZeroNeighbor = false;

   for (int dx = -1; dx <= 1; dx++)
   {
       for (int dy = -1; dy <= 1; dy++)
       {
           if (dx == 0 && dy == 0) continue;

           int nx = x + dx;
           int ny = y + dy;

           if (nx >= 0 && nx < width && ny >= 0 && ny < height && heightmap[nx, ny] == 0)
           {
               hasZeroNeighbor = true;
               break;
           }
       }

       if (hasZeroNeighbor && currentSkip == 0)
       {
           borderCoordinates.Add(new Point(x, y));
           break;
       }

       currentSkip++;
       if (currentSkip == skipCount) currentSkip = 0;
   }
}
}
}

return borderCoordinates;
}*/

// Calculate distance between two coordinates
    static double Distance(IGraphNode coord1, IGraphNode coord2)
    {
        int dx = coord1.Point.X - coord2.Point.X;
        int dy = coord1.Point.Y - coord2.Point.Y;
        return Math.Sqrt(dx * dx + dy * dy);
    }

    public static Vector2 GetVectorFromPoints(Point point1, Point point2)
    {
        float x = (float)point1.X - (float)point2.X;
        float y = (float)point1.Y - (float)point2.Y;
        return new Vector2(x, y);
    }

  //static double Distance(int x1, int x2, int y1, int y2)
  //{
  //    int dx = x1 - x2;
  //    int dy = y1 - y2;
  //    return Math.Sqrt(dx * dx + dy * dy);
  //}
    static double Distance(Point p1, Point p2)
    {
        int dx = p1.X - p2.X;
        int dy = p1.Y - p2.Y;
        return Math.Sqrt(dx * dx + dy * dy);
    }


    List<DirectedNode> SelectRandomBorderCoordsAsDirectedNodes(List<DirectedNode> borderCoords,
        double selectionPercentage,
        double prunePercentage)
    {
        int randomCount = (int)(borderCoords.Count * selectionPercentage);
        int prunedCount = (int)(randomCount * (1 - prunePercentage));

        // Select random border coordinates
        List<DirectedNode> randomBorderCoords = new List<DirectedNode>();
        while (randomBorderCoords.Count < randomCount)
        {
            int index = Random.Next(borderCoords.Count);
            var coord = borderCoords[index];
            if (!randomBorderCoords.Contains(coord))
            {
                randomBorderCoords.Add(coord);
            }
        }

        // Prune randomBorderCoords
        while (randomBorderCoords.Count > prunedCount)
        {
            double minDistance = double.MaxValue;
            int index1 = -1;
            int index2 = -1;

            // Find the two closest coordinates
            for (int i = 0; i < randomBorderCoords.Count; i++)
            {
                for (int j = i + 1; j < randomBorderCoords.Count; j++)
                {
                    double distance = Distance(randomBorderCoords[i], randomBorderCoords[j]);
                    if (distance < minDistance)
                    {
                        minDistance = distance;
                        index1 = i;
                        index2 = j;
                    }
                }
            }

            // Remove the two closest coordinates
            if (index1 >= 0 && index2 >= 0)
            {
                randomBorderCoords.RemoveAt(index2);
                randomBorderCoords.RemoveAt(index1);
            }
        }

        return randomBorderCoords;
    }

// ...

// In the Main() method or wherever you


    List<DirectedNode> FindNewRiverNodes(float[,] heightmap, DirectedNode drainNode, double distance,
        double variance,
        float minHeight, int maxTries, List<DirectedNode> riverNodes, float minDistToOtherNodes)
    {
        int width = heightmap.GetLength(0);
        int height = heightmap.GetLength(1);


        RiverAction action;
        var nextActionSelector = Random.NextDouble();
        if (drainNode.Type == NodeType.Border || drainNode.Type == NodeType.RiverMouth ||
            nextActionSelector < parameters.ContinueProbability)
        {
            action = RiverAction.Continue;
        }
        else if (nextActionSelector < parameters.ContinueProbability + parameters.AsymmetricBranchProbability)
        {
            action = RiverAction.AsymmetricBranch;
        }
        else
        {
            action = RiverAction.SymmetricBranch;
        }

        int x, y;
        double currentHeight = heightmap[drainNode.X, drainNode.Y];


        var newNodes = new List<DirectedNode>();

        for (int i = 0; i < maxTries; i++)
        {
            double randomDistance = distance + Random.NextDouble() * variance * 2 - variance;
            double angle = Random.NextDouble() * 2 * Math.PI;

            x = (int)Math.Round(drainNode.X + randomDistance * Math.Cos(angle));
            y = (int)Math.Round(drainNode.Y + randomDistance * Math.Sin(angle));
            var checkingPoint = new Point(x, y);
            if (x >= 0 && x < width && y >= 0 && y < height && heightmap[x, y] >= minHeight &&
                heightmap[x, y] > currentHeight)
            {
                bool isTooClose = false;
                foreach (DirectedNode riverNode in riverNodes)
                {
                    if (riverNode.Equals(drainNode)) continue;

                    if (Distance(checkingPoint, riverNode.Point) < minDistToOtherNodes)
                    {
                        isTooClose = true;
                        break;
                    }
                }

                if (!isTooClose)
                {
                    newNodes.Add(new DirectedNode(checkingPoint, NodeType.River,
                        GetVectorFromPoints(checkingPoint, drainNode.Point)));
                    if (newNodes.Count > 2 && action != RiverAction.Continue
                        || action == RiverAction.Continue)
                    {
                        return newNodes;
                    }
                }
            }
        }

        return newNodes;
    }

/*
    public static List<RiverEdge> GenerateVoronoiEdges(List<DirectedNode> riverNodes)
    {
        var voronoi = VoronoiMesh.Create<DirectedNode, DefaultTriangulationCell<DirectedNode>>(riverNodes);

        var edges = new List<RiverEdge>();
        foreach (var edge in voronoi.Edges)
        {
            var source = edge.Source.Vertices[0] as DirectedNode?;
            var target = edge.Target.Vertices[0] as DirectedNode?;
            edges.Add(new RiverEdge(source.Value, target.Value));
        }

        return edges;
    }*/

    public static List<GraphEdge> GenerateVoronoiEdges(List<DirectedNode> riverNodes)
    {
        var voroObject = new Voronoi(0.1);

        double[] xVal = riverNodes.Select(x => (double)x.X).ToArray();
        double[] yVal = riverNodes.Select(x => (double)x.Y).ToArray();

        List<GraphEdge> voronoi = voroObject.generateVoronoi(xVal, yVal, 0, 512, 0, 512);

        return voronoi;
    }

/*
    public static List<RiverEdge> GenerateVoronoiEdges(List<DirectedNode> riverNodes)
    {
        var geometryFactory = new GeometryFactory();
        var coordinates = riverNodes.Select(rn => new Coordinate(rn.X, rn.Y)).ToArray();
        var multiPoint = geometryFactory.CreateMultiPointFromCoords(coordinates);

        var voronoiBuilder = new VoronoiDiagramBuilder();
        voronoiBuilder.SetSites(multiPoint);

        var voronoiDiagram = voronoiBuilder.GetDiagram(geometryFactory);
        var voronoiEdges = new List<RiverEdge>();

        foreach (var geometry in voronoiDiagram.Geometries)
        {
            if (geometry is Polygon polygon)
            {
                var vertices = polygon.Coordinates;
                for (int i = 0; i < vertices.Length - 1; i++)
                {
                    var p1 = FindClosestDirectedNode(riverNodes, vertices[i]);
                    var p2 = FindClosestDirectedNode(riverNodes, vertices[i + 1]);
                    voronoiEdges.Add(new RiverEdge(p1, p2));
                }
            }
        }

        return voronoiEdges;
    }*/

/*
    public static List<RiverEdge> GenerateVoronoiEdges(List<DirectedNode> riverNodes)
    {
        var input = new Polygon(new LinearRing(new Coordinate[]{new Coordinate(0,0)}));

        foreach (var riverNode in riverNodes)
        {
            input.Coordinate.Add(new Vertex(riverNode.X, riverNode.Y));
        }

        var quality = new QualityOptions() { MinimumAngle = 25.0 };
        var mesh = (Mesh)input. .Triangulate(quality);

        var voronoi = mesh.GenerateVoronoi();

        var riverNodeDict = riverNodes.ToDictionary(rn => new Point((int)rn.X, (int)rn.Y));
        var edges = new List<RiverEdge>();

        foreach (var edge in voronoi.Edges)
        {
            if (edge.IsBoundary || edge.P0 == null || edge.P1 == null) continue;

            var source = riverNodeDict[new Point((int)edge.P0.X, (int)edge.P0.Y)];
            var target = riverNodeDict[new Point((int)edge.P1.X, (int)edge.P1.Y)];

            edges.Add(new RiverEdge(source, target));
        }

        return edges;
    }
    private static DirectedNode FindClosestDirectedNode(List<DirectedNode> riverNodes, Coordinate coordinate)
    {
        return riverNodes.OrderBy(rn => Math.Sqrt(Math.Pow(rn.X - coordinate.X, 2) + Math.Pow(rn.Y - coordinate.Y, 2))).First();
    }*/

    public void SaveMapAsPNG(HydrologyMap map, string outputPath)
    {
        var heightmap = map.HeightMap;
        int width = heightmap.GetLength(0);
        int height = heightmap.GetLength(1);


        using (Bitmap bitmap = new Bitmap(width, height))
        {
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    int grayValue = (int)(heightmap[x, y] * 255);
                    bitmap.SetPixel(x, y, Color.FromArgb(grayValue, grayValue, grayValue));
                }
            }


            foreach (var edge in map.RiverEdges)
            {
                DrawLineOnBitmap(bitmap, edge.P1.Point, edge.P2.Point, Color.Blue, 1);
            }

               foreach (var edge in map.VoronoiEdges)
               {
                   DrawLineOnBitmap(bitmap, new Point((int)edge.x1, (int)edge.y1),
                       new Point((int)edge.x2, (int)edge.y2),
                       Color.Yellow, 1);
               }


            foreach (var riverMouth in map.AllRiverNodes)
            {
                bitmap.SetPixel(riverMouth.X, riverMouth.Y, Color.Red);

                bitmap.SetPixel(riverMouth.X + 1, riverMouth.Y - 1, Color.Red);
                bitmap.SetPixel(riverMouth.X + 1, riverMouth.Y + 1, Color.Red);

                bitmap.SetPixel(riverMouth.X - 1, riverMouth.Y - 1, Color.Red);

                bitmap.SetPixel(riverMouth.X - 1, riverMouth.Y + 1, Color.Red);
            }


            bitmap.Save(outputPath, ImageFormat.Png);
        }
    }


    static void DrawLineOnBitmap(Bitmap bitmap, Point coord1, Point coord2, Color lineColor, int lineWidth)
    {
        using (Graphics graphics = Graphics.FromImage(bitmap))
        {
            using (Pen pen = new Pen(lineColor, lineWidth))
            {
                graphics.DrawLine(pen, coord1.X, coord1.Y, coord2.X, coord2.Y);
            }
        }
    }
}