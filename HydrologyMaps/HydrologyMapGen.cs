using System.Drawing;
using System.Drawing.Imaging;
using System.Numerics;
using Voronoi2;

namespace HydrologyMaps;

public class HydrologyMapGen
{
    private readonly HydrologyParameters parameters;
    Random Random = new Random();

    public HydrologyMapGen(HydrologyParameters parameters)
    {
        this.parameters = parameters;
    }


    public HydrologyMap GenerateIsland(int width, int height, int seed)
    {
        if (seed == -1) seed = Random.Next();
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

        (List<DirectedNode> riverMouthCandidates, List<GraphNode> extraVoronoiPoints) =
            GetRiverMouthCandidates(heightmap, parameters.SpaceBetweenRiverMouthCandidates);

        if (riverMouthCandidates.Count == 0)
        {
            return new HydrologyMap(heightmap, riverMouthCandidates, new List<RiverEdge>(), new List<GraphEdge>());
        }

        List<DirectedNode>
            allRiverNodes =
                riverMouthCandidates.Select(x => x)
                    .ToList(); // SelectRandomBorderCoordsAsDirectedNodes(borderCoords, 0.2, 0.5);
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
                node.Type = NodeType.RiverMouth;
                riverNodes.Enqueue(newDirectedNode);
                allRiverNodes.Add(newDirectedNode);
                riverEdges.Add(new RiverEdge(node, newDirectedNode));
            }

//
            if (riverNodes.Count == 0) break;
        }


        List<IGraphNode> allPointsForVoronoi = new List<IGraphNode>();
        allRiverNodes.ForEach(x => allPointsForVoronoi.Add(x));
        extraVoronoiPoints.ForEach(x => allPointsForVoronoi.Add(x));
        List<GraphEdge> voronoiEdges = IslandVoronoi.GenerateVoronoiEdges(allPointsForVoronoi);
        voronoiEdges = IslandVoronoi.GenerateExtraEdges(riverMouthCandidates, voronoiEdges, heightmap);
        //   voronoiEdges = voronoiEdges.Where(x => Distance((int)x.x1, (int)x.x2, (int)x.y1, (int)x.y2) < 51.0).ToList();
        
        allRiverNodes.RemoveAll(x => x.Type == NodeType.Border);
        
        return new HydrologyMap(heightmap, allRiverNodes, riverEdges, voronoiEdges);
        //return new HydrologyMap(heightmap, borderCoords, new List<RiverEdge>(), new List<GraphEdge>());
    }

    (List<DirectedNode> borderCoordinates, List<GraphNode> extraVoronoiPoints) GetRiverMouthCandidates(
        float[,] heightmap, int skipCount, int maxConcaveness = 3)
    {
        // Initialize variables
        int width = heightmap.GetLength(0);
        int height = heightmap.GetLength(1);
        List<DirectedNode> borderCoordinates = new List<DirectedNode>();
        List<GraphNode> extraVoronoiPoints = new List<GraphNode>();

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
            return (borderCoordinates, extraVoronoiPoints); // Return empty lists if no border point is found
        }

        // Crawl over boundary of shape to find border coordinates, adding every skipCount to a list as DirectedNode objects
        int skipCounter = 0;

        Point currentPoint = firstBorderPoint;
        Vector2 oceanDirection = new Vector2(-1, 0);

        var facing = new Point(0, -1); // up

        int k = 0;

        int concaveness = 0; // negative value means we are turning left more and in a concave

        Point lastAddedNode = Point.Zero;

        // used for detecting concaves and adding additional nodes
        int trailingPointCount = skipCount;
        List<Point> trailingPoints = new List<Point>(trailingPointCount);


        while (true)
        {
            /*
                ↑ ▮   ▮▮         ▮
                | ▮    ↰▮   ⮣▮   ⮣▮ 
                1       2     3    4
                If there's a wall ahead on the right and no wall blocking the path, move forward one cell (1).
                If there's a wall ahead on the right and a wall blocking the path, add a point and turn left (2).
                Otherwise add a point and turn right (3) (also covers 4)
                Note that Y = 0 is the top
                Note that we do point math as if we are crawling in the wall, not on the space next to it is changed a bit
             */

            Point left = new Point(facing.Y, -facing.X);
            Point right = -left;

            // if there is more land ahead
            if (Blocked(currentPoint + facing))
            {
                if (Blocked(currentPoint + facing + left))
                {
                    // If there was also land on the left, then we turn left so stay on the coast
                    currentPoint += left; // turn left
                    Point leftBack = left - facing;
                    facing = left;
                    oceanDirection = new Vector2(leftBack.X, leftBack.Y);
                    concaveness--;
                }
                else
                {
                    // No land on the left so we go forward and will stay on the coast
                    currentPoint += facing; // forward
                    oceanDirection = new Vector2(left.X, left.Y);
                }
            }
            else
            {
                // no land ahead, got to turn right
                currentPoint += right;
                facing = right;
                oceanDirection = new Vector2(left.X, left.Y);
                concaveness++;
            }

            trailingPoints.Add(currentPoint);

            if (currentPoint == firstBorderPoint) break; // Stop when we are back where we started


            if (skipCounter == 0)
            {
                // First we need to check if we want to add any additional additional nodes before current to cover concave shape
                if (trailingPoints.Count >= trailingPointCount && borderCoordinates.Count > 0)
                {
                    // recursively add midway points while the midway point is in a significantly convex area
                    int startIndex = 0;
                    int endIndex = trailingPoints.Count - 1;
                    int midTrailingPoint = endIndex / 2;
                    AddAdditionalNodesIfNeeded(borderCoordinates, extraVoronoiPoints, trailingPoints, startIndex,
                        endIndex,
                        midTrailingPoint);

                    trailingPoints.Clear();
                }

                //  var riverNodePoint = new Point(currentPoint.X - (int)oceanDirection.X*2, currentPoint.Y- (int)oceanDirection.Y*2);
                concaveness = 0;
                oceanDirection /= oceanDirection.Length(); // normalize
                borderCoordinates.Add(new DirectedNode(currentPoint, NodeType.Border, oceanDirection));
                lastAddedNode = currentPoint;
            }

            skipCounter = (skipCounter + 1) % skipCount;

            bool Blocked(Point p) => heightmap[p.X, p.Y] > 0;
        }


        // recursively add midway points while the midway point is in a significantly convex area
        int startIndexF = 0;
        int endIndexF = trailingPoints.Count - 1;
        int midIndexF = endIndexF / 2;
        AddAdditionalNodesIfNeeded(borderCoordinates, extraVoronoiPoints, trailingPoints, startIndexF, endIndexF,
            midIndexF);

        // Return result list
        return (borderCoordinates, extraVoronoiPoints);
    }

    private void AddAdditionalNodesIfNeeded(List<DirectedNode> concaveExtraNodes, List<GraphNode> convexExtraNodes,
        List<Point> trailingPoints,
        int startIndex, int endIndex, int midIndex)
    {
        var midPoint = trailingPoints[midIndex];
        var startPoint = trailingPoints[startIndex];
        var endPoint = trailingPoints[endIndex];
        (double dist, bool isLeft) = Vector2D.DistanceBetweenPointAndLineSegment(startPoint, endPoint, midPoint);

        if (dist > 1)
        {
            int newMid1 = startIndex + (midIndex - startIndex) / 2;
            int newMid2 = midIndex + (endIndex - midIndex) / 2;
            AddAdditionalNodesIfNeeded(concaveExtraNodes, convexExtraNodes, trailingPoints, startIndex, midIndex,
                newMid1);
            AddAdditionalNodesIfNeeded(concaveExtraNodes, convexExtraNodes, trailingPoints, midIndex, endIndex,
                newMid2);

            if (isLeft)
                concaveExtraNodes.Add(new DirectedNode(midPoint, NodeType.Border, Vector2.One));
            else
                convexExtraNodes.Add(new GraphNode(midPoint, NodeType.Border));
        }
    }

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
}