﻿using System.Drawing;
using System.Drawing.Imaging;
using System.Numerics;
using NetTopologySuite.Geometries;
using NetTopologySuite.Index.KdTree;
using Voronoi2;

namespace HydrologyMaps;

public class HydrologyMapGen
{
    private readonly HydrologyParameters parameters;
    Random Random = new Random();

    public static List<IEdge> DebugGraphEdges1 = new();
    public static List<IEdge> DebugGraphEdges2 = new();

    public HydrologyMapGen(HydrologyParameters parameters)
    {
        this.parameters = parameters;
    }


    public HydrologyMap GenerateIsland(int width, int seed, int k)
    {
        if (seed == -1)
        {
            seed = Random.Next();
            Console.WriteLine("Seed: " + seed);
        }

        Random = new Random(seed);


        float[,] heightmap = new float[width, width];

        HydrologyTerrainFormer.FillBaseHeightmap(heightmap, width, Random, parameters);

        
      // return new HydrologyMap(heightmap, new List<DirectedNode>(), new List<RiverEdge>(), new GridCell[0, 0],
      //     new List<GraphEdge>(),
      //     new List<IGraphNode>());
        
        (List<DirectedNode> riverMouthCandidates, List<GraphNode> extraVoronoiPoints) =
            GetRiverMouthCandidates(heightmap, parameters.SpaceBetweenRiverMouthCandidates);

        if (riverMouthCandidates.Count == 0)
        {
            return new HydrologyMap(heightmap, riverMouthCandidates, new List<RiverEdge>(), new GridCell[0, 0],
                new List<GraphEdge>(),
                new List<IGraphNode>());
        }

        
        HydrologyKdTree2D borderKdTree = new HydrologyKdTree2D();
        riverMouthCandidates.ForEach(n => borderKdTree.Insert(n));
        extraVoronoiPoints.ForEach(n => borderKdTree.Insert(n));
        
        HydrologyKdTree2D hydrologyKdTree = new HydrologyKdTree2D();
        (List<DirectedNode> allRiverNodes, List<RiverEdge> riverEdges) =
            ExpandRiverGrid(heightmap, riverMouthCandidates, k);


        allRiverNodes.ForEach(n => hydrologyKdTree.Insert(n));

        List<IGraphNode> allPointsForVoronoi = new List<IGraphNode>();
        allRiverNodes.ForEach(x => allPointsForVoronoi.Add(x));
        extraVoronoiPoints.ForEach(x => allPointsForVoronoi.Add(x));

        List<GraphEdge> voronoiEdges = IslandVoronoi.GenerateVoronoiEdges(allPointsForVoronoi);
        voronoiEdges = IslandVoronoi.GenerateExtraEdges(riverMouthCandidates, voronoiEdges, heightmap);
        //   voronoiEdges = voronoiEdges.Where(x => Distance((int)x.x1, (int)x.x2, (int)x.y1, (int)x.y2) < 51.0).ToList();

        allRiverNodes.RemoveAll(x => x.Type == NodeType.Border);

        GridCell[,] gridCells = new GridCell[width, width];


        allRiverNodes.ForEach(x => hydrologyKdTree.Insert(x));

        FindAreaForNodes(heightmap, hydrologyKdTree, allRiverNodes, width, width);
        // flowrate = 0.42 · A^0.69 formula from hydrology paper. Before this point FlowRate holds the number of pixels for which the closest node is this
        allRiverNodes.ForEach(n => n.FlowRate = 0.42 * Math.Pow(n.FlowRate, 0.69));
        riverMouthCandidates.ForEach(n => SetCumulativeFlowRate(n));

        GenerateRiverNodeHeights(riverMouthCandidates);

        HydrologyTerrainFormer.InterpolateHeightMapSimple(heightmap, borderKdTree, allRiverNodes, width, 10, 5f);

        HydrologyTerrainFormer.CarveRivers(heightmap, riverMouthCandidates/*.Skip(1).Take(1).ToList()*/, 0.5f, 0.35f, 0.3f, 1f, 8f);

        return new HydrologyMap(heightmap, allRiverNodes, riverEdges, gridCells, voronoiEdges, allPointsForVoronoi);
    }

    private double SetCumulativeFlowRate(DirectedNode node)
    {
        if (node.Children.Count > 0)
            node.FlowRate = node.Children.Sum(n => SetCumulativeFlowRate((DirectedNode)n));
        return node.FlowRate;
    }

    private void GenerateRiverNodeHeights(List<DirectedNode> riverMouthCandidates)
    {
        riverMouthCandidates.ForEach(n => GenerateRiverNodeHeightsRecursive(n, 0));

        void GenerateRiverNodeHeightsRecursive(DirectedNode node, float height)
        {
            node.Height = height;
            if (node.Children.Count > 0)
            {
                float newHeight = node.FlowRate > parameters.ProperRiverMinimumFlowRate
                    ? 0
                    : Math.Min(1, height + 0.1f * (float)Random.NextDouble());
                node.Children.ForEach(
                    n => GenerateRiverNodeHeightsRecursive((DirectedNode)n, newHeight));
            }
        }
    }


    (List<DirectedNode> borderCoordinates, List<GraphNode> extraVoronoiPoints) GetRiverMouthCandidates(
        float[,] heightmap, int skipCount)
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

            currentPoint.X = Math.Clamp(currentPoint.X, 0, width);
            currentPoint.Y = Math.Clamp(currentPoint.Y, 0, width);

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
                        midTrailingPoint, parameters.MaxNodePriority);
                }

                //  var riverNodePoint = new Point(currentPoint.X - (int)oceanDirection.X*2, currentPoint.Y- (int)oceanDirection.Y*2);
                concaveness = 0;
                oceanDirection /= oceanDirection.Length(); // normalize
                borderCoordinates.Add(new DirectedNode(currentPoint, NodeType.Border, null, oceanDirection,
                    Random.Next(parameters.MaxNodePriority, parameters.MaxNodePriority)));
                lastAddedNode = currentPoint;
                trailingPoints.Clear();
            }

            skipCounter = (skipCounter + 1) % skipCount;

            bool Blocked(Point p) => (p.X >= width || p.X < 0 || p.Y >= width || p.Y < 0) || heightmap[p.X, p.Y] > 0;
        }


        // recursively add midway points while the midway point is in a significantly convex area
        int startIndexF = 0;
        int endIndexF = trailingPoints.Count - 1;
        int midIndexF = endIndexF / 2;
        AddAdditionalNodesIfNeeded(borderCoordinates, extraVoronoiPoints, trailingPoints, startIndexF, endIndexF,
            midIndexF, Random.Next(2, parameters.MaxNodePriority));

        // Return result list
        return (borderCoordinates, extraVoronoiPoints);
    }


    private void FindAreaForNodes(float[,] heightmap, HydrologyKdTree2D kdTree, List<DirectedNode> allNodes, int width,
        int height)
    {
        // Iterate through the grid points
        for (int y = 0; y < height; y++)
        {
            for (int x = 0; x < width; x++)
            {
                if (heightmap[x, y] == 0) continue;

                // Add this pixel to the area of the closest rivernode (FlowRate is used for area at this point)
                List<GraphNode> nearestNeighbors = kdTree.FindClosest(x, y, 1);
                ((DirectedNode)nearestNeighbors.First()).FlowRate++;
            }
        }
    }

    private void AddAdditionalNodesIfNeeded(List<DirectedNode> concaveExtraNodes, List<GraphNode> convexExtraNodes,
        List<Point> trailingPoints,
        int startIndex, int endIndex, int midIndex, int newNodePriority)
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
                newMid1, newNodePriority);

            if (isLeft)
                concaveExtraNodes.Add(new DirectedNode(midPoint, NodeType.Border, null, Vector2.One, newNodePriority));
            else
                convexExtraNodes.Add(new GraphNode(midPoint, NodeType.Border, null));

            AddAdditionalNodesIfNeeded(concaveExtraNodes, convexExtraNodes, trailingPoints, midIndex, endIndex,
                newMid2, newNodePriority);
        }
    }

    private (List<DirectedNode>, List<RiverEdge>) ExpandRiverGrid(
        float[,] heightmap,
        List<DirectedNode> riverMouthCandidates, int k)
    {
        List<DirectedNode> riverGrid = new List<DirectedNode>();
        List<RiverEdge> riverEdges = new List<RiverEdge>();

        PriorityQueue<DirectedNode, int> riverExpandQueue = new();
        //Queue<DirectedNode> riverNodes = new Queue<DirectedNode>();
        riverMouthCandidates.ForEach(n =>
        {
            riverExpandQueue.Enqueue(n, parameters.MaxNodePriority - n.Priority);
            riverGrid.Add(n);
        });

        for (int i = 0; i < k; i++)
        {
            var node = riverExpandQueue.Dequeue();
            List<DirectedNode> newRiverNodes =
                ContinueRiver(heightmap, node, riverGrid, riverEdges);

            foreach (var newDirectedNode in newRiverNodes)
            {
                node.Type = NodeType.River;
                riverExpandQueue.Enqueue(newDirectedNode, parameters.MaxNodePriority - newDirectedNode.Priority);
                riverGrid.Add(newDirectedNode);
                riverEdges.Add(new RiverEdge(node, newDirectedNode, newDirectedNode.Priority));
            }


            if (riverExpandQueue.Count == 0) break;
        }

        return (riverGrid, riverEdges);
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


    List<DirectedNode> ContinueRiver(float[,] heightmap, DirectedNode drainNode,
        List<DirectedNode> riverNodes, List<RiverEdge> riverEdgesSoFar)
    {
        var newNodes = new List<DirectedNode>();

        int width = heightmap.GetLength(0);
        int height = heightmap.GetLength(1);

        Vector2 directionVector = drainNode.Direction;
        directionVector /= directionVector.Length();

        int branch1Priority = drainNode.Priority;
        int branch2Priority = 0;

        RiverAction action;
        var nextActionSelector = Random.NextDouble();
        if (drainNode.Type == NodeType.Border || drainNode.Type == NodeType.RiverMouth ||
            nextActionSelector < parameters.ContinueProbability ||
            drainNode.Priority == 1)
        {
            action = RiverAction.Continue;
        }
        else if (nextActionSelector < parameters.ContinueProbability + parameters.AsymmetricBranchProbability)
        {
            action = RiverAction.AsymmetricBranch;
            branch2Priority = Random.Next(1, drainNode.Priority - 1);
        }
        else
        {
            action = RiverAction.SymmetricBranch;
            branch1Priority = drainNode.Priority - 1;
            branch2Priority = drainNode.Priority - 1;
        }


        double currentHeight = heightmap[drainNode.X, drainNode.Y];


        //    double randomDistance1 = distance + Random.NextDouble() * variance * 2 - variance;
        //    double angle1 = Random.NextDouble() * 2 * Math.PI;
        //    
        //    if (action == RiverAction.Continue)
        //    {
        //        
        //    }
        //    else
        //    {
        //        double randomDistance2 = distance + Random.NextDouble() * variance * 2 - variance;
        //        double angle1 = Random.NextDouble() * 2 * Math.PI;
        //    }

        double randomDistance1 = parameters.InterNodeDist + Random.NextDouble() * parameters.InterNodeDistVar * 2 -
                                 parameters.InterNodeDistVar;


        double angleBetweenBranches = RandomInRange(Random, parameters.MinAngleBetweenBranches, 170);

        double angle1 = action == RiverAction.Continue
            ? RandomInRange(Random, -90, 90)
            : RandomInRange(Random, -90, Math.Min(0, 90 - angleBetweenBranches));
        double angle1Rad = angle1 * (Math.PI / 180);

        Point newPoint = new(
            drainNode.X + (int)(directionVector.X * randomDistance1 * Math.Cos(angle1Rad) -
                                directionVector.Y * randomDistance1 * Math.Sin(angle1Rad)),
            drainNode.Y + (int)(directionVector.X * randomDistance1 * Math.Sin(angle1Rad) +
                                directionVector.Y * randomDistance1 * Math.Cos(angle1Rad)));


        for (int i = 0; i < parameters.NodeExpansionMaxTries; i++)
        {
            if (ValidateNewPoint(newPoint, riverEdgesSoFar))
            {
                newNodes.Add(new DirectedNode(newPoint, NodeType.River, drainNode,
                    GetVectorFromPoints(newPoint, drainNode.Point), Math.Max(1, branch1Priority)));
                break;
            }

            randomDistance1 = parameters.InterNodeDist + Random.NextDouble() * parameters.InterNodeDistVar * 2 -
                              parameters.InterNodeDistVar;
            angle1 = drainNode.Type == NodeType.Border ? Random.NextDouble() * 360 : Random.NextDouble() * -90;
            angle1Rad = angle1 * (Math.PI / 180);
            directionVector = drainNode.Direction;
            directionVector /= directionVector.Length();
            newPoint = new(
                drainNode.X + (int)(directionVector.X * randomDistance1 * Math.Cos(angle1Rad) -
                                    directionVector.Y * randomDistance1 * Math.Sin(angle1Rad)),
                drainNode.Y + (int)(directionVector.X * randomDistance1 * Math.Sin(angle1Rad) +
                                    directionVector.Y * randomDistance1 * Math.Cos(angle1Rad)));
        }


        if (action != RiverAction.Continue)
        {
            double randomDistance2 = parameters.InterNodeDist + Random.NextDouble() * parameters.InterNodeDistVar * 2 -
                                     parameters.InterNodeDistVar;


            // var minAngle2 = angle1 + parameters.MinAngleBetweenBranches;


            for (int i = 0; i < parameters.NodeExpansionMaxTries * 3; i++)
            {
                double angle2Rad = (angle1 + angleBetweenBranches) * (Math.PI / 180);

                newPoint = new(
                    drainNode.X + (int)(directionVector.X * randomDistance2 * Math.Cos(angle2Rad) -
                                        directionVector.Y * randomDistance2 * Math.Sin(angle2Rad)),
                    drainNode.Y + (int)(directionVector.X * randomDistance2 * Math.Sin(angle2Rad) +
                                        directionVector.Y * randomDistance2 * Math.Cos(angle2Rad)));

                if (ValidateNewPoint(newPoint, riverEdgesSoFar))
                {
                    newNodes.Add(new DirectedNode(newPoint, NodeType.River, drainNode,
                        GetVectorFromPoints(newPoint, drainNode.Point), branch2Priority));
                    break;
                }

                randomDistance2 = parameters.InterNodeDist + Random.NextDouble() * parameters.InterNodeDistVar * 2 -
                                  parameters.InterNodeDistVar;

                angleBetweenBranches = RandomInRange(Random, parameters.MinAngleBetweenBranches, 90);
            }
        }

        return newNodes;

        double RandomInRange(Random random, double min, double max)
        {
            var result = random.NextDouble() * (max - min) + min;
            return result;
        }


        bool ValidateNewPoint(Point p, List<RiverEdge> edgesSoFar)
        {
            bool valid = false;
            var x = p.X;
            var y = p.Y;
            if (x >= 0 && x < width && y >= 0 && y < width && heightmap[x, y] >= parameters.NodeExpansionMinHeight &&
                y < height /* && 
                heightmap[x, y] >= currentHeight*/)
            {
                bool isTooClose = false;
                foreach (DirectedNode riverNode in riverNodes)
                {
                    if (riverNode.Equals(drainNode)) continue;

                    if (Point.Distance(p, riverNode.Point) < parameters.MinDistanceBetweenRiverNodes)
                    {
                        isTooClose = true;
                        break;
                    }
                }

                if (!isTooClose)
                {
                    valid = true;
                }
            }

            // finally check we dont cross other river edges
            if (valid)
            {
                var candidateEdge = new GraphEdge(p.X, p.Y, drainNode.X, drainNode.Y);
                RiverEdge? collide = edgesSoFar.FirstOrDefault(e =>
                    e.P2 != drainNode && new GraphEdge(e.P1.X, e.P1.Y, e.P2.X, e.P2.Y).Intersection(candidateEdge) !=
                    null);

                valid = edgesSoFar.All(e =>
                    e.P2 == drainNode ||
                    new GraphEdge(e.P1.X, e.P1.Y, e.P2.X, e.P2.Y).Intersection(new GraphEdge(p.X, p.Y, drainNode.X,
                        drainNode.Y)) == null);
            }

            return valid;
        }
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