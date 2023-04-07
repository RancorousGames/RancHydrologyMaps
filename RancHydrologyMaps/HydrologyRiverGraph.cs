using System.Numerics;
using HydrologyMaps;
using Voronoi2;
using Point = HydrologyMaps.Point;

namespace RancHydrologyMaps;

public class HydrologyRiverGraph
{
    private readonly HydrologyParameters _parameters;
    private readonly Random _random;

    public HydrologyRiverGraph(HydrologyParameters parameters, Random random)
    {
        _parameters = parameters;
        _random = random;
    }

    public double SetCumulativeFlowRate(DirectedNode node)
    {
        if (node.Children.Count > 0)
            node.FlowRate = node.Children.Sum(n => SetCumulativeFlowRate((DirectedNode)n));
        return node.FlowRate;
    }

    public void GenerateRiverNodeHeights(List<DirectedNode> riverMouthCandidates)
    {
        riverMouthCandidates.ForEach(n => GenerateRiverNodeHeightsRecursive(n, 0));

        void GenerateRiverNodeHeightsRecursive(DirectedNode node, float height)
        {
            node.Height = height;
            if (node.Children.Count > 0)
            {
                float newHeight = node.FlowRate > _parameters.ProperRiverMinimumFlowRate
                    ? 0
                    : Math.Min(1, height + 0.1f * (float)_random.NextDouble());
                node.Children.ForEach(
                    n => GenerateRiverNodeHeightsRecursive((DirectedNode)n, newHeight));
            }
        }
    }


    public (List<DirectedNode> borderCoordinates, List<GraphNode> extraVoronoiPoints)
        GetRiverMouthCandidates(float[,] heightmap, int skipCount, KdTree<Point> allBorderPointKdTree)
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

        Point currentPoint = firstBorderPoint;
        Vector2 oceanDirection = new Vector2(-1, 0);

        var facing = new Point(0, -1); // up

        int concaveness = 0; // negative value means we are turning left more and in a concave

        Point lastAddedNode = Point.Zero;

        // saving the last skipCount borderPoints used for detecting concaves and adding additional nodes
        int trailingPointCount = skipCount;
        List<Point> trailingPoints = new List<Point>(trailingPointCount);

        int i = 0;
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

            if (i % 8 == 0)
                allBorderPointKdTree.Insert(currentPoint);

            if (currentPoint == firstBorderPoint) break; // Stop when we are back where we started


            if (i % skipCount == 0)
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
                        midTrailingPoint, _parameters.MaxNodePriority);
                }

                //  var riverNodePoint = new Point(currentPoint.X - (int)oceanDirection.X*2, currentPoint.Y- (int)oceanDirection.Y*2);
                concaveness = 0;
                oceanDirection /= oceanDirection.Length(); // normalize
                borderCoordinates.Add(new DirectedNode(currentPoint, NodeType.Border, null, oceanDirection,
                    _random.Next(_parameters.MaxNodePriority, _parameters.MaxNodePriority)));
                lastAddedNode = currentPoint;
                trailingPoints.Clear();
            }

            ++i;

            bool Blocked(Point p) => (p.X >= width || p.X < 0 || p.Y >= width || p.Y < 0) || heightmap[p.X, p.Y] > 0;
        }


        // recursively add midway points while the midway point is in a significantly convex area
        int startIndexF = 0;
        int endIndexF = trailingPoints.Count - 1;
        int midIndexF = endIndexF / 2;
        AddAdditionalNodesIfNeeded(borderCoordinates, extraVoronoiPoints, trailingPoints, startIndexF, endIndexF,
            midIndexF, _random.Next(2, _parameters.MaxNodePriority));

        // Return result list
        return (borderCoordinates, extraVoronoiPoints);
    }


    public void FindAreaForNodes(float[,] heightmap, KdTree<GraphNode> graphNodeKdTree,
        List<DirectedNode> allNodes, int width,
        int height)
    {
        // Iterate through the grid points
        for (int y = 0; y < height; y++)
        {
            for (int x = 0; x < width; x++)
            {
                if (heightmap[x, y] == 0) continue;

                // Add this pixel to the area of the closest rivernode (FlowRate is used for area at this point)
                List<GraphNode> nearestNeighbors = graphNodeKdTree.FindClosest(x, y, 1);
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

    public (List<DirectedNode>, List<RiverEdge>) ExpandRiverGrid(
        float[,] heightmap,
        List<DirectedNode> riverMouthCandidates, int k)
    {
        List<DirectedNode> riverGrid = new List<DirectedNode>();
        List<RiverEdge> riverEdges = new List<RiverEdge>();

        PriorityQueue<DirectedNode, int> riverExpandQueue = new();
        //Queue<DirectedNode> riverNodes = new Queue<DirectedNode>();
        riverMouthCandidates.ForEach(n =>
        {
            riverExpandQueue.Enqueue(n, _parameters.MaxNodePriority - n.Priority);
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
                riverExpandQueue.Enqueue(newDirectedNode, _parameters.MaxNodePriority - newDirectedNode.Priority);
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


    List<DirectedNode> ContinueRiver(float[,] heightmap, DirectedNode drainNode,
        List<DirectedNode> riverNodes, List<RiverEdge> riverEdgesSoFar)
    {
        var newNodes = new List<DirectedNode>();

        int width = heightmap.GetLength(0);

        Vector2 directionVector = drainNode.Direction;
        directionVector /= directionVector.Length();

        int branch1Priority = drainNode.Priority;
        int branch2Priority = 0;

        RiverAction action;
        var nextActionSelector = _random.NextDouble();
        if (drainNode.Type == NodeType.Border || drainNode.Type == NodeType.RiverMouth ||
            nextActionSelector < _parameters.ContinueProbability ||
            drainNode.Priority == 1)
        {
            action = RiverAction.Continue;
        }
        else if (nextActionSelector < _parameters.ContinueProbability + _parameters.AsymmetricBranchProbability)
        {
            action = RiverAction.AsymmetricBranch;
            branch2Priority = _random.Next(1, drainNode.Priority - 1);
        }
        else
        {
            action = RiverAction.SymmetricBranch;
            branch1Priority = drainNode.Priority - 1;
            branch2Priority = drainNode.Priority - 1;
        }


        double currentHeight = heightmap[drainNode.X, drainNode.Y];


        double _randomDistance1 = _parameters.InterNodeDist + _random.NextDouble() * _parameters.InterNodeDistVar * 2 -
                                  _parameters.InterNodeDistVar;


        double angleBetweenBranches = RandomInRange( _parameters.MinAngleBetweenBranches, 170);

        double angle1 = action == RiverAction.Continue
            ? RandomInRange(-90, 90)
            : RandomInRange(-90, Math.Min(0, 90 - angleBetweenBranches));
        double angle1Rad = angle1 * (Math.PI / 180);

        Point newPoint = new(
            drainNode.X + (int)(directionVector.X * _randomDistance1 * Math.Cos(angle1Rad) -
                                directionVector.Y * _randomDistance1 * Math.Sin(angle1Rad)),
            drainNode.Y + (int)(directionVector.X * _randomDistance1 * Math.Sin(angle1Rad) +
                                directionVector.Y * _randomDistance1 * Math.Cos(angle1Rad)));


        for (int i = 0; i < _parameters.NodeExpansionMaxTries; i++)
        {
            if (ValidateNewPoint(newPoint, riverEdgesSoFar))
            {
                newNodes.Add(new DirectedNode(newPoint, NodeType.River, drainNode,
                    GetVectorFromPoints(newPoint, drainNode.Point), Math.Max(1, branch1Priority)));
                break;
            }

            _randomDistance1 = _parameters.InterNodeDist + _random.NextDouble() * _parameters.InterNodeDistVar * 2 -
                               _parameters.InterNodeDistVar;
            angle1 = drainNode.Type == NodeType.Border ? _random.NextDouble() * 360 : _random.NextDouble() * -90;
            angle1Rad = angle1 * (Math.PI / 180);
            directionVector = drainNode.Direction;
            directionVector /= directionVector.Length();
            newPoint = new(
                drainNode.X + (int)(directionVector.X * _randomDistance1 * Math.Cos(angle1Rad) -
                                    directionVector.Y * _randomDistance1 * Math.Sin(angle1Rad)),
                drainNode.Y + (int)(directionVector.X * _randomDistance1 * Math.Sin(angle1Rad) +
                                    directionVector.Y * _randomDistance1 * Math.Cos(angle1Rad)));
        }


        if (action != RiverAction.Continue)
        {
            double _randomDistance2 = _parameters.InterNodeDist +
                                      _random.NextDouble() * _parameters.InterNodeDistVar * 2 -
                                      _parameters.InterNodeDistVar;


            // var minAngle2 = angle1 + _parameters.MinAngleBetweenBranches;


            for (int i = 0; i < _parameters.NodeExpansionMaxTries * 3; i++)
            {
                double angle2Rad = (angle1 + angleBetweenBranches) * (Math.PI / 180);

                newPoint = new(
                    drainNode.X + (int)(directionVector.X * _randomDistance2 * Math.Cos(angle2Rad) -
                                        directionVector.Y * _randomDistance2 * Math.Sin(angle2Rad)),
                    drainNode.Y + (int)(directionVector.X * _randomDistance2 * Math.Sin(angle2Rad) +
                                        directionVector.Y * _randomDistance2 * Math.Cos(angle2Rad)));

                if (ValidateNewPoint(newPoint, riverEdgesSoFar))
                {
                    newNodes.Add(new DirectedNode(newPoint, NodeType.River, drainNode,
                        GetVectorFromPoints(newPoint, drainNode.Point), branch2Priority));
                    break;
                }

                _randomDistance2 = _parameters.InterNodeDist + _random.NextDouble() * _parameters.InterNodeDistVar * 2 -
                                   _parameters.InterNodeDistVar;

                angleBetweenBranches = RandomInRange(_parameters.MinAngleBetweenBranches, 90);
            }
        }

        return newNodes;

        double RandomInRange(double min, double max)
        {
            var result = _random.NextDouble() * (max - min) + min;
            return result;
        }


        bool ValidateNewPoint(Point p, List<RiverEdge> edgesSoFar)
        {
            bool valid = false;
            var x = p.X;
            var y = p.Y;

            bool validCoordinate = x >= 0 && x < width && y >= 0 && y < width;
            float pHeight = validCoordinate ? heightmap[x, y] : 0;
            if (validCoordinate && 
                pHeight  >= _parameters.NodeExpansionMinHeight &&
                pHeight > currentHeight - 0.06f && 
                (pHeight > _parameters.NodeExpansionMinHeight * 3 || pHeight > currentHeight))
            {
                bool isTooClose = false;
                foreach (DirectedNode riverNode in riverNodes)
                {
                    if (riverNode.Equals(drainNode)) continue;

                    if (Point.Distance(p, riverNode.Point) < _parameters.MinDistanceBetweenRiverNodes)
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
}