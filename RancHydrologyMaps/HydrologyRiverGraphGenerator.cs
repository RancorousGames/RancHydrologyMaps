using System.Numerics;
using HydrologyMaps;
using Voronoi2;
using Point = HydrologyMaps.Point;

namespace RancHydrologyMaps;

// responsible for generating the graph of rivers given a heightmap without rivers or lakes
public class HydrologyRiverGraphGenerator
{
    private readonly HydrologyParameters _parameters;
    private readonly Random _random;

    public HydrologyRiverGraphGenerator(HydrologyParameters parameters, Random random)
    {
        _parameters = parameters;
        _random = random;
    }

    // Find potential points along the border of a shape that could become river mouths. We simply take every n'th border point and call it a mouth where n is skipCount
    public (List<DirectedNode> borderCoordinates, List<GraphNode> extraVoronoiPoints)
        GetRiverMouthCandidates(float[,] heightmap, int skipCount, KdTree<Point> allBorderPointKdTree)
    {
        int width = heightmap.GetLength(0);
        int height = heightmap.GetLength(1);
        List<DirectedNode> borderCoordinates = new List<DirectedNode>();
        List<GraphNode> extraVoronoiPoints = new List<GraphNode>();

        // Find the first border point by doing a quick scan through the middle of the heightmap
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

        // Crawl over boundary of the shape to find border coordinates, adding every skipCount coordinate to a list as DirectedNode river mouth objects
        Point currentPoint = firstBorderPoint;

        // A vector pointing towards the ocean for each node, currently just pointing in an approximate not very precise direction
        Vector2 oceanDirection = new Vector2(-1, 0);


        // saving the last skipCount borderPoints used for detecting concaves and adding additional nodes
        int trailingPointCount = skipCount;
        List<Point> trailingPoints = new List<Point>(trailingPointCount);

        // Used for the crawling algorithm
        var facing = new Point(0, -1); // up

        int i = 0;
        while (true)
        {
            /* Crawling algorithm, "wall" means land here.
                4 rules for crawling the outline of the shape:
                ↑ ▮   ▮▮         ▮
                | ▮    ↰▮   ⮣▮   ⮣▮ 
                1       2     3    4
                If there's a wall ahead on the right and no wall blocking the path, move forward one cell (1).
                If there's a wall ahead on the right and a wall blocking the path, add a point and turn left (2).
                Otherwise add a point and turn right (3) (also covers 4)
                Note that lower Y value means up
                Note that we do point math as if we are crawling in/on the wall, not on the space next to it
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
            }

            currentPoint.X = Math.Clamp(currentPoint.X, 0, width);
            currentPoint.Y = Math.Clamp(currentPoint.Y, 0, width);

            trailingPoints.Add(currentPoint);

            if (i % 8 == 0)
                allBorderPointKdTree.Insert(currentPoint);

            if (currentPoint == firstBorderPoint) break; // Stop when we are back where we started

            // check if we need to add another potential river point for convex parts of the graph since we want higher chance of rivers here
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
                oceanDirection /= oceanDirection.Length(); // normalize
                borderCoordinates.Add(new DirectedNode(currentPoint, NodeType.Border, null, oceanDirection,
                    _random.Next(_parameters.MaxNodePriority, _parameters.MaxNodePriority)));
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

    //Check if we want to add any additional additional nodes before current to cover concave shape
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

    // Continue to add more river nodes based on the river-mouths until we have created a full river graph
    public (List<DirectedNode>, List<RiverEdge>) ExpandRiverGrid(
        float[,] heightmap,
        List<DirectedNode> riverMouthCandidates, 
        int stepIndex // Step value is used for showing step-by-step progress of the river graph building, default to 1 for not showing any intermediate steps
        )
    {
        List<DirectedNode> riverGrid = new List<DirectedNode>();
        List<RiverEdge> riverEdges = new List<RiverEdge>();

        PriorityQueue<DirectedNode, int> riverExpandQueue = new();

        riverMouthCandidates.ForEach(n =>
        {
            riverExpandQueue.Enqueue(n, _parameters.MaxNodePriority - n.Priority);
            riverGrid.Add(n);
        });

        // Keep taking the highest priority node from the queue and try to expand that river node further
        for (int i = 0; i < stepIndex; i++)
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
    
    // Try to expand a rivernode further
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
            // asymmetric branches has the main river continue with the same width/priority and a smaller offshoot
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


        double angleBetweenBranches = RandomInRange(_parameters.MinAngleBetweenBranches, 170);

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
                    Point.GetVectorBetweenPoints(newPoint, drainNode.Point), Math.Max(1, branch1Priority)));
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
            double randomDistance2 = _parameters.InterNodeDist +
                                      _random.NextDouble() * _parameters.InterNodeDistVar * 2 -
                                      _parameters.InterNodeDistVar;


            // var minAngle2 = angle1 + _parameters.MinAngleBetweenBranches;


            for (int i = 0; i < _parameters.NodeExpansionMaxTries * 3; i++)
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
                        Point.GetVectorBetweenPoints(newPoint, drainNode.Point), branch2Priority));
                    break;
                }

                randomDistance2 = _parameters.InterNodeDist + _random.NextDouble() * _parameters.InterNodeDistVar * 2 -
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


        // Points should be far enough from other rivernodes and should fill certain height requirements
        bool ValidateNewPoint(Point p, List<RiverEdge> edgesSoFar)
        {
            bool valid = false;
            var x = p.X;
            var y = p.Y;

            bool validCoordinate = x >= 0 && x < width && y >= 0 && y < width;
            float pHeight = validCoordinate ? heightmap[x, y] : 0;
            if (validCoordinate &&
                pHeight >= _parameters.NodeExpansionMinHeight &&
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

    public void FindFlowRateForAllNodes(float[,] heightmap, KdTree<GraphNode> graphNodeKdTree, int width,
        List<DirectedNode> allNodes)
    {
        // Iterate through the grid points
        for (int y = 0; y < width; y++)
        {
            for (int x = 0; x < width; x++)
            {
                if (heightmap[x, y] == 0) continue;

                // Add this pixel to the area of the closest rivernode (FlowRate is used for area at this point)
                List<GraphNode> nearestNeighbors = graphNodeKdTree.FindClosest(x, y, 1);
                ((DirectedNode)nearestNeighbors.First()).FlowRate++;
            }
        }

        // flowrate = 0.42 · A^0.69 formula from hydrology paper. Before this point FlowRate holds the number of pixels for which the closest node is this
        allNodes.ForEach(n => n.FlowRate = 0.42 * Math.Pow(n.FlowRate, 0.69));

        allNodes.ForEach(n => SetCumulativeFlowRate(n));
    }


    // For a given river-node, find the combined flow amount of all its children and assign it to node.FlowRate (before this node.FlowRate is assumed to hold a number representing the area that drains into that node)
    private double SetCumulativeFlowRate(DirectedNode node)
    {
        if (node.Children.Count > 0)
            node.FlowRate = node.Children.Sum(n => SetCumulativeFlowRate((DirectedNode)n));
        return node.FlowRate;
    }

    // For a given river-node, set the height of each node to a random number that is greater than it's parent. This is meant to support generating a realistic heightmap FROM the river nodes but currently we simply put the river nodes ON a realistic heightmap tso this height calculation is not currently useful
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
}