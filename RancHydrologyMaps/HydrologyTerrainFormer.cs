using System.Numerics;
using RancHydrologyMaps;

namespace HydrologyMaps;

public class HydrologyTerrainFormer
{
    private static FastNoiseLite _noiseGenerator = new();

    public HydrologyTerrainFormer()
    {
        _noiseGenerator.SetNoiseType(FastNoiseLite.NoiseType.Perlin);
        _noiseGenerator.SetFrequency(0.013f);
        _noiseGenerator.SetFractalGain(1);
    }


    public static void CarveRivers(float[,] heightMap, List<DirectedNode> riverRoots, float baseRiverWidth,
        float riverDepth, float falloff, float flowRateInfluence)
    {
        float[,] accumulatedDepthMap = new float[heightMap.GetLength(0), heightMap.GetLength(1)];

        foreach (var root in riverRoots)
        {
            AccumulateRiverSegmentsBFS(accumulatedDepthMap, root, baseRiverWidth, riverDepth, falloff,
                flowRateInfluence);
        }

        ApplyAccumulatedDepthMap(heightMap, accumulatedDepthMap);
    }

    private static void AccumulateRiverSegmentsBFS(float[,] accumulatedDepthMap, DirectedNode root,
        float baseRiverWidth, float riverDepth, float falloff, float flowRateInfluence)
    {
        Queue<DirectedNode> queue = new Queue<DirectedNode>();
        queue.Enqueue(root);

        while (queue.Count > 0)
        {
            DirectedNode currentNode = queue.Dequeue();

            foreach (var childGraphNode in currentNode.Children)
            {
                var childNode = (DirectedNode)childGraphNode;
                queue.Enqueue(childNode);
            }

            if (currentNode.Parent != null)
            {
                DirectedNode parentNode = (DirectedNode)currentNode.Parent;

                List<Point> fullPath = new List<Point>
                    { parentNode.Point, currentNode.Point }; //BuildFullPath(parentNode);
                List<Point> smoothedPath = GenerateSmoothedPath(parentNode.Point, currentNode.Point, new Random());

                float startFlowRate = (float)parentNode.FlowRate / 512;
                float endFlowRate = (float)currentNode.FlowRate / 512;

                AccumulateRiverSegment(accumulatedDepthMap, smoothedPath, startFlowRate, endFlowRate, baseRiverWidth,
                    riverDepth, falloff, flowRateInfluence);
            }

            // If the current node is a leaf node, process it here.
            //  if (currentNode.Children.Count == 0 && currentNode.Parent != null)
            //  {
            //      DirectedNode parentNode = (DirectedNode)currentNode.Parent;
//
            //      List<Point> fullPath = new List<Point> { parentNode.Point, currentNode.Point };
            //      List<Point> smoothedPath =  GenerateSmoothedPath(parentNode.Point, currentNode.Point, new Random());
//
            //      float startFlowRate = (float)parentNode.FlowRate / 512;
            //      float endFlowRate = (float)currentNode.FlowRate / 512;
//
            //      AccumulateRiverSegment(accumulatedDepthMap, smoothedPath, startFlowRate, endFlowRate, baseRiverWidth,
            //          riverDepth, falloff, flowRateInfluence);
            //  }
        }
    }

    private static List<Point> BuildFullPath(DirectedNode currentNode)
    {
        List<DirectedNode> fullPath = new List<DirectedNode> { currentNode };

        DirectedNode parentNode = (DirectedNode)currentNode.Parent;
        while (parentNode != null)
        {
            fullPath.Insert(0, parentNode);
            parentNode = (DirectedNode)parentNode.Parent;
        }

        return fullPath.Select(node => node.Point).ToList();
    }

    private static List<Point> GenerateSmoothedPath(Point start, Point end, Random random)
    {
        List<Point> path = new List<Point> { start, end };
        List<Point> smoothedPath = new List<Point>();

        Vector2 direction = (end - start).ToVector();
        direction /= direction.Length();

        // Add random control points to the path for twists and turns
        int numControlPoints = random.Next(1, 4);
        for (int i = 0; i < numControlPoints; i++)
        {
            float t = (float)i / (numControlPoints + 1);
            Point controlPoint = Point.Lerp(start, end, t);

            // Create a random angle based on the current direction (between -45 and 45 degrees)
            float angle = random.Next(-90, 90) * (float)Math.PI / 180f;

            // Calculate the new direction with the random angle applied
            Vector2 newDirection = new Vector2(
                direction.X * (float)Math.Cos(angle) - direction.Y * (float)Math.Sin(angle),
                direction.X * (float)Math.Sin(angle) + direction.Y * (float)Math.Cos(angle)
            );

            // Calculate the random offset using the new direction
            Vector2 randomOffset = newDirection * random.Next(1, 5);

            controlPoint += Point.FromVector(randomOffset);
            path.Insert(1 + i, controlPoint);
        }

        // Smooth the path using Catmull-Rom splines
        for (float t = 0; t < path.Count - 1; t += 0.5f)
        {
            int index = (int)Math.Floor(t);
            Point p0 = path[Math.Clamp(index - 1, 0, path.Count - 1)];
            Point p1 = path[index];
            Point p2 = path[index + 1];
            Point p3 = path[Math.Clamp(index + 2, 0, path.Count - 1)];

            Vector2 point = CatmullRomSpline(p0.ToVector(), p1.ToVector(), p2.ToVector(), p3.ToVector(), t - index);
            smoothedPath.Add(Point.FromVector(point));
        }

        smoothedPath.Add(path[^1]);
        return smoothedPath;
    }

    // Calculate the Catmull-Rom spline for a set of control points and a parameter t
    private static Vector2 CatmullRomSpline(Vector2 p0, Vector2 p1, Vector2 p2, Vector2 p3, float t)
    {
        Vector2 a = 2f * p1;
        Vector2 b = p2 - p0;
        Vector2 c = 2f * p0 - 5f * p1 + 4f * p2 - p3;
        Vector2 d = -p0 + 3f * p1 - 3f *
            p2 + p3;

        return 0.5f * (a + (b * t) + (c * t * t) + (d * t * t * t));
    }

    private static void AccumulateRiverSegment(float[,] accumulatedDepthMap, List<Point> smoothedPath,
        float startFlowRate, float endFlowRate, float baseRiverWidth, float riverDepth, float fallOffSharpness,
        float flowRateInfluence)
    {
        for (int i = 0; i < smoothedPath.Count - 1; i++)
        {
            Point start = smoothedPath[i];
            Point end = smoothedPath[i + 1];

            float t = (float)i / (smoothedPath.Count - 1);


            // formula is baseRiverWidth + sigmoid formula centered on 1
            float startFlowRateFactor =
                (float)(1 + flowRateInfluence * ((1 / (1 + Math.Pow(Math.E, -startFlowRate * 10)) - 0.5)));
            float endFlowRateFactor =
                (float)(1 + flowRateInfluence * ((1 / (1 + Math.Pow(Math.E, -endFlowRate * 10)) - 0.5)));
            float startWidth = baseRiverWidth * startFlowRateFactor;
            float endWidth = baseRiverWidth * endFlowRateFactor;

            float
                startDepth =
                    riverDepth; //*(float)Math.Pow(startFlowRateFactor, 0.6f);// 0.5f * riverDepth + 0.02f * riverDepth * startFlowRate;
            float
                endDepth = riverDepth; //*(float)Math.Pow(endFlowRateFactor, 0.6f);//0.5f *riverDepth + 0.02f * riverDepth* endFlowRate;


            CarveRiverSegment(accumulatedDepthMap, start, end, startDepth, endDepth, startWidth, endWidth,
                fallOffSharpness);
        }
    }

    private static void CarveRiverSegment(float[,] accumulatedDepthMap, Point start, Point end, float startDepth,
        float endDepth, float startWidth, float endWidth, float fallOffSharpness)
    {
        // Calculate the line equation between start and end points
        Vector2 direction = (end - start).ToVector();
        int numSteps = (int)direction.Length();
        direction /= numSteps;

        // Iterate over the points along the line between start and end points
        for (int i = 0; i <= numSteps; i++)
        {
            float t = (float)i / numSteps;

            Point position = start + Point.FromVector(direction * i);
            float currentDepth = Lerp(startDepth, endDepth, t);
            float currentWidth = Lerp(startWidth, endWidth, t);

            int widthInt = (int)Math.Ceiling(currentWidth);

            // Iterate over the points in the width direction
            for (int xOffset = -widthInt; xOffset <= widthInt; xOffset++)
            {
                for (int yOffset = -widthInt; yOffset <= widthInt; yOffset++)
                {
                    Point offsetPosition = position + new Point(xOffset, yOffset);

                    if (!IsInBounds(offsetPosition, accumulatedDepthMap)) continue;

                    float distance = new Vector2(xOffset, yOffset).Length();

                    if (distance <= currentWidth)
                    {
                        float scaledDistance = distance / currentWidth;
                        //float depthFalloff = (float)(1 - (3 * Math.Pow(scaledDistance, fallOffSharpness) - 2 * Math.Pow(scaledDistance, fallOffSharpness + 1)));
                        //    double depthFalloff = 1 / (1+Math.Pow(Math.E, fallOffSharpness*(distance-(currentWidth*0.8f))));

                        //double depthFalloff = 1-(distance/currentWidth);

                        double depthFalloff = 1 / (1 + Math.Pow(Math.E,
                            fallOffSharpness * (distance - (currentWidth * (0.7)))));
                        //depthFalloff = 1-scaledDistance;
                        float adjustedDepth = currentDepth * (float)depthFalloff;
                        accumulatedDepthMap[offsetPosition.X, offsetPosition.Y] =
                            Math.Max(accumulatedDepthMap[offsetPosition.X, offsetPosition.Y], adjustedDepth);
                    }
                }
            }
        }
    }

    private static bool IsInBounds(Point position, float[,] map)
    {
        int width = map.GetLength(0);
        int height = map.GetLength(1);

        return position.X >= 0 && position.X < width && position.Y >= 0 && position.Y < height;
    }

    private static void ApplyAccumulatedDepthMap(float[,] heightMap, float[,] accumulatedDepthMap)
    {
        for (int x = 0; x < heightMap.GetLength(0); x++)
        {
            for (int y = 0; y < heightMap.GetLength(1); y++)
            {
                heightMap[x, y] = Math.Max(heightMap[x, y] - accumulatedDepthMap[x, y], 0);
            }
        }
    }

    private static float Lerp(float start, float end, float t)
    {
        return start + (end - start) * t;
    }

    public static void InterpolateHeightMapSimple(float[,] heightmap, KdTree<Point> borderGraphNodeKdTree,
        int width, int numNeighbors, float maxSlope)
    {
        for (int y = 0; y < width; y++)
        {
            for (int x = 0; x < width; x++)
            {
                if (heightmap[x, y] == 0) continue;
                var closestBorderPoint = borderGraphNodeKdTree.FindClosest(x, y, 1).First();
                float dist = Math.Min(1,(float)Point.Distance(new Point(x, y), new Point(closestBorderPoint.X, closestBorderPoint.Y))/180);
           //     heightmap[x, y] = heightmap[x, y] *dist * Math.Min(1, 0.7f + _noiseGenerator.GetNoise(x * 3, y * 3) * 0.3f);
                heightmap[x, y] = dist*Math.Min(1, 0.7f + _noiseGenerator.GetNoise(x, y) * 0.3f);
            }
        }
    }


    // ========== Base heightmap formation ===========
    public static void InterpolateHeightMap(float[,] heightmap, KdTree<GraphNode> graphNodeKdTree, List<DirectedNode> allNodes,
        int width, int height, int numNeighbors, float maxSlope)
    {
        // Iterate through the grid points
        for (int y = 0; y < height; y++)
        {
            for (int x = 0; x < width; x++)
            {
                if (heightmap[x, y] == 0) continue;

                // Query k-d tree for k nearest neighbors
                List<GraphNode> nearestNeighbors = graphNodeKdTree.FindClosest(x, y, numNeighbors);

                // Calculate the interpolated height using IDW or other methods
                float interpolatedHeight = 0;
                double totalWeight = 0;
                float minHeight = nearestNeighbors.Max(x => x.Height);
                double distanceSum = 0;
                foreach (var neighbor in nearestNeighbors)
                {
                    double distance = Point.Distance(new Point(x, y), neighbor.Point);
                    distanceSum += distance;
                    if (distance == 0) distance = 1;

                    double weight = 1 / distance; // You can raise this to a power to control smoothness

                    totalWeight += weight;
                    interpolatedHeight += (float)weight * neighbor.Height;
                }

                float minBaseHeight = 0.1f; // Math.Min(0.3f, (float)Math.Pow(distanceSum, 1.1) * 0.001f);

                interpolatedHeight = Math.Max(minBaseHeight, minHeight + (float)(interpolatedHeight / totalWeight));

                // Adjust the interpolated height based on the closest node and maxSlope if needed
                // ...

                heightmap[x, y] = interpolatedHeight;
            }
        }

        // Second pass: adjust the height values based on the maxSlope constraint
        for (int y = 0; y < height; y++)
        {
            for (int x = 0; x < width; x++)
            {
                float currentHeight = heightmap[x, y];

                if (currentHeight == 0) continue;

                // Iterate through the 8 neighboring pixels
                for (int dy = -1; dy <= 1; dy++)
                {
                    for (int dx = -1; dx <= 1; dx++)
                    {
                        if (dx == 0 && dy == 0) continue;

                        int nx = x + dx;
                        int ny = y + dy;

                        // Check if the neighboring pixel is within the bounds of the heightmap
                        if (nx >= 0 && nx < width && ny >= 0 && ny < height)
                        {
                            float neighborHeight = heightmap[nx, ny];
                            float heightDifference = Math.Abs(currentHeight - neighborHeight);

                            // Check if the height difference exceeds the maxSlope constraint
                            if (heightDifference > maxSlope)
                            {
                                float adjustedHeightDifference = maxSlope * (heightDifference / maxSlope);
                                if (currentHeight > neighborHeight)
                                {
                                    heightmap[nx, ny] = currentHeight - adjustedHeightDifference;
                                }
                                else
                                {
                                    heightmap[x, y] = neighborHeight - adjustedHeightDifference;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    public static void GetBaseHeightmapMask(float[,] heightmap, int width, Random random, HydrologyParameters parameters)
    {
        FastNoiseLite noiseGenerator = new FastNoiseLite();
        noiseGenerator.SetNoiseType(FastNoiseLite.NoiseType.Perlin);
        noiseGenerator.SetSeed(random.Next());
        noiseGenerator.SetFrequency(0.01f);
        noiseGenerator.SetFractalGain(0);


        int borderSize = 5;

        float beachWidth = parameters.BeachWidth;

        for (int y = 0; y < width; y++)
        {
            for (int x = 0; x < width; x++)
            {
                float radialGradient = GetIslandFalloffFactor(x, y, width, parameters.IslandShapeChaos);

                float perlinValue =
                    (noiseGenerator.GetNoise(x, y) + 1) / 2 * radialGradient; // Normalize Perlin noise to the range [0, 1]

           //     if (perlinValue < 0.2)
           //     {
                    if (perlinValue < 0.8)
                    {
                        perlinValue = 0;
                    }
                    else
                    {
                        perlinValue = 1;
                    }
           //         else // if (x > beachWidth && x < width - beachWidth  && y > beachWidth && y < width - beachWidth)
           //         {
                /*        bool isNearOcean = false;
                        for (int yo = -(int)beachWidth; yo < beachWidth; yo++)
                        {
                            for (int xo = -(int)beachWidth; xo < beachWidth; xo++)
                            {
                                float offsetRadialGradient = GetIslandFalloffFactor(x + xo, y + yo, width, parameters.IslandFalloffWidth);
                                float offSetPerlinValue =
                                    (noiseGenerator.GetNoise(x + xo, y + yo) + 1f) * offsetRadialGradient;

                                if (offSetPerlinValue < 0.15)
                                {
                                   // isNearOcean = true;
                                    break;
                                }
                            }

                            if (isNearOcean) break;
                        }

                        if (isNearOcean) perlinValue = 0.1f;
           //         }*/
            //    }


                //  else if (perlinValue < 0.165) perlinValue = 0.15f;
                //else perlinValue = 0.5f;
                heightmap[x, y] = Math.Min(1, perlinValue);
            }

        }
        
        // maxDiameterOfArtifacts, higher value is more certain to remove artifacts but increases computation time of the algorithm
        HeightmapCleaner.RemoveArtifacts(heightmap, maxDiameterOfArtifacts: 100);
    }
    
   //static float GetRadialGradient(int x, int y, int width, double falloffWidthFromEdge)
   //{
   //    float center = width / 2f;
   //    
   //    float distanceX = x - center;
   //    float distanceY = y - center;
   //    float distanceFromCenter = (float)Math.Sqrt(distanceX * distanceX + distanceY * distanceY);
   //    float distanceFromEdge = width -  distanceFromCenter;
   //    float falloffFactor = Math.Min(1f, distanceFromEdge / (width / 2.5f));
   //    return distanceFromEdge > falloffWidthFromEdge ? 1 : falloffFactor;
   //}
    
    
    static float GetIslandFalloffFactor(int x, int y, int width, float islandShapeChaos)
    {
        float center = width / 2f;
        
        float distanceX = x - center;
        float distanceY = y - center;
        float distanceToCenter = (float)Math.Sqrt(distanceX * distanceX + distanceY * distanceY);
        float distanceFromEdge = width / 2f -  distanceToCenter;
        const float minDistFromWall = 3;
        float chaos = 50 * islandShapeChaos;
        float factor = (distanceFromEdge - minDistFromWall + chaos) / (minDistFromWall + chaos);

        return factor; // (1f - Math.Min(1f, normalizedDistance)) / 2;
    }
}