using System.Numerics;

namespace HydrologyMaps;

public class HydrologyTerrainFormer
{
    private static FastNoiseLite _noiseGenerator = new();

    public HydrologyTerrainFormer()
    {
        _noiseGenerator.SetNoiseType(FastNoiseLite.NoiseType.Perlin);
        _noiseGenerator.SetFrequency(0.003f);
        _noiseGenerator.SetFractalGain(0);
    }
    

    public static void CarveRivers(float[,] heightMap, List<DirectedNode> riverRoots, float baseRiverWidth,
        float riverDepth, float smoothness, float falloff, float flowRateInfluence)
    {
        float[,] accumulatedDepthMap = new float[heightMap.GetLength(0), heightMap.GetLength(1)];

        foreach (var root in riverRoots)
        {
            AccumulateRiverSegmentsBFS(accumulatedDepthMap, root, baseRiverWidth, riverDepth, smoothness, falloff, flowRateInfluence);
        }

        ApplyAccumulatedDepthMap(heightMap, accumulatedDepthMap);
    }

    private static void AccumulateRiverSegmentsBFS(float[,] accumulatedDepthMap, DirectedNode root,
        float baseRiverWidth, float riverDepth, float smoothness, float falloff, float flowRateInfluence)
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

                List<Point> fullPath = new List<Point> { parentNode.Point, currentNode.Point };//BuildFullPath(parentNode);
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
        float startFlowRate, float endFlowRate, float baseRiverWidth, float riverDepth, float fallOffSharpness, float flowRateInfluence)
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

            float startDepth =riverDepth;//*(float)Math.Pow(startFlowRateFactor, 0.6f);// 0.5f * riverDepth + 0.02f * riverDepth * startFlowRate;
            float endDepth =  riverDepth;//*(float)Math.Pow(endFlowRateFactor, 0.6f);//0.5f *riverDepth + 0.02f * riverDepth* endFlowRate;
            
            
            CarveRiverSegment(accumulatedDepthMap, start, end, startDepth, endDepth, startWidth, endWidth, fallOffSharpness);
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
                    
                        double depthFalloff = 1 / (1+Math.Pow(Math.E, fallOffSharpness*(distance-(currentWidth*(0.7)))));
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
}