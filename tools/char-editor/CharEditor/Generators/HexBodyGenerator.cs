using System.Windows.Media;
using System.Windows.Media.Media3D;
using CharEditor.Models;

namespace CharEditor.Generators;

/// <summary>
/// Generates a voxel-style humanoid from hex prisms.
/// Matches the hex-planets aesthetic: blocky, hex-shaped body parts.
/// </summary>
public static class HexBodyGenerator
{
    // Hex prism: 6-sided polygon extruded vertically
    private static MeshGeometry3D CreateHexPrism(
        Point3D center, double radius, double height,
        double widthScale = 1.0, double depthScale = 1.0)
    {
        var mesh = new MeshGeometry3D();

        double halfH = height / 2.0;
        int sides = 6;

        // Generate top and bottom hex vertices
        var topVerts = new Point3D[sides];
        var botVerts = new Point3D[sides];

        for (int i = 0; i < sides; i++)
        {
            double angle = Math.PI / 3.0 * i + Math.PI / 6.0; // flat-top hex
            double x = center.X + radius * Math.Cos(angle) * widthScale;
            double z = center.Z + radius * Math.Sin(angle) * depthScale;
            topVerts[i] = new Point3D(x, center.Y + halfH, z);
            botVerts[i] = new Point3D(x, center.Y - halfH, z);
        }

        // Top face (fan from center) — CCW winding when viewed from above
        var topCenter = new Point3D(center.X, center.Y + halfH, center.Z);
        for (int i = 0; i < sides; i++)
        {
            int next = (i + 1) % sides;
            AddTriangle(mesh, topCenter, topVerts[next], topVerts[i]);
        }

        // Bottom face — CCW winding when viewed from below
        var botCenter = new Point3D(center.X, center.Y - halfH, center.Z);
        for (int i = 0; i < sides; i++)
        {
            int next = (i + 1) % sides;
            AddTriangle(mesh, botCenter, botVerts[i], botVerts[next]);
        }

        // Side faces (quads as 2 triangles each) — CCW winding when viewed from outside
        for (int i = 0; i < sides; i++)
        {
            int next = (i + 1) % sides;
            AddTriangle(mesh, topVerts[i], topVerts[next], botVerts[next]);
            AddTriangle(mesh, topVerts[i], botVerts[next], botVerts[i]);
        }

        return mesh;
    }

    /// <summary>
    /// Hex prism extruded along Z axis (facing forward), for eyeball-like shapes.
    /// The hex face is in the XY plane, extruded in +Z/-Z.
    /// </summary>
    private static MeshGeometry3D CreateHexPrismZ(
        Point3D center, double radius, double depth,
        double widthScale = 1.0, double heightScale = 1.0)
    {
        var mesh = new MeshGeometry3D();

        double halfD = depth / 2.0;
        int sides = 6;

        var frontVerts = new Point3D[sides];
        var backVerts = new Point3D[sides];

        for (int i = 0; i < sides; i++)
        {
            double angle = Math.PI / 3.0 * i + Math.PI / 6.0; // flat-top hex
            double x = center.X + radius * Math.Cos(angle) * widthScale;
            double y = center.Y + radius * Math.Sin(angle) * heightScale;
            frontVerts[i] = new Point3D(x, y, center.Z + halfD);
            backVerts[i] = new Point3D(x, y, center.Z - halfD);
        }

        // Front face (+Z) — CCW when viewed from front
        var frontCenter = new Point3D(center.X, center.Y, center.Z + halfD);
        for (int i = 0; i < sides; i++)
        {
            int next = (i + 1) % sides;
            AddTriangle(mesh, frontCenter, frontVerts[i], frontVerts[next]);
        }

        // Back face (-Z) — CCW when viewed from back
        var backCenter = new Point3D(center.X, center.Y, center.Z - halfD);
        for (int i = 0; i < sides; i++)
        {
            int next = (i + 1) % sides;
            AddTriangle(mesh, backCenter, backVerts[next], backVerts[i]);
        }

        // Side faces — CCW when viewed from outside
        for (int i = 0; i < sides; i++)
        {
            int next = (i + 1) % sides;
            AddTriangle(mesh, frontVerts[next], frontVerts[i], backVerts[i]);
            AddTriangle(mesh, frontVerts[next], backVerts[i], backVerts[next]);
        }

        return mesh;
    }

    private static void AddTriangle(MeshGeometry3D mesh, Point3D p0, Point3D p1, Point3D p2)
    {
        int idx = mesh.Positions.Count;
        mesh.Positions.Add(p0);
        mesh.Positions.Add(p1);
        mesh.Positions.Add(p2);

        // Flat shading normal
        var v1 = p1 - p0;
        var v2 = p2 - p0;
        var normal = Vector3D.CrossProduct(v1, v2);
        normal.Normalize();
        mesh.Normals.Add(normal);
        mesh.Normals.Add(normal);
        mesh.Normals.Add(normal);

        mesh.TriangleIndices.Add(idx);
        mesh.TriangleIndices.Add(idx + 1);
        mesh.TriangleIndices.Add(idx + 2);
    }

    /// <summary>
    /// Merge multiple meshes into one.
    /// </summary>
    private static MeshGeometry3D Merge(params MeshGeometry3D[] meshes)
    {
        var result = new MeshGeometry3D();
        foreach (var m in meshes)
        {
            int offset = result.Positions.Count;
            foreach (var p in m.Positions) result.Positions.Add(p);
            foreach (var n in m.Normals) result.Normals.Add(n);
            foreach (var i in m.TriangleIndices) result.TriangleIndices.Add(i + offset);
        }
        return result;
    }

    /// <summary>
    /// Generate the full humanoid body from appearance parameters.
    /// Returns (primaryMesh, accentMesh) for two-tone coloring.
    /// </summary>
    public static (MeshGeometry3D primary, MeshGeometry3D accent, MeshGeometry3D eyes) Generate(AgentAppearance app)
    {
        double scale = app.Height;
        double stocky = 0.5 + app.Stockiness * 0.5; // 0.5 - 1.0

        double hexR = 0.15; // base hex radius

        // Torso — chunky default
        double torsoW = stocky * app.TorsoWidth;
        double torsoR = hexR * (0.8 + torsoW * 0.8); // wider base
        double torsoH = 0.4 * scale;
        double torsoY = 0.25 * scale;

        // Pelvis — wide flat hex at bottom of torso
        double pelvisR = torsoR * 1.05;
        double pelvisH = 0.06 * scale;
        double pelvisY = torsoY - torsoH / 2.0 - pelvisH / 2.0;

        // Neck
        double neckR = hexR * 0.5;
        double neckH = 0.06 * scale;
        double neckY = torsoY + torsoH / 2.0 + neckH / 2.0;

        // Head
        double headR = hexR * 1.1 * app.HeadScale;
        double headH = 0.22 * scale * app.HeadScale;
        double headY = neckY + neckH / 2.0 + headH / 2.0;

        // Shoulder pads — two hex blocks sitting on top corners of torso
        double padR = hexR * 0.65;
        double padH = 0.1 * scale;
        double padSpacing = torsoR * 0.85;
        double padY = torsoY + torsoH / 2.0 - padH * 0.3;

        // Arms (hang from shoulder pads)
        double armR = hexR * 0.55;
        double armH = 0.38 * scale * app.ArmLength;
        double armSpacing = padSpacing + padR * 0.3;
        double armY = padY - padH / 2.0 - armH / 2.0;

        // Legs
        double legR = hexR * 0.55;
        double legH = 0.4 * scale * app.LegLength;
        double legSpacing = hexR * stocky * 0.75;
        double legY = pelvisY - pelvisH / 2.0 - legH / 2.0;

        // Anchor feet to Y=0
        double feetY = legY - legH / 2.0;
        double yOff = -feetY;
        torsoY += yOff; pelvisY += yOff; neckY += yOff; headY += yOff;
        padY += yOff; armY += yOff; legY += yOff;

        // ---- Primary color: torso, legs, pelvis, shoulder pads ----
        var torso = CreateHexPrism(
            new Point3D(0, torsoY, 0),
            torsoR, torsoH, stocky, 0.8);

        var leftShoulderPad = CreateHexPrism(
            new Point3D(-padSpacing, padY, 0),
            padR, padH, 0.8, 0.8);

        var rightShoulderPad = CreateHexPrism(
            new Point3D(padSpacing, padY, 0),
            padR, padH, 0.8, 0.8);

        var pelvis = CreateHexPrism(
            new Point3D(0, pelvisY, 0),
            pelvisR, pelvisH, stocky, 0.7);

        var leftLeg = CreateHexPrism(
            new Point3D(-legSpacing, legY, 0),
            legR, legH, stocky * 0.65, 0.65);

        var rightLeg = CreateHexPrism(
            new Point3D(legSpacing, legY, 0),
            legR, legH, stocky * 0.65, 0.65);

        var primaryMesh = Merge(torso, leftShoulderPad, rightShoulderPad, pelvis, leftLeg, rightLeg);

        // ---- Accent color: head, neck, arms ----
        var neck = CreateHexPrism(
            new Point3D(0, neckY, 0),
            neckR, neckH, 1.0, 1.0);

        var head = CreateHexPrism(
            new Point3D(0, headY, 0),
            headR, headH, 1.0, 1.0);

        var leftArm = CreateHexPrism(
            new Point3D(-armSpacing, armY, 0),
            armR, armH, 0.55, 0.55);

        var rightArm = CreateHexPrism(
            new Point3D(armSpacing, armY, 0),
            armR, armH, 0.55, 0.55);

        var accentMesh = Merge(neck, head, leftArm, rightArm);

        // ---- Eyes: small hex prisms extruded along Z (facing forward) ----
        double eyeR = headR * 0.22;
        double eyeDepth = headR * 0.3; // protrude slightly from face
        double eyeSpacing = headR * 0.45;
        double eyeY = headY + headH * 0.08; // slightly above head center
        double eyeZ = headR * 0.7; // on the front face of the head

        var leftEye = CreateHexPrismZ(
            new Point3D(-eyeSpacing, eyeY, eyeZ),
            eyeR, eyeDepth);

        var rightEye = CreateHexPrismZ(
            new Point3D(eyeSpacing, eyeY, eyeZ),
            eyeR, eyeDepth);

        var eyesMesh = Merge(leftEye, rightEye);

        return (primaryMesh, accentMesh, eyesMesh);
    }

    /// <summary>
    /// Individual body parts for animation. Each mesh is generated at origin
    /// with a pivot point for rotation.
    /// </summary>
    public record BodyParts
    {
        public MeshGeometry3D Torso = null!;       // includes shoulder pads + pelvis
        public MeshGeometry3D Head = null!;         // head + neck + eyes
        public MeshGeometry3D HeadAccent = null!;   // head accent color
        public MeshGeometry3D Eyes = null!;
        public MeshGeometry3D LeftArm = null!;
        public MeshGeometry3D RightArm = null!;
        public MeshGeometry3D LeftLeg = null!;
        public MeshGeometry3D RightLeg = null!;

        // Pivot points (world Y positions for rotation centers)
        public double HeadPivotY;       // neck base
        public double ArmPivotY;        // shoulder top
        public double ArmSpacing;       // X offset for arms
        public double LegPivotY;        // hip joint
        public double LegSpacing;       // X offset for legs
        public double BodyCenterY;      // torso center
    }

    public static BodyParts GenerateParts(AgentAppearance app)
    {
        double scale = app.Height;
        double stocky = 0.5 + app.Stockiness * 0.5;
        double hexR = 0.15;

        double torsoW = stocky * app.TorsoWidth;
        double torsoR = hexR * (0.8 + torsoW * 0.8);
        double torsoH = 0.4 * scale;
        double torsoY = 0.25 * scale;

        double pelvisR = torsoR * 1.05;
        double pelvisH = 0.06 * scale;
        double pelvisY = torsoY - torsoH / 2.0 - pelvisH / 2.0;

        double neckR = hexR * 0.5;
        double neckH = 0.06 * scale;
        double neckY = torsoY + torsoH / 2.0 + neckH / 2.0;

        double headR = hexR * 1.1 * app.HeadScale;
        double headH = 0.22 * scale * app.HeadScale;
        double headY = neckY + neckH / 2.0 + headH / 2.0;

        double padR = hexR * 0.65;
        double padH = 0.1 * scale;
        double padSpacing = torsoR * 0.85;
        double padY = torsoY + torsoH / 2.0 - padH * 0.3;

        double armR = hexR * 0.55;
        double armH = 0.38 * scale * app.ArmLength;
        double armSpacing = padSpacing + padR * 0.3;
        double armY = padY - padH / 2.0 - armH / 2.0;

        double legR = hexR * 0.55;
        double legH = 0.4 * scale * app.LegLength;
        double legSpacing = hexR * stocky * 0.75;
        double legY = pelvisY - pelvisH / 2.0 - legH / 2.0;

        // Anchor feet to Y=0: shift everything up by the bottom of the legs
        double feetY = legY - legH / 2.0;
        double yOff = -feetY; // shift up so feet sit at Y=0

        torsoY += yOff;
        pelvisY += yOff;
        neckY += yOff;
        headY += yOff;
        padY += yOff;
        armY += yOff;
        legY += yOff;

        double eyeR = headR * 0.22;
        double eyeDepth = headR * 0.3;
        double eyeSpacing = headR * 0.45;
        double eyeY = headY + headH * 0.08;
        double eyeZ = headR * 0.7;

        return new BodyParts
        {
            Torso = Merge(
                CreateHexPrism(new Point3D(0, torsoY, 0), torsoR, torsoH, stocky, 0.8),
                CreateHexPrism(new Point3D(-padSpacing, padY, 0), padR, padH, 0.8, 0.8),
                CreateHexPrism(new Point3D(padSpacing, padY, 0), padR, padH, 0.8, 0.8),
                CreateHexPrism(new Point3D(0, pelvisY, 0), pelvisR, pelvisH, stocky, 0.7)),
            HeadAccent = Merge(
                CreateHexPrism(new Point3D(0, neckY, 0), neckR, neckH, 1.0, 1.0),
                CreateHexPrism(new Point3D(0, headY, 0), headR, headH, 1.0, 1.0)),
            Eyes = Merge(
                CreateHexPrismZ(new Point3D(-eyeSpacing, eyeY, eyeZ), eyeR, eyeDepth),
                CreateHexPrismZ(new Point3D(eyeSpacing, eyeY, eyeZ), eyeR, eyeDepth)),
            LeftArm = CreateHexPrism(new Point3D(-armSpacing, armY, 0), armR, armH, 0.55, 0.55),
            RightArm = CreateHexPrism(new Point3D(armSpacing, armY, 0), armR, armH, 0.55, 0.55),
            LeftLeg = CreateHexPrism(new Point3D(-legSpacing, legY, 0), legR, legH, stocky * 0.65, 0.65),
            RightLeg = CreateHexPrism(new Point3D(legSpacing, legY, 0), legR, legH, stocky * 0.65, 0.65),

            HeadPivotY = neckY - neckH / 2.0,
            ArmPivotY = padY,
            ArmSpacing = armSpacing,
            LegPivotY = pelvisY - pelvisH / 2.0,
            LegSpacing = legSpacing,
            BodyCenterY = torsoY,
        };
    }

    public static Color ToColor(double r, double g, double b)
    {
        return Color.FromRgb(
            (byte)(Math.Clamp(r, 0, 1) * 255),
            (byte)(Math.Clamp(g, 0, 1) * 255),
            (byte)(Math.Clamp(b, 0, 1) * 255));
    }
}
