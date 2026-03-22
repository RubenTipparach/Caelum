using System.Windows.Media.Media3D;

namespace CharEditor.Generators;

public enum AnimState
{
    Idle,
    Walk,
    Jump,
    Talk,
    PickUp,
    Place,
    Wave,
    Celebrate,
}

/// <summary>
/// Procedural animation via per-body-part transforms.
/// Returns offsets/rotations for each part given a time value.
/// </summary>
public static class AnimationController
{
    public record BodyPose
    {
        public Vector3D BodyOffset;
        public double BodyTilt;          // forward/back pitch in degrees
        public double HeadTilt;          // head nod degrees
        public double HeadTurn;          // head yaw degrees
        public double LeftArmPitch;      // forward/back swing degrees
        public double RightArmPitch;
        public double LeftArmRaise;      // sideways raise degrees
        public double RightArmRaise;
        public double LeftLegPitch;      // forward/back swing degrees
        public double RightLegPitch;
    }

    public static BodyPose Evaluate(AnimState state, double t)
    {
        return state switch
        {
            AnimState.Idle => EvalIdle(t),
            AnimState.Walk => EvalWalk(t),
            AnimState.Jump => EvalJump(t),
            AnimState.Talk => EvalTalk(t),
            AnimState.PickUp => EvalPickUp(t),
            AnimState.Place => EvalPlace(t),
            AnimState.Wave => EvalWave(t),
            AnimState.Celebrate => EvalCelebrate(t),
            _ => new BodyPose()
        };
    }

    private static double Sin(double t) => Math.Sin(t);
    private static double Cos(double t) => Math.Cos(t);

    private static BodyPose EvalIdle(double t)
    {
        // Gentle breathing bob
        double breathe = Sin(t * 2.0) * 0.003;
        return new BodyPose
        {
            BodyOffset = new Vector3D(0, breathe, 0),
            LeftArmPitch = Sin(t * 1.5) * 2,
            RightArmPitch = Sin(t * 1.5 + 0.5) * 2,
        };
    }

    private static BodyPose EvalWalk(double t)
    {
        double speed = 6.0;
        double stride = Sin(t * speed);
        double bob = Math.Abs(Sin(t * speed)) * 0.015;
        return new BodyPose
        {
            BodyOffset = new Vector3D(0, bob, 0),
            BodyTilt = 3,
            LeftLegPitch = stride * 30,
            RightLegPitch = -stride * 30,
            LeftArmPitch = -stride * 25,
            RightArmPitch = stride * 25,
            HeadTilt = Sin(t * speed * 2) * 2,
        };
    }

    private static BodyPose EvalJump(double t)
    {
        // Loop: crouch (0-0.3), launch (0.3-0.5), air (0.5-0.8), land (0.8-1.0)
        double phase = (t * 1.5) % 2.0;
        double yOff = 0, tilt = 0, armRaise = 0;

        if (phase < 0.5)
        {
            // Crouch
            double p = phase / 0.5;
            yOff = -p * 0.04;
            tilt = p * 8;
        }
        else if (phase < 1.0)
        {
            // Launch + air
            double p = (phase - 0.5) / 0.5;
            yOff = Sin(p * Math.PI) * 0.12;
            tilt = (1 - p) * 8;
            armRaise = Sin(p * Math.PI) * 40;
        }
        else if (phase < 1.5)
        {
            // Landing
            double p = (phase - 1.0) / 0.5;
            yOff = -Sin(p * Math.PI * 0.5) * 0.02;
            tilt = Sin(p * Math.PI) * 3;
        }

        return new BodyPose
        {
            BodyOffset = new Vector3D(0, yOff, 0),
            BodyTilt = tilt,
            LeftArmRaise = armRaise,
            RightArmRaise = armRaise,
            LeftLegPitch = tilt * 1.5,
            RightLegPitch = tilt * 1.5,
        };
    }

    private static BodyPose EvalTalk(double t)
    {
        // Head bobbing, slight body sway, small arm gestures
        return new BodyPose
        {
            HeadTilt = Sin(t * 4) * 8,
            HeadTurn = Sin(t * 2.3) * 10,
            BodyTilt = Sin(t * 1.8) * 2,
            LeftArmPitch = Sin(t * 3.5) * 10 + 5,
            RightArmPitch = Sin(t * 2.7 + 1) * 12,
            LeftArmRaise = Sin(t * 3.0) * 5,
        };
    }

    private static BodyPose EvalPickUp(double t)
    {
        // Cyclic: lean forward to reach ground, grab, stand back up
        double phase = (t * 0.8) % 3.0;
        double amount = 0;

        if (phase < 1.0)
        {
            // Lean down
            amount = Sin(phase * Math.PI * 0.5);
        }
        else if (phase < 1.5)
        {
            // Hold — grabbing
            amount = 1.0;
        }
        else if (phase < 2.5)
        {
            // Stand back up
            amount = 1.0 - (phase - 1.5) / 1.0;
        }

        // Positive BodyTilt = lean forward, negative arm pitch = arms swing forward
        return new BodyPose
        {
            BodyTilt = amount * 30,
            BodyOffset = new Vector3D(0, -amount * 0.04, 0),
            LeftArmPitch = -amount * 40,   // arms reach forward/down
            RightArmPitch = -amount * 40,
            LeftLegPitch = -amount * 10,    // slight knee bend back
            RightLegPitch = -amount * 10,
            HeadTilt = amount * 10,         // compensate — look at ground, not straight down
        };
    }

    private static BodyPose EvalPlace(double t)
    {
        // Arms forward and up, slight lean
        double phase = (t * 1.5) % 2.0;
        double reach = 0;

        if (phase < 1.0)
        {
            reach = Sin(phase * Math.PI) * 1.0;
        }

        return new BodyPose
        {
            BodyTilt = reach * -5,
            LeftArmPitch = -30 * reach,
            RightArmPitch = -30 * reach,
            LeftArmRaise = 15 * reach,
            RightArmRaise = 15 * reach,
            HeadTilt = -10 * reach,
        };
    }

    private static BodyPose EvalWave(double t)
    {
        double wave = Sin(t * 6) * 25;
        return new BodyPose
        {
            RightArmRaise = 140 + wave,  // side to side swing
            HeadTurn = 10,
            HeadTilt = Sin(t * 3) * 5,
            BodyOffset = new Vector3D(0, Sin(t * 3) * 0.005, 0),
        };
    }

    private static BodyPose EvalCelebrate(double t)
    {
        double bounce = Math.Abs(Sin(t * 5)) * 0.03;
        double armWiggle = Sin(t * 8) * 10;
        return new BodyPose
        {
            BodyOffset = new Vector3D(0, bounce, 0),
            LeftArmRaise = 75 + armWiggle,
            RightArmRaise = 75 - armWiggle,
            HeadTilt = Sin(t * 6) * 8,
            HeadTurn = Sin(t * 4) * 12,
            LeftLegPitch = Sin(t * 5) * 10,
            RightLegPitch = -Sin(t * 5) * 10,
        };
    }
}
