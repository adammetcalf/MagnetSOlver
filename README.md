This is the magnetic tentacle solver

---
title: Tentacle Solver
---

```mermaid
classDiagram
    class Joint {
        -double a
        -double d
        -double theta
        -double alpha
        -double[][] Frame
        +Joint(theta, alpha, a, d)
        +UpdateAngles(Angles)
        +getFrame() double[][]
        -evaluateFrame()
    }
    class Tentacle{
        -Joints[] Joints
        -double[][] HGMs
        -double[] Jacobean
        -double[] Links
        -double[][] MagneticMoments
        -double MomentStrength
        +Tentacle(LinkLength,Angles[],Magnetisation)
        +UpdateAngles(Angles)
        +getHGMs() double[][]
        +getJacobean() double
        +getLinks() double[]
        +getMagneticMoments() double[][]
        -EvaluateHGM()
        -EvaluateJacobean()
        -EvaluateMagMoments()
        -GetCOM()
    }
    class World {
        +double Gravity
        +Tentacle Tentacle
        +double MagForceTorque
        +double PlotLength
        +double mu0
        +double threshold
        +double HomegeneousField
        +logical IncludeMultiPole
        +double x
        +double y
        +double z
        +double xIncr
        +double yIncr
        +double zIncr
        +double Bx
        +double By
        +double Bz
        +World(LinkLength, Angles, Magnetisation, HomogeneousField, MultipoleActive)
        +UpDateAngles(Angles)
        +plotWorld(PlotOrientationOn, MagField)
        -PlotMagenticMoments(LinkPos, Moments)
        -DetermineMagneticContribution()
        -calcFieldContribution(Pos, Moment) double double double
        -InitMagField()
        -EvaluateMagneticForces()
        -findClosestGridPoint(Pos) idx idy idz
    }
    World "1" *-- "1" Tentacle
    Tentacle "1" *-- "1..n" Joint
```
