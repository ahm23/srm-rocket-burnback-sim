# Solid Rocket Motor Burnback Simulator

MATLAB prototype for simulating solid rocket motor (SRM) chamber pressure and propellant burnback. The code models pressure evolution and thrust curve as a balance between propellant gas generation and choked nozzle outflow, with geometry-dependent burning area and combustion chamber volume.

The repository currently contains both simple analytical grain models and a generic numerical solver using computational geometry algorithms and an ODE solver to simulate the burnback of complex geometries.

## Files

| File | Purpose |
| --- | --- |
| `main.m` | Main experiment script. Defines propellant/motor parameters, builds an example finocyl geometry, calls the ODE solver, and plots pressure and grain-state histories. |
| `RocketSystem_Annular.m` | ODE setup for an annular grain whose inner port grows while grain length regresses. |
| `RocketSystem_Cylinder.m` | ODE setup for a cylindrical/end-burning phase with pressure and remaining length as states. |
| `RocketSystem_Generic.m` | ODE setup intended to use arbitrary grain geometry/distance-field data. |
| `genDistanceField.m` | Helper for generating a signed distance field from an outer casing radius and inner-port polygon. |
| `finocyl_pointy.m` | Generates a 2D finocyl-style inner-port polygon with pointed fins. |
| `configure.m` | Configuration stub for distance-field resolution settings. |

## Requirements

- MATLAB
- Aerospace Toolbox, if using `atmosisa` from `main.m`
- An `inpoly2` implementation on the MATLAB path, if using `genDistanceField.m`

## Running

From MATLAB, open the repository root and run:

```matlab
main
```

The intended output is a set of figures showing:

- chamber pressure versus time
- inner grain radius versus time
- outer grain radius versus time
- grain length versus time
- taper angle versus time

## Safety Note

This is an engineering simulation prototype, not a validated motor-design tool. Use independent analysis, conservative margins, and qualified review before applying results to any physical rocket motor.
