# Flow Around a Sphere: Numerical Simulations with Adaptive Mesh Refinement

This directory contains calculations for compressible flow around a sphere using the shockFluid solver with adaptive mesh refinement (AMR).

## Running the Verification Problem (Shock Tube)

To execute the shock tube verification case, run:

./Allrun

To modify the spatial resolution level, edit the file `constant/dynamicMeshDict` and adjust the refinement parameters.

To obtain the analytical solution for comparison, execute the Python script located in the `python/` directory.

## Running the Main Problem (Flow Around Sphere)

Set necessary velocity in the file 0/U to obtain required Mach number.

Execute the baseline calculation using:

./Allrun

The initial phase uses the **upwind scheme interpolation** to establish the flow field and shock structure. Once the main features of the flow develop, continue the simulation with the **van Albada limiter** by running:

./Allrun-continue

## Adaptive Mesh Refinement Strategy

It is recommended to manually adjust the density gradient threshold for refinement at each iteration. As the shock develops and becomes steeper, the gradient magnitude increases, allowing for more accurate shock localization and mesh adaptation.

## Special Considerations for Mach 4 Simulations

At Mach 4, particular care must be taken with the mesh near the sphere surface. Second-order schemes can produce spurious oscillations in this region, which may destabilize the calculation.

Based on experimental investigation, the following mesh parameters are recommended:

- **stlSurface refinement level:** 3
- **Boundary layer mesh:** Enable mesh layers around the sphere to resolve the flow near the sphere

These settings prevent numerical oscillations and ensure stable, accurate computation of the complex shock structure at high Mach numbers.

