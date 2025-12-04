// ====================================================================
// GMSH Geometry File: Sphere with Inlet, Outlet, Walls
// ====================================================================

R_sphere = 0.5;
L_domain = 10.0;
X_length = 20;
X_min = -5;
mesh_sphere = 0.1;
mesh_far = 0.5;

SetFactory("OpenCASCADE");

// --------------------------------------------------------------------
// Geometry Construction
// --------------------------------------------------------------------

Sphere(1) = {0, 0, 0, R_sphere};
Box(2) = {X_min, -L_domain/2, -L_domain/2, X_length, L_domain, L_domain};

// Fixed: Proper OpenCASCADE BooleanDifference syntax
BooleanDifference(3) = { Volume{2}; Delete; }{ Volume{1}; };

// Synchronize geometry
Synchronize;

// --------------------------------------------------------------------
// Physical Groups - Define Inlet, Outlet, Walls, Sphere
// --------------------------------------------------------------------

// Get all surfaces
all_surfs() = Surface{:};

// Find sphere surface (smallest area)
minArea = 1e9;
sphere_surf = 0;
For i In {0:#all_surfs()-1}
  area = Abs(Area{Surface{all_surfs(i)};});
  If (area < minArea)
    minArea = area;
    sphere_surf = all_surfs(i);
  EndIf
EndFor

// Identify cube faces by location
inlet_surf = 0;
outlet_surf = 0;
wall_surfs() = {};

For i In {0:#all_surfs()-1}
  If (all_surfs(i) != sphere_surf)
    // Get surface center coordinates
    bbox[] = BoundingBox Surface{all_surfs(i)};
    x_center = (bbox[0] + bbox[3]) / 2;
    
    If (Fabs(x_center - X_min) < 0.01)
      inlet_surf = all_surfs(i);
    ElseIf (Fabs(x_center - (X_min + X_length)) < 0.01)
      outlet_surf = all_surfs(i);
    Else
      wall_surfs(#wall_surfs()) = all_surfs(i);
    EndIf
  EndIf
EndFor

// Define physical surfaces
Physical Surface("sphere") = {sphere_surf};
Physical Surface("inlet") = {inlet_surf};
Physical Surface("outlet") = {outlet_surf};
Physical Surface("walls") = {wall_surfs()};
Physical Volume("fluid") = {3};

// --------------------------------------------------------------------
// Mesh Control
// --------------------------------------------------------------------

Field[1] = Distance;
Field[1].SurfacesList = {sphere_surf};
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = mesh_sphere;
Field[2].SizeMax = mesh_far;
Field[2].DistMin = R_sphere;
Field[2].DistMax = 3*R_sphere;
Field[2].Sigmoid = 1;
Background Field = 2;

Mesh.MshFileVersion = 2.2;
Mesh.ElementOrder = 1;
Mesh.Algorithm3D = 1;
Mesh.OptimizeNetgen = 1;
Mesh.Smoothing = 5;
Mesh.Binary = 0;

// ====================================================================

