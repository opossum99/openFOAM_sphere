
import numpy as np


def write_vtk_openfoam_style(filename, x, p, rho, u, gamma):
    """
    Write solution to VTK format compatible with OpenFOAM visualization.

    Uses CELLS format like OpenFOAM's rhoCentralFoam output.
    Creates 1D line elements (VTK_LINE type 3) for proper visualization.

    Parameters:
    -----------
    filename : str
        Output filename (without .vtk extension)
    x : array
        Position array [m]
    p : array
        Pressure [Pa]
    rho : array
        Density [kg/mÂ³]
    u : array
        Velocity [m/s]
    gamma : float
        Heat capacity ratio
    """
    npts = len(x)
    ncells = npts - 1  # number of 1D line cells

    # Calculate derived quantities
    e = p / (rho * (gamma - 1))      # internal energy
    a = np.sqrt(gamma * p / rho)     # sound speed
    Ma = np.abs(u) / a               # Mach number
    E = e + 0.5 * u**2               # total energy

    with open(f"{filename}.vtk", "w") as f:
        # VTK Header (exactly like OpenFOAM)
        f.write("# vtk DataFile Version 2.0\n")
        f.write(f"Sod Shock Tube Solution\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")

        # ========== POINTS SECTION ==========
        f.write(f"POINTS {npts} float\n")
        for xi in x:
            f.write(f"{xi:.8e} 0 0\n")

        # ========== CELLS SECTION (OpenFOAM style) ==========
        # Each cell is a line (VTK type 3) with 2 points
        # Format: CELLS <num_cells> <total_connectivity_size>
        # Connectivity size = sum of (points_per_cell + 1) for all cells
        #                   = (2+1) * ncells = 3 * ncells
        connectivity_size = 3 * ncells

        f.write(f"\nCELLS {ncells} {connectivity_size}\n")
        for i in range(ncells):
            # Format: <num_points> <point_id_1> <point_id_2>
            f.write(f"2 {i} {i+1}\n")

        # ========== CELL TYPES SECTION (OpenFOAM style) ==========
        # Type 3 = VTK_LINE (1D cell)
        f.write(f"\nCELL_TYPES {ncells}\n")
        for i in range(ncells):
            f.write("3\n")

        # ========== POINT DATA SECTION ==========
        f.write(f"\nPOINT_DATA {npts}\n")

        # Pressure
        f.write("SCALARS pressure float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for val in p:
            f.write(f"{val:.8e}\n")

        # Density
        f.write("\nSCALARS density float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for val in rho:
            f.write(f"{val:.8e}\n")

        # Velocity
        f.write("\nSCALARS velocity float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for val in u:
            f.write(f"{val:.8e}\n")

        # Internal energy
        f.write("\nSCALARS internal_energy float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for val in e:
            f.write(f"{val:.8e}\n")

        # Sound speed
        f.write("\nSCALARS sound_speed float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for val in a:
            f.write(f"{val:.8e}\n")

        # Mach number
        f.write("\nSCALARS mach_number float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for val in Ma:
            f.write(f"{val:.8e}\n")

        # Total energy
        f.write("\nSCALARS total_energy float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for val in E:
            f.write(f"{val:.8e}\n")

    print(f"âœ“ OpenFOAM-style VTK file saved: {filename}.vtk")
    print(f"  Format: UNSTRUCTURED_GRID with CELLS")
    print(f"  Points: {npts}")
    print(f"  Cells: {ncells} (VTK_LINE type)")
    print(f"  Variables: 7 (pressure, density, velocity, internal_energy, sound_speed, mach_number, total_energy)")
    return True


def write_vtk_openfoam_cell_data(filename, x, p, rho, u, gamma):
    """
    Write solution with CELL_DATA instead of POINT_DATA (alternative OpenFOAM style).

    Data is stored at cell centers (between grid points).

    Parameters:
    -----------
    filename : str
        Output filename (without .vtk extension)
    x : array
        Position array [m]
    p : array
        Pressure [Pa]
    rho : array
        Density [kg/mÂ³]
    u : array
        Velocity [m/s]
    gamma : float
        Heat capacity ratio
    """
    npts = len(x)
    ncells = npts - 1

    # Calculate derived quantities at cell centers
    e = p / (rho * (gamma - 1))
    a = np.sqrt(gamma * p / rho)
    Ma = np.abs(u) / a
    E = e + 0.5 * u**2

    # Average to cell centers
    p_cell = 0.5 * (p[:-1] + p[1:])
    rho_cell = 0.5 * (rho[:-1] + rho[1:])
    u_cell = 0.5 * (u[:-1] + u[1:])
    e_cell = 0.5 * (e[:-1] + e[1:])
    a_cell = 0.5 * (a[:-1] + a[1:])
    Ma_cell = 0.5 * (Ma[:-1] + Ma[1:])
    E_cell = 0.5 * (E[:-1] + E[1:])

    with open(f"{filename}.vtk", "w") as f:
        # VTK Header
        f.write("# vtk DataFile Version 2.0\n")
        f.write(f"Sod Shock Tube Solution (Cell Data)\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")

        # Points
        f.write(f"POINTS {npts} float\n")
        for xi in x:
            f.write(f"{xi:.8e} 0 0\n")

        # Cells
        connectivity_size = 3 * ncells
        f.write(f"\nCELLS {ncells} {connectivity_size}\n")
        for i in range(ncells):
            f.write(f"2 {i} {i+1}\n")

        # Cell types
        f.write(f"\nCELL_TYPES {ncells}\n")
        for i in range(ncells):
            f.write("3\n")

        # ========== CELL DATA SECTION ==========
        f.write(f"\nCELL_DATA {ncells}\n")

        # Pressure (at cell centers)
        f.write("SCALARS pressure float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for val in p_cell:
            f.write(f"{val:.8e}\n")

        # Density
        f.write("\nSCALARS density float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for val in rho_cell:
            f.write(f"{val:.8e}\n")

        # Velocity
        f.write("\nSCALARS velocity float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for val in u_cell:
            f.write(f"{val:.8e}\n")

        # Internal energy
        f.write("\nSCALARS internal_energy float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for val in e_cell:
            f.write(f"{val:.8e}\n")

        # Sound speed
        f.write("\nSCALARS sound_speed float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for val in a_cell:
            f.write(f"{val:.8e}\n")

        # Mach number
        f.write("\nSCALARS mach_number float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for val in Ma_cell:
            f.write(f"{val:.8e}\n")

        # Total energy
        f.write("\nSCALARS total_energy float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for val in E_cell:
            f.write(f"{val:.8e}\n")

    print(f"âœ“ OpenFOAM-style VTK file saved (CELL_DATA): {filename}.vtk")
    print(f"  Format: UNSTRUCTURED_GRID with CELLS")
    print(f"  Points: {npts}")
    print(f"  Cells: {ncells} (VTK_LINE type)")
    print(f"  Data: CELL_DATA (at cell centers)")
    return True


if __name__ == '__main__':
    import sodshock
    import matplotlib.pyplot as plt

    gamma = 1.4
    dustFrac = 0.0
    npts = 1000
    t = 0.2
    left_state = (1, 1, 0)
    right_state = (0.1, 0.125, 0.)

    print("="*70)
    print("SOD SHOCK TUBE - OpenFOAM-Style VTK Export")
    print("="*70)

    # Solve Riemann problem
    positions, regions, values = sodshock.solve(
        left_state=left_state,
        right_state=right_state,
        geometry=(-1, 1., 0.0),
        t=t,
        gamma=gamma,
        npts=npts,
        dustFrac=dustFrac
    )

    print(f"\nRiemann Problem Parameters:")
    print(f"  gamma = {gamma}")
    print(f"  t = {t} s")
    print(f"  grid points = {npts}")
    print(f"  left state (p, rho, u) = {left_state}")
    print(f"  right state (p, rho, u) = {right_state}")

    # Print solution info
    print(f"\nWave Positions:")
    for desc, vals in positions.items():
        print(f"  {desc:20} : {vals}")

    print(f"\nIntermediate State:")
    for region, vals in sorted(regions.items()):
        print(f"  {region:20} : p={vals[2]}, rho={vals[1]}, u={vals[0]}")

    # Extract arrays
    x = values['x']
    p = values['p']
    rho = values['rho']
    u = values['u']

    print("\n" + "="*70)
    print("Exporting to VTK format (OpenFOAM style)...")
    print("="*70)

    # Export both variants
    write_vtk_openfoam_style("sod_shock_tube_point_data", x, p, rho, u, gamma)
    write_vtk_openfoam_cell_data("sod_shock_tube_cell_data", x, p, rho, u, gamma)

    print("\n" + "="*70)
    print("Files ready for comparison with OpenFOAM:")
    print("\n1. sod_shock_tube_point_data.vtk")
    print("   â””â”€ Data at grid points (POINT_DATA)")
    print("   â””â”€ Use this to compare with OpenFOAM")

    print("\n2. sod_shock_tube_cell_data.vtk")
    print("   â””â”€ Data at cell centers (CELL_DATA)")
    print("   â””â”€ Alternative format")

    print("\nOpenFOAM Comparison in ParaView:")
    print("  1. File â†’ Open â†’ your OpenFOAM result")
    print("  2. File â†’ Open â†’ sod_shock_tube_point_data.vtk")
    print("  3. View both simultaneously (Window â†’ Tile)")
    print("  4. Color by same variable for comparison")
    print("  5. Use color scale range for alignment")
    print("="*70)

    # Matplotlib plot for verification
    print("\nGenerating matplotlib comparison...")
    f, axarr = plt.subplots(3, sharex=True, figsize=(13, 9))

    axarr[0].plot(x, p, linewidth=2, color='b', label='Sodshock (analytical)')
    axarr[0].set_ylabel('Pressure [Pa]', fontweight='bold', fontsize=11)
    axarr[0].grid(True, alpha=0.3, linestyle='--')
    axarr[0].legend(fontsize=10)

    axarr[1].plot(x, rho, linewidth=2, color='r', label='Sodshock (analytical)')
    axarr[1].set_ylabel('Density [kg/mÂ³]', fontweight='bold', fontsize=11)
    axarr[1].grid(True, alpha=0.3, linestyle='--')
    axarr[1].legend(fontsize=10)

    axarr[2].plot(x, u, linewidth=2, color='g', label='Sodshock (analytical)')
    axarr[2].set_ylabel('Velocity [m/s]', fontweight='bold', fontsize=11)
    axarr[2].set_xlabel('Position [m]', fontweight='bold', fontsize=11)
    axarr[2].grid(True, alpha=0.3, linestyle='--')
    axarr[2].legend(fontsize=10)

    f.suptitle(f'Sod Shock Tube Analytical Solution (sodshock)\nt={t}s, Î³={gamma}',
               fontweight='bold', fontsize=13)
    plt.tight_layout()

    print("Plot displayed. Close to continue.\n")
    plt.show()

    print("="*70)
    print("Next step: Compare VTK files with OpenFOAM results in ParaView")
    print("="*70)