"""
 Coplanar Waveguide (CPW) Line Test — CPW Width Sweep

 Verifies CPWPort setup and S-parameter extraction for a finite-length CPW
 transmission line on RO4350B substrate with 50-Ohm port terminations.

 Ground-via spacing is held constant at 2000 um and CPW_gap is held
 constant at 200 um. This version sweeps the CPW centre-conductor width
 from 580 um to 660 um in 10 um steps, running a full FDTD simulation for
 each width value and saving the S11 and S21 curves as separate PNG
 images per width value:

     S11_CPW_WIDTH=<width>um.png
     S21_CPW_WIDTH=<width>um.png

 into ./cpw_width_sweep/

 NOTE: this runs 9 complete openEMS FDTD simulations back-to-back (one
 per width value), so total runtime will be roughly 9x a single run.

 Pass criteria (per run, printed but not asserted):
   max(dB(S11)) < -20 dB   (low reflection)
   min(dB(S21)) > -0.5 dB  (near-lossless — substrate is modelled lossless)
   max(dB(S21)) < +0.01 dB (sign-error guard)

 Tested with
  - python 3.14
  - openEMS v0.0.36+

 (c) 2026 Thorsten Liebig <thorsten.liebig@gmx.de>

"""

import os, tempfile
import pathlib
import numpy as np
from math import floor

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *
from openEMS.ports import CPWPort
# from scipy.constants import epsilon_0
epsilon_0 = 8.8541878128e-12

import matplotlib
matplotlib.use('Agg')  # non-interactive backend so figures can be saved in a batch loop
import matplotlib.pyplot as plt


def run_simulation(cpw_width):
    """Build geometry, run FDTD, and return (f, s11_dB, s21_dB) for a given CPW_width (um).
    via_spacing is held fixed at 2000 um and CPW_gap fixed at 200 um for this sweep."""

    print(f'epsilon_0 is {epsilon_0}')

    ### Setup the simulation
    Sim_Path = os.path.join(tempfile.gettempdir(), f'CPW_Line_width_{int(round(cpw_width))}um')

    unit                = 1e-6   # drawing unit in um
    CPW_length          = 40000
    CPW_port_length     = 10000
    CPW_width           = cpw_width  # swept parameter (um)
    CPW_gap             = 200      # held constant for this sweep (um)
    via_spacing         = 2000     # held constant for this sweep (um)
    substrate_thickness = 510 # 512 um chaging to
    substrate_width     = 10000
    substrate_epr       = 3.48
    f_max               = 10e9
    substrate_tand      = 0.0031 # 2.5GHz 23degC
    kappa_estimate      = 2 * np.pi * f_max/2 * epsilon_0 * substrate_epr * substrate_tand
    air_spacing         = 7000

    feed_R = 50  # lumped port termination resistance

    ### Setup FDTD parameters & excitation
    FDTD = openEMS(EndCriteria=1e-4)
    FDTD.SetGaussExcite(f_max/2, f_max/2)
    FDTD.SetBoundaryCond(['PMC']*6)

    ### Setup CSXCAD geometry & mesh
    CSX = ContinuousStructure()
    FDTD.SetCSX(CSX)
    mesh = CSX.GetGrid()
    mesh.SetDeltaUnit(unit)

    print(f'C0 is {C0}')
    resolution = C0 / (f_max * np.sqrt(substrate_epr)) / unit / 30
    edge_res   = 40

    # x-mesh: fine at port transitions, coarse in middle and beyond ports
    mesh.AddLine('x', [0])
    mesh.AddLine('x', [-CPW_length/2, -CPW_length/2 + CPW_port_length,
                        CPW_length/2 - CPW_port_length, CPW_length/2])
    mesh.AddLine('x', [-CPW_length/2 - air_spacing, CPW_length/2 + air_spacing])
    mesh.SmoothMeshLines('x', resolution, ratio=1.5)

    # y-mesh: fine at CPW signal and gap edges, coarse elsewhere
    third_mesh = np.array([-2/3, 1/3]) * edge_res
    mesh.AddLine('y', [0])
    mesh.AddLine('y',  CPW_width/2 + third_mesh)
    mesh.AddLine('y',  CPW_width/2 + CPW_gap - third_mesh)
    mesh.SmoothMeshLines('y', edge_res*1.5, ratio=1.5)
    y_pos = mesh.GetLines('y')
    mesh.AddLine('y', np.concatenate([-y_pos,
                                       [-substrate_width/2,  substrate_width/2],
                                       [-substrate_width/2 - air_spacing,
                                         substrate_width/2 + air_spacing]]))
    mesh.SmoothMeshLines('y', resolution, ratio=1.3)

    # z-mesh: fine inside substrate, coarse in air
    mesh.AddLine('z', np.linspace(0, substrate_thickness, 5))
    mesh.AddLine('z', [-air_spacing, substrate_thickness + air_spacing])
    mesh.SmoothMeshLines('z', resolution)

    ### Substrate
    substrate = CSX.AddMaterial('RO4350B', epsilon=substrate_epr, kappa=kappa_estimate)
    start = [-CPW_length/2, -substrate_width/2, 0]
    stop  = [ CPW_length/2,  substrate_width/2, substrate_thickness]
    substrate.AddBox(start, stop)

    ### CPW ports (include the port metal)
    cpw_port_metal = CSX.AddMetal('CPW_PORT')

    portstart = [-CPW_length/2,                  -CPW_width/2, substrate_thickness]
    portstop  = [-CPW_length/2 + CPW_port_length,  CPW_width/2, substrate_thickness]
    port1 = CPWPort(CSX, 1, cpw_port_metal, portstart, portstop, 'x', 'z', CPW_gap,
                    excite=1, priority=999,
                    MeasPlaneShift=CPW_port_length, Feed_R=feed_R)

    portstart = [ CPW_length/2,                  -CPW_width/2, substrate_thickness]
    portstop  = [ CPW_length/2 - CPW_port_length,  CPW_width/2, substrate_thickness]
    port2 = CPWPort(CSX, 2, cpw_port_metal, portstart, portstop, 'x', 'z', CPW_gap,
                    priority=999, MeasPlaneShift=CPW_port_length, Feed_R=feed_R)

    ports = [port1, port2]

    ### CPW centre conductor between the two ports
    cpw = CSX.AddMetal('CPW')
    start = [-CPW_length/2 + CPW_port_length, -CPW_width/2, substrate_thickness]
    stop  = [ CPW_length/2 - CPW_port_length,  CPW_width/2, substrate_thickness]
    cpw.AddBox(start, stop, priority=999)

    ### CPW ground planes (left and right of the gap)
    gnd = CSX.AddMetal('GND')
    start = [-CPW_length/2, -CPW_width/2 - CPW_gap, substrate_thickness]
    stop  = [ CPW_length/2, -substrate_width/2,      substrate_thickness]
    gnd.AddBox(start, stop, priority=999)

    start = [-CPW_length/2,  CPW_width/2 + CPW_gap, substrate_thickness]
    stop  = [ CPW_length/2,  substrate_width/2,      substrate_thickness]
    gnd.AddBox(start, stop, priority=999)

    bottom_and_vias = CSX.AddMetal('BOTTOM_AND_VIAS')
    start = [-CPW_length/2, -substrate_width/2, 0]
    stop = [ CPW_length/2, substrate_width/2, 0]
    # gnd.AddBox(start, stop, priority=999)
    bottom_and_vias.AddBox(start, stop, priority=999)

    via_radias = 500 #um
    # via_spacing held constant at 2000um for this sweep (see top of function)
    distance_from_centre = 840 #um

    half_num_vias_wide = floor(((CPW_length/2)+2*via_radias) / via_spacing)
    print(f"half_num_vias_wide is {half_num_vias_wide}")
    print(f"num_vias_wide is {2*half_num_vias_wide}")

    half_num_vias_tall = floor(((substrate_width/2)+2*via_radias - distance_from_centre) / via_spacing)
    print(f"half_num_vias_tall is {2*half_num_vias_tall}")
    print(f"num_vias_tall is {2*half_num_vias_tall}")

    print(f"num_vias_wide is {half_num_vias_wide}")
    for wide_via_idx in range(half_num_vias_wide):
        for tall_via_idx in range(half_num_vias_tall):
            start = [wide_via_idx*via_spacing + via_radias/2, tall_via_idx*via_spacing + distance_from_centre + via_radias/2, substrate_thickness]
            stop  = [wide_via_idx*via_spacing -via_radias/2, tall_via_idx*via_spacing + distance_from_centre - via_radias/2, 0]
            gnd.AddBox(start, stop, priority=999)
            # bottom_and_vias.AddBox(start, stop, priority=999)

            start = [-wide_via_idx*via_spacing + via_radias/2, tall_via_idx*via_spacing + distance_from_centre + via_radias/2, substrate_thickness]
            stop  = [-wide_via_idx*via_spacing -via_radias/2, tall_via_idx*via_spacing + distance_from_centre - via_radias/2, 0]
            gnd.AddBox(start, stop, priority=999)
            # bottom_and_vias.AddBox(start, stop, priority=999)

            start = [wide_via_idx*via_spacing + via_radias/2, -tall_via_idx*via_spacing - distance_from_centre + via_radias/2, substrate_thickness]
            stop  = [wide_via_idx*via_spacing - via_radias/2, -tall_via_idx*via_spacing - distance_from_centre - via_radias/2, 0]
            gnd.AddBox(start, stop, priority=999)
            # bottom_and_vias.AddBox(start, stop, priority=999)

            start = [-wide_via_idx*via_spacing + via_radias/2, -tall_via_idx*via_spacing -distance_from_centre + via_radias/2, substrate_thickness]
            stop  = [-wide_via_idx*via_spacing - via_radias/2, -tall_via_idx*via_spacing -distance_from_centre - via_radias/2, 0]
            gnd.AddBox(start, stop, priority=999)
            # bottom_and_vias.AddBox(start, stop, priority=999)

    # Write geometry and open AppCSXCAD
    # (absolute path for the same reason out_dir is resolved up front in __main__)
    simdir = pathlib.Path("./simulation").resolve()
    xmlname = pathlib.Path("simulation.xml")

    if not simdir.exists():
        simdir.mkdir(parents=True)

    xmlpath = simdir / xmlname

    # CSX.Write2XML(str(xmlpath))
    # os.system(f'~/opt/openEMS/bin/AppCSXCAD "{xmlpath}"')

    ### Run the simulation
    FDTD.Run(Sim_Path, cleanup=True)

    ### Post-processing
    f = np.linspace(1e6, f_max, 1601)
    for port in ports:
        port.CalcPort(Sim_Path, f, ref_impedance=50)

    s11 = ports[0].uf_ref / ports[0].uf_inc
    s21 = ports[1].uf_ref / ports[0].uf_inc

    s11_dB = 20 * np.log10(np.abs(s11))
    s21_dB = 20 * np.log10(np.abs(s21))

    return f, s11_dB, s21_dB


def save_plot(f, data_dB, label, cpw_width, out_dir):
    """Save a single S-parameter curve to <label>_CPW_WIDTH=<width>um.png"""
    fig, axis = plt.subplots(num=f'{label}_width_{int(round(cpw_width))}um', tight_layout=True)
    sub = label[1:]  # e.g. '11' or '21'
    axis.plot(f/1e9, data_dB, 'k-', linewidth=2, label=f'${label[0]}_{{{sub}}}$')
    axis.grid()
    axis.set_xmargin(0)
    axis.set_xlabel('Frequency (GHz)')
    axis.set_ylabel(f'{label} (dB)')
    axis.set_title(f'CPW_width = {int(round(cpw_width))} um  (CPW_gap = 200 um, via_spacing = 2000 um)')
    axis.legend()

    fname = out_dir / f'{label}_CPW_WIDTH={int(round(cpw_width))}um.png'
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f'  saved {fname}')


if __name__ == '__main__':
    # Resolve to an absolute path up front — FDTD.Run() changes the process's
    # working directory internally (to launch the openEMS binary) and does not
    # restore it, so a relative path here would silently point somewhere else
    # by the time save_plot() is called.
    out_dir = pathlib.Path('./cpw_width_sweep').resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Sweep CPW_width from 580um to 660um in 10um steps (inclusive of both ends -> 9 runs)
    cpw_widths_um = np.arange(580, 900 + 1, 10)

    for cpw_width in cpw_widths_um:
        print(f'\n=== Running simulation for CPW_width = {cpw_width} um (CPW_gap = 200 um, via_spacing = 2000 um) ===')
        f, s11_dB, s21_dB = run_simulation(float(cpw_width))

        ### Pass / fail checks (printed only, not asserted)
        mask = f > 100e6
        print(f'  max(dB(S11)) = {np.max(s11_dB[mask]):.1f} dB')
        print(f'  min(dB(S21)) = {np.min(s21_dB[mask]):.1f} dB,  max(dB(S21)) = {np.max(s21_dB[mask]):.2f} dB')

        save_plot(f, s11_dB, 'S11', cpw_width, out_dir)
        save_plot(f, s21_dB, 'S21', cpw_width, out_dir)

    print('\nAll simulations complete. PNGs saved in ./cpw_width_sweep/')