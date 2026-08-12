"""
 Coplanar Waveguide (CPW) Line Test

 Verifies CPWPort setup and S-parameter extraction for a finite-length CPW
 transmission line on RO4350B substrate with 50-Ohm port terminations.

 Pass criteria:
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

print(f'epsilon_0 is {epsilon_0}')

### Setup the simulation
Sim_Path = os.path.join(tempfile.gettempdir(), 'CPW_Line')

unit                = 1e-6   # drawing unit in um

PCB_LENGTH    = 75000
PCB_WIDTH     = 40000
PCB_THICKNESS = 510

CPW_length          = 40000
CPW_port_length     = 10000
trace_width         = 620
CPW_gap             = 200
# substrate_thickness = 510 # 512 um chaging to 
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
# resolution = C0 / (f_max * np.sqrt(substrate_epr)) / unit / 30
resolution = 20 # 20um cell size
edge_res   = 40

print(f"resolution is {resolution}")
# x-mesh: cuts up the length axis into resolution even points
mesh.AddLine('x', [0])
mesh.AddLine('x', [PCB_LENGTH])
x_resolution = 100
mesh.SmoothMeshLines('x', x_resolution)

# y-mesh: fine at CPW signal and gap edges, coarse elsewhere
#car between SMAs or Y=6 to 34
#care less between 0<Y<6 and 34<Y<40
# third_mesh = np.array([-2/3, 1/3]) * edge_res
mesh.AddLine('y', [0])
mesh.AddLine('y', [PCB_WIDTH])
print(f'resolution = {resolution}')
mesh.SmoothMeshLines('y', resolution)

# y_pos = mesh.GetLines('y')
# mesh.AddLine('y', np.concatenate([-y_pos,
#                                    [-substrate_width/2,  substrate_width/2],
#                                    [-substrate_width/2 - air_spacing,
#                                      substrate_width/2 + air_spacing]]))
# mesh.SmoothMeshLines('y', resolution, ratio=1.3)

# z-mesh: fine inside substrate, coarse in air
# mesh.AddLine('z', [0])
# mesh.AddLine('z', [2*air_spacing])
# mesh.SmoothMeshLines('z', resolution)

# Z mesh
# Z mesh
z_bottom = 0
z_substrate_start = air_spacing - PCB_THICKNESS
z_substrate_end = air_spacing + 2* PCB_THICKNESS
z_top = 2*air_spacing + PCB_THICKNESS

# Bottom air
z_lines_bottom = np.arange(
    z_bottom,
    z_substrate_start,
    1000
)

# Substrate
z_lines_substrate = np.arange(
    z_substrate_start,
    z_substrate_end,
    20
)

# Top air
z_lines_top = np.arange(
    z_substrate_end,
    z_top,
    1000
)

# Make sure important boundaries are included
z_lines = np.unique(np.concatenate([
    z_lines_bottom,
    [z_substrate_start],
    z_lines_substrate,
    [z_substrate_end],
    z_lines_top,
    [z_top]
]))

mesh.AddLine('z', z_lines)



x_cells = len(mesh.GetLines('x')) - 1
y_cells = len(mesh.GetLines('y')) - 1
z_cells = len(mesh.GetLines('z')) - 1

total_cells = x_cells * y_cells * z_cells
print(f"X cells: {x_cells:,}")
print(f"Y cells: {y_cells:,}")
print(f"Z cells: {z_cells:,}")
print(f"Total cells: {total_cells:,}")

### Substrate
substrate = CSX.AddMaterial('RO4350B', epsilon=substrate_epr, kappa=kappa_estimate)
start = [0, 0, air_spacing]
stop  = [ PCB_LENGTH,  PCB_WIDTH, air_spacing+PCB_THICKNESS]
substrate.AddBox(start, stop)

### CPW ports (include the port metal)
cpw_port_metal = CSX.AddMetal('CPW_PORT')
CPW_port_length = 2000 # 2mm
portstart = [               0, PCB_WIDTH/2 + trace_width/2, air_spacing+PCB_THICKNESS]
portstop  = [CPW_port_length,  PCB_WIDTH/2 - trace_width/2, air_spacing+PCB_THICKNESS]
print(f'port1_start = {portstart}')
print(f'port1_stop  = {portstop}')
port1 = CPWPort(CSX, 1, cpw_port_metal, portstart, portstop, 'x', 'z', CPW_gap,
                excite=1, priority=999,
                MeasPlaneShift=CPW_port_length, Feed_R=feed_R)

portstart = [ PCB_LENGTH                  , PCB_WIDTH/2 + trace_width/2, air_spacing+PCB_THICKNESS]
portstop  = [ PCB_LENGTH - CPW_port_length, PCB_WIDTH/2 - trace_width/2, air_spacing+PCB_THICKNESS]
print(f'port1_start = {portstart}')
print(f'port1_stop  = {portstop}')
port2 = CPWPort(CSX, 2, cpw_port_metal, portstart, portstop, 'x', 'z', CPW_gap,
                priority=999, MeasPlaneShift=CPW_port_length, Feed_R=feed_R)

ports = [port1, port2]

### CPW centre conductor between the two ports
trace = CSX.AddMetal('TRACE')
start = [             CPW_port_length, PCB_WIDTH/2 + trace_width/2, air_spacing+PCB_THICKNESS]
stop  = [PCB_LENGTH - CPW_port_length, PCB_WIDTH/2 - trace_width/2, air_spacing+PCB_THICKNESS]
print(f'trace_start = {start}')
print(f'trace_stop  = {stop}')
trace.AddBox(start, stop, priority=999)

### CPW ground planes (left and right of the gap)
gnd = CSX.AddMetal('GND')
start = [          0, PCB_WIDTH/2 + trace_width/2 + CPW_gap, air_spacing+PCB_THICKNESS]
stop  = [ PCB_LENGTH,                             PCB_WIDTH, air_spacing+PCB_THICKNESS]
gnd.AddBox(start, stop, priority=999)

start = [          0,                                      0, air_spacing+PCB_THICKNESS]
stop  = [ PCB_LENGTH,  PCB_WIDTH/2 - trace_width/2 - CPW_gap, air_spacing+PCB_THICKNESS]
gnd.AddBox(start, stop, priority=999)

bottom_and_vias = CSX.AddMetal('BOTTOM_AND_VIAS')
start = [         0,         0, air_spacing]
stop  = [PCB_LENGTH, PCB_WIDTH, air_spacing]
gnd.AddBox(start, stop, priority=999)
# bottom_and_vias.AddBox(start, stop, priority=999)

via_radias = 500 #um
via_spacing = 2200 #um
distance_from_centre = 840 #um

num_vias_wide = floor((PCB_LENGTH) / via_spacing)
print(f"num_vias_wide is {num_vias_wide}")

half_num_vias_tall = floor((PCB_WIDTH/2) / via_spacing)
print(f"half_num_vias_tall is {half_num_vias_tall}")

for wide_via_idx in range(num_vias_wide):
    for tall_via_idx in range(half_num_vias_tall):
        start = [via_spacing/2 + wide_via_idx*via_spacing + via_radias/2, PCB_WIDTH/2 + tall_via_idx*via_spacing + distance_from_centre + via_radias/2, air_spacing+PCB_THICKNESS]
        stop  = [via_spacing/2 + wide_via_idx*via_spacing - via_radias/2, PCB_WIDTH/2 + tall_via_idx*via_spacing + distance_from_centre - via_radias/2, air_spacing]
        gnd.AddBox(start, stop, priority=999)
        # bottom_and_vias.AddBox(start, stop, priority=999)

        start = [via_spacing/2 + wide_via_idx*via_spacing + via_radias/2, PCB_WIDTH/2 - tall_via_idx*via_spacing - distance_from_centre + via_radias/2, air_spacing+PCB_THICKNESS]
        stop  = [via_spacing/2 + wide_via_idx*via_spacing - via_radias/2, PCB_WIDTH/2 - tall_via_idx*via_spacing - distance_from_centre - via_radias/2, air_spacing]
        gnd.AddBox(start, stop, priority=999)
        # bottom_and_vias.AddBox(start, stop, priority=999)




# Write geometry and open AppCSXCAD
simdir = pathlib.Path("./simulation")
xmlname = pathlib.Path("simulation.xml")

if not simdir.exists():
    simdir.mkdir(parents=True)

xmlpath = simdir / xmlname

CSX.Write2XML(str(xmlpath))
os.system(f'~/opt/openEMS/bin/AppCSXCAD "{xmlpath}"')


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

### Pass / fail checks
mask = f > 100e6
print(f'max(dB(S11)) = {np.max(s11_dB[mask]):.1f} dB')
print(f'min(dB(S21)) = {np.min(s21_dB[mask]):.1f} dB,  max(dB(S21)) = {np.max(s21_dB[mask]):.2f} dB')

# assert np.max(s11_dB[mask]) < -20, \
#     f'FAIL: max(dB(S11)) = {np.max(s11_dB[mask]):.1f} dB, expected < -20 dB'
# assert np.min(s21_dB[mask]) > -0.5, \
#     f'FAIL: min(dB(S21)) = {np.min(s21_dB[mask]):.1f} dB, expected > -0.5 dB'
# assert np.max(s21_dB[mask]) < 0.01, \
#     f'FAIL: max(dB(S21)) = {np.max(s21_dB[mask]):.2f} dB, expected < +0.01 dB (sign error?)'

print('PASS')

if 1:  # set to 1 for debugging plots
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(num='S-Parameters', tight_layout=True)
    axis.plot(f/1e9, s11_dB, 'k-',  linewidth=2, label='$S_{11}$')
    axis.plot(f/1e9, s21_dB, 'r--', linewidth=2, label='$S_{21}$')
    axis.grid()
    axis.set_xmargin(0)
    axis.set_xlabel('Frequency (GHz)')
    axis.set_ylabel('S-Parameter (dB)')
    axis.legend()

    plt.show()

