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

import os
import pathlib
import argparse
import re
import numpy as np
from math import floor

parser = argparse.ArgumentParser(description='CPW line simulation')
parser.add_argument('--plot', action='store_true',
                     help='Open the AppCSXCAD geometry viewer and show the S-parameter plot window')
args = parser.parse_args()

from CSXCAD  import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import *
from openEMS.ports import CPWPort
# from scipy.constants import epsilon_0
epsilon_0 = 8.8541878128e-12

print(f'epsilon_0 is {epsilon_0}')

import os
os.environ["VTK_SILENCE_DEPRECATION_WARNINGS"] = "1"

from CSXCAD import ContinuousStructure

### Setup the simulation
# Anchor output paths to this script's own location rather than the
# current working directory. openEMS's FDTD.Run() changes the process's
# cwd internally and does not restore it, so any later "./relative" paths
# silently end up somewhere unexpected (e.g. under /tmp) instead of next
# to this script -- which is why plots could go missing in CI.
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
print(f'Script directory (anchor for outputs): {SCRIPT_DIR}')
print(f'Current working directory at startup: {pathlib.Path.cwd()}')

Sim_Path = str(SCRIPT_DIR / 'pcb_large_sim')

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
mesh.AddLine('x', [-air_spacing, PCB_LENGTH + air_spacing])

num_resistors = 3
#2512 resistor pad boundaries
res_pad_width  = 1225   # um
res_pad_height = 3350   # um
res_left_x     = [20037.5, 47037.5, 47037.5] 
res_right_x    = [25962.5, 52962.5, 52962.5]
res_y          = [20000.0, 10000, 30000]

mesh.AddLine('x', [
    res_left_x[0]  - res_pad_width/2,
    res_left_x[0]  + res_pad_width/2,
    res_right_x[0] - res_pad_width/2,
    res_right_x[0] + res_pad_width/2
])
x_resolution = 100
mesh.SmoothMeshLines('x', x_resolution)

mesh.AddLine('x', [
    res_left_x[1]  - res_pad_width/2,
    res_left_x[1]  + res_pad_width/2,
    res_right_x[1] - res_pad_width/2,
    res_right_x[1] + res_pad_width/2
])
x_resolution = 100
mesh.SmoothMeshLines('x', x_resolution)

# y-mesh: fine at CPW signal and gap edges, coarse elsewhere
#car between SMAs or Y=6 to 34
#care less between 0<Y<6 and 34<Y<40
third_mesh = np.array([-2/3, 1/3]) * edge_res
mesh.AddLine('y', [0])
mesh.AddLine('y', [PCB_WIDTH])
mesh.AddLine('y', PCB_WIDTH/2 + trace_width/2 + third_mesh)
mesh.AddLine('y', PCB_WIDTH/2 + trace_width/2 + CPW_gap - third_mesh)
mesh.AddLine('y', PCB_WIDTH/2 - trace_width/2 - third_mesh)
mesh.AddLine('y', PCB_WIDTH/2 - trace_width/2 - CPW_gap + third_mesh)
mesh.AddLine('y', 10000 + trace_width/2 + third_mesh)
mesh.AddLine('y', 10000 + trace_width/2 + CPW_gap - third_mesh)
mesh.AddLine('y', 10000 - trace_width/2 - third_mesh)
mesh.AddLine('y', 10000 - trace_width/2 - CPW_gap + third_mesh)

mesh.AddLine('y', 30000 + trace_width/2 + third_mesh)
mesh.AddLine('y', 30000 + trace_width/2 + CPW_gap - third_mesh)
mesh.AddLine('y', 30000 - trace_width/2 - third_mesh)
mesh.AddLine('y', 30000 - trace_width/2 - CPW_gap + third_mesh)

# mesh.AddLine('y', [
#     res_y - res_pad_height/2,
#     res_y + res_pad_height/2
# ])

mesh.SmoothMeshLines('y', edge_res*1.5, ratio=1.5)
print(f'resolution = {resolution}')
resolution = C0 / (f_max * np.sqrt(substrate_epr)) / unit / 30
mesh.SmoothMeshLines('y', resolution, ratio=1.3)

# z-mesh: fine inside substrate, coarse in air
# mesh.AddLine('z', [0])
# mesh.AddLine('z', [2*air_spacing])
# mesh.SmoothMeshLines('z', resolution)

# Z mesh
# Z mesh
z_bottom = 0
z_substrate_start = air_spacing
z_substrate_end = air_spacing + PCB_THICKNESS
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
portstart = [               0, PCB_WIDTH/2 - trace_width/2, air_spacing+PCB_THICKNESS]
portstop  = [CPW_port_length,  PCB_WIDTH/2 + trace_width/2, air_spacing+PCB_THICKNESS]
print(f'port1_start = {portstart}')
print(f'port1_stop  = {portstop}')
port1 = CPWPort(CSX, 1, cpw_port_metal, portstart, portstop, 'x', 'z', CPW_gap,
                excite=1, priority=999,
                MeasPlaneShift=CPW_port_length, Feed_R=feed_R)

portstart = [ PCB_LENGTH                  , 30000 - trace_width/2, air_spacing+PCB_THICKNESS]
portstop  = [ PCB_LENGTH - CPW_port_length, 30000 + trace_width/2, air_spacing+PCB_THICKNESS]
print(f'port2_start = {portstart}')
print(f'port2_stop  = {portstop}')
port2 = CPWPort(CSX, 2, cpw_port_metal, portstart, portstop, 'x', 'z', CPW_gap,
                priority=999, MeasPlaneShift=CPW_port_length, Feed_R=feed_R)

portstart = [ PCB_LENGTH                  , 10000 - trace_width/2, air_spacing+PCB_THICKNESS]
portstop  = [ PCB_LENGTH - CPW_port_length, 10000 + trace_width/2, air_spacing+PCB_THICKNESS]
print(f'port3_start = {portstart}')
print(f'port3_stop  = {portstop}')
port3 = CPWPort(CSX, 3, cpw_port_metal, portstart, portstop, 'x', 'z', CPW_gap,
                priority=999, MeasPlaneShift=CPW_port_length, Feed_R=feed_R)

ports = [port1, port2, port3]

### CPW centre conductor between the two ports
trace_in = CSX.AddMetal('TRACE_IN')
start = [             CPW_port_length, PCB_WIDTH/2 + trace_width/2, air_spacing+PCB_THICKNESS]
stop  = [19.5*1000, PCB_WIDTH/2 - trace_width/2, air_spacing+PCB_THICKNESS]
print(f'trace_start = {start}')
print(f'trace_stop  = {stop}')
trace_in.AddBox(start, stop, priority=999)


### CPW centre conductor between the two ports
trace_middle = CSX.AddMetal('TRACE_MIDDLE')
currDir = os.getcwd()
trace_middle_poly = trace_middle.AddPolyhedronReader(os.path.join(currDir, 'middle_trace.stl'), priority=9900)
trace_middle_poly.ReadFile()
trace_middle_poly.AddTransform('Scale', [1000, -1000, 1000]) #mm to um, y is inverted in kicad
trace_middle_poly.AddTransform('Translate', [0, 0, air_spacing + PCB_THICKNESS]) # lifting up


trace_out2 = CSX.AddMetal('TRACE_OUT2')
start = [PCB_LENGTH-CPW_port_length, 30000 - trace_width/2, air_spacing+PCB_THICKNESS]
stop  = [                     53000, 30000 + trace_width/2, air_spacing+PCB_THICKNESS]
print(f'trace_start = {start}')
print(f'trace_stop  = {stop}')
trace_out2.AddBox(start, stop, priority=999)

trace_out3 = CSX.AddMetal('TRACE_OUT3')
start = [PCB_LENGTH-CPW_port_length, 10000 - trace_width/2, air_spacing+PCB_THICKNESS]
stop  = [                     53000, 10000 + trace_width/2, air_spacing+PCB_THICKNESS]
print(f'trace_start = {start}')
print(f'trace_stop  = {stop}')
trace_out3.AddBox(start, stop, priority=999)


### 2512 resistor
resistor_R = 16.0  # Ohms
res_z = air_spacing + PCB_THICKNESS

for i in range(num_resistors):
    res_pad1 = CSX.AddMetal(f"RESISTOR{i}_PAD1")
    # Left pad
    res_pad1.AddBox(
        [res_left_x[i] - res_pad_width/2,
        res_y[i] - res_pad_height/2,
        res_z],
        [res_left_x[i] + res_pad_width/2,
        res_y[i] + res_pad_height/2,
        res_z],
        priority=999
    )

    res_pad2 = CSX.AddMetal(f"RESISTOR{i}_PAD2")
    # Right pad
    res_pad2.AddBox(
        [res_right_x[i] - res_pad_width/2,
        res_y[i] - res_pad_height/2,
        res_z],
        [res_right_x[i] + res_pad_width/2,
        res_y[i] + res_pad_height/2,
        res_z],
        priority=999
    )

    # # 16-ohm lumped element across the gap
    resistor = CSX.AddLumpedElement(
        f"RESISTOR_{i}_16R",
        ny=0,
        R=resistor_R
    )
    res_height = 550
    resistor.AddBox(
        [res_left_x[i],
        res_y[i] - res_pad_height/2,
        res_z],
        [res_right_x[i],
        res_y[i] + res_pad_height/2,
        res_z+res_height],
        priority=1000
    )

### CPW ground planes
gnd = CSX.AddMetal('GND')

# Top grounds
start = [          0, 19.49*1000, air_spacing+PCB_THICKNESS]
stop  = [ 19.28*1000,          0, air_spacing+PCB_THICKNESS]
gnd.AddBox(start, stop, priority=999)

start = [          0, 20.51*1000, air_spacing+PCB_THICKNESS]
stop  = [ 19.28*1000,    40*1000, air_spacing+PCB_THICKNESS]
gnd.AddBox(start, stop, priority=999)

start = [ 19.28*1000, 18.12*1000, air_spacing+PCB_THICKNESS]
stop  = [ 26.83*1000,          0, air_spacing+PCB_THICKNESS]
gnd.AddBox(start, stop, priority=999)

start = [ 19.28*1000, 21.88*1000, air_spacing+PCB_THICKNESS]
stop  = [ 26.83*1000, 40.00*1000, air_spacing+PCB_THICKNESS]
gnd.AddBox(start, stop, priority=999)

start = [ 26.83*1000, 19.49*1000, air_spacing+PCB_THICKNESS]
stop  = [ 35.2*1000,           0, air_spacing+PCB_THICKNESS]
gnd.AddBox(start, stop, priority=999)

start = [ 26.83*1000, 20.51*1000, air_spacing+PCB_THICKNESS]
stop  = [ 35.2*1000,  40.00*1000, air_spacing+PCB_THICKNESS]
gnd.AddBox(start, stop, priority=999)

start = [ 53.78*1000, 09.49*1000, air_spacing+PCB_THICKNESS]
stop  = [ PCB_LENGTH,  0.00*1000, air_spacing+PCB_THICKNESS]
gnd.AddBox(start, stop, priority=999)

start = [ 53.78*1000, 30.51*1000, air_spacing+PCB_THICKNESS]
stop  = [ PCB_LENGTH,  40.00*1000, air_spacing+PCB_THICKNESS]
gnd.AddBox(start, stop, priority=999)

start = [ 53.78*1000, 10.51*1000, air_spacing+PCB_THICKNESS]
stop  = [ PCB_LENGTH, 29.49*1000, air_spacing+PCB_THICKNESS]
gnd.AddBox(start, stop, priority=999)

start = [ 46.21*1000, 31.88*1000, air_spacing+PCB_THICKNESS]
stop  = [ 53.78*1000, 40.00*1000, air_spacing+PCB_THICKNESS]
gnd.AddBox(start, stop, priority=999)

start = [ 46.21*1000, 00.00*1000, air_spacing+PCB_THICKNESS]
stop  = [ 53.78*1000, 08.12*1000, air_spacing+PCB_THICKNESS]
gnd.AddBox(start, stop, priority=999)

start = [ 46.21*1000, 11.88*1000, air_spacing+PCB_THICKNESS]
stop  = [ 53.78*1000, 28.13*1000, air_spacing+PCB_THICKNESS]
gnd.AddBox(start, stop, priority=999)

start = [ 44.90*1000, 00.00*1000, air_spacing+PCB_THICKNESS]
stop  = [ 46.21*1000, 09.49*1000, air_spacing+PCB_THICKNESS]
gnd.AddBox(start, stop, priority=999)

start = [ 44.90*1000, 30.51*1000, air_spacing+PCB_THICKNESS]
stop  = [ 46.21*1000, 40.00*1000, air_spacing+PCB_THICKNESS]
gnd.AddBox(start, stop, priority=999)

start = [ 45.42*1000, 10.51*1000, air_spacing+PCB_THICKNESS]
stop  = [ 46.21*1000, 29.49*1000, air_spacing+PCB_THICKNESS]
gnd.AddBox(start, stop, priority=999)

points = [
    [45.42*1000, 45.42*1000,  36.20*1000],
    [10.51*1000, 29.49*1000,  20.00*1000],
]
gnd.AddPolygon(
    points=points,
    norm_dir='z',
    elevation=air_spacing+PCB_THICKNESS,
    priority=999
)

points = [
    [35.20*1000, 35.20*1000,  44.96*1000, 44.96*1000],
    [00.00*1000, 19.49*1000,  09.49*1000, 0*1000],
]
gnd.AddPolygon(
    points=points,
    norm_dir='z',
    elevation=air_spacing+PCB_THICKNESS,
    priority=999
)

points = [
    [35.20*1000, 35.20*1000,  44.96*1000, 44.96*1000],
    [40.00*1000, 20.51*1000,  30.51*1000, 40*1000],
]
gnd.AddPolygon(
    points=points,
    norm_dir='z',
    elevation=air_spacing+PCB_THICKNESS,
    priority=999
)


# Bottom ground
start = [         0,         0, air_spacing]
stop  = [PCB_LENGTH, PCB_WIDTH, air_spacing]
gnd.AddBox(start, stop, priority=999)

### Via / drill-hole geometry — imported from the KiCad drill file
#
# Parses every tool in the drill file (not just the GND stitching vias),
# so the connectors' component holes (T2/T3) come in with their own real
# diameters too, instead of a single hard-coded size.

DRILL_FILE = (SCRIPT_DIR / '../pcb/Resistive Splitter-PTH.drl').resolve() 


def parse_drill_file(drl_path):
    """
    Parse a full KiCad/Excellon drill file (metric, decimal, absolute
    format -- the KiCad default) and return every tool's diameter,
    aperture function, and hole coordinates.

    Returns: dict tool_num -> {
        'diameter': float (mm),
        'function': str or None   (e.g. 'ViaDrill', 'ComponentDrill'),
        'holes':    [(x_mm, y_mm), ...],
    }
    """
    with open(drl_path) as f:
        drl_lines = f.read().splitlines()

    # --- Tool table: map tool number -> {diameter_mm, function} ---
    tools = {}
    current_function = None
    for line in drl_lines:
        line = line.strip()
        if line.startswith('; #@! TA.AperFunction'):
            current_function = line.split(',')[-1]
            continue
        m = re.match(r'^T(\d+)C([\d.]+)$', line)
        if m:
            tools[int(m.group(1))] = {
                'diameter': float(m.group(2)),
                'function': current_function,
                'holes': [],
            }
            continue
        if line == '%':
            break  # end of header / tool table

    # --- Coordinate blocks: a bare "Tn" line selects the active tool ---
    active_tool = None
    coord_re = re.compile(r'^X([-\d.]+)Y([-\d.]+)$')
    for line in drl_lines:
        line = line.strip()
        tm = re.match(r'^T(\d+)$', line)
        if tm:
            active_tool = int(tm.group(1))
            continue
        cm = coord_re.match(line)
        if cm and active_tool is not None and active_tool in tools:
            tools[active_tool]['holes'].append((float(cm.group(1)), float(cm.group(2))))
    return tools


def holes_to_boxes(holes_mm, diameter_mm, z_bottom, z_top, pad_margin_um=0.0):
    """
    Convert a list of (x_mm, y_mm) drill hits -- all sharing one drill
    diameter -- into [x_start, x_end, y_start, y_end, z_start, z_end]
    boxes in the sim's um-scaled coordinates. Each hole becomes a solid
    square post the size of the finished drill diameter (plus an
    optional pad_margin_um, to approximate the pad/annular ring instead
    of just the bare drilled hole).

    Coordinate mapping from Excellon -> sim frame (matches this board):
      x_sim =  x_drill_mm * 1000   (drill X already runs 0..PCB_LENGTH)
      y_sim = -y_drill_mm * 1000   (drill Y is negative-down in the file;
                                     sim Y runs 0..PCB_WIDTH, positive-up)
    Coincident duplicate holes (can happen on a real board) are
    collapsed to a single box.
    """
    r_um = diameter_mm * 1000.0 / 2.0 + pad_margin_um
    boxes = []
    seen = set()
    n_dupes = 0
    for x_mm, y_mm in holes_mm:
        x_um = x_mm * 1000.0
        y_um = -y_mm * 1000.0
        key = (round(x_um, 3), round(y_um, 3))
        if key in seen:
            n_dupes += 1
            continue
        seen.add(key)
        boxes.append([x_um + r_um, x_um - r_um, y_um + r_um, y_um - r_um, z_top, z_bottom])
    if n_dupes:
        print(f'    skipped {n_dupes} duplicate/coincident hole(s)')
    return boxes


print(f'Reading drill data from: {DRILL_FILE}')
drill_tools = parse_drill_file(DRILL_FILE)
for tnum, info in sorted(drill_tools.items()):
    print(f"  T{tnum}: {info['function'] or '?':<15s} dia={info['diameter']:.3f} mm  "
          f"x {len(info['holes'])} holes")

# T1 -- plated stitching vias (ViaDrill) -> GND
VIA_TOOL = next(t for t, info in drill_tools.items() if info['function'] == 'ViaDrill')
via_locations = holes_to_boxes(
    drill_tools[VIA_TOOL]['holes'], drill_tools[VIA_TOOL]['diameter'],
    z_bottom=air_spacing, z_top=air_spacing + PCB_THICKNESS,
)
print(f'{len(via_locations)} via boxes after de-duplication '
      f'(T{VIA_TOOL}, {drill_tools[VIA_TOOL]["diameter"]} mm)')

# T2/T3 -- SMA connector holes (ComponentDrill), each kept at its own real
# diameter. These sit right on top of the CPW ports (0/PCB_LENGTH, y=10/20/30 mm),
# consistent with a connector's grounded pins/legs, so by default they're
# folded into the same GND metal as the stitching vias below. If any of
# these tools are actually non-plated mechanical/mounting holes on your
# board rather than grounded pins, route that tool's boxes to their own
# CSX.AddMetal(...) (or drop them) instead of appending to via_locations.
COMPONENT_TOOLS = sorted(t for t, info in drill_tools.items() if info['function'] == 'ComponentDrill')
connector_locations = []
for t in COMPONENT_TOOLS:
    info = drill_tools[t]
    boxes = holes_to_boxes(info['holes'], info['diameter'],
                            z_bottom=air_spacing, z_top=air_spacing + PCB_THICKNESS)
    print(f'{len(boxes)} connector-pin boxes from T{t} ({info["diameter"]} mm)')
    connector_locations += boxes

via_locations += connector_locations
print(f'{len(via_locations)} total boxes (vias + connector pins) going into GND')

gnd_vias = CSX.AddMetal('GND_VIAS') 
for x_start, x_end, y_start, y_end, z_start, z_end in via_locations:
    start = [x_start, y_start, z_start]
    stop  = [x_end, y_end, z_end]
    # gnd.AddBox(start, stop, priority=999)
    gnd_vias.AddBox(start, stop, priority=999)


# ---------------------------------------------------------------------------
# E-FIELD DUMP FOR PARAVIEW
# ---------------------------------------------------------------------------
# Save the time-domain E-field (one .vtr per timestep) as VTK, on a single
# horizontal plane at mid-substrate-thickness spanning the full PCB
# footprint.
#
# openEMS uses:
#   dump_type=0  : time-domain E-field
#   file_type=0  : VTK output
#   dump_mode=2  : cell-based interpolation
dumpBoxList = {}

dumpboxName = "efield_efield"
dumpBoxList[dumpboxName] = CSX.AddDump(
    dumpboxName, dump_type=0, file_type=0, dump_mode=2, sub_sampling=[2, 2, 2]
)
dumpboxStart = [0,          0,         air_spacing + PCB_THICKNESS/2]
dumpboxStop  = [PCB_LENGTH, PCB_WIDTH, air_spacing + PCB_THICKNESS/2]
dumpBoxList[dumpboxName].AddBox(dumpboxStart, dumpboxStop)





# Write geometry and open AppCSXCAD
simdir = SCRIPT_DIR / "simulation"
xmlname = pathlib.Path("simulation.xml")

if not simdir.exists():
    simdir.mkdir(parents=True)

xmlpath = simdir / xmlname

print(f'Writing CSXCAD geometry to: {xmlpath.resolve()}')
CSX.Write2XML(str(xmlpath))
if args.plot:
    os.system(f'~/opt/openEMS/bin/AppCSXCAD "{xmlpath}"')

run_sim = 1

if run_sim == 1:
    ### Run the simulation
    print(f'Current working directory before FDTD.Run: {pathlib.Path.cwd()}')
    FDTD.Run(Sim_Path, cleanup=True)
    # Note: FDTD.Run() is known to os.chdir() into Sim_Path internally and
    # not restore the original cwd afterward. Printing it here makes that
    # visible instead of silently breaking any "./relative" paths below.
    print(f'Current working directory after FDTD.Run: {pathlib.Path.cwd()}')

    ### Post-processing
    f = np.linspace(1e6, f_max, 1601)
    for port in ports:
        port.CalcPort(Sim_Path, f, ref_impedance=50)

    s11 = ports[0].uf_ref / ports[0].uf_inc
    s21 = ports[1].uf_ref / ports[0].uf_inc
    s31 = ports[2].uf_ref / ports[0].uf_inc

    print('ports[0].uf_ref is')
    print(ports[0].uf_ref)
    print('ports[1].uf_ref is')
    print(ports[1].uf_ref)
    print('ports[2].uf_ref is')
    print(ports[2].uf_ref)

    s11_dB = 20 * np.log10(np.abs(s11))
    s21_dB = 20 * np.log10(np.abs(s21))
    s31_dB = 20 * np.log10(np.abs(s31))

    ### Pass / fail checks
    mask = f > 100e6
    print(f'max(dB(S11)) = {np.max(s11_dB[mask]):.1f} dB')
    print(f'min(dB(S21)) = {np.min(s21_dB[mask]):.1f} dB,  max(dB(S21)) = {np.max(s21_dB[mask]):.2f} dB')
    print(f'min(dB(S21)) = {np.min(s31_dB[mask]):.1f} dB,  max(dB(S31)) = {np.max(s31_dB[mask]):.2f} dB')

    # assert np.max(s11_dB[mask]) < -20, \
    #     f'FAIL: max(dB(S11)) = {np.max(s11_dB[mask]):.1f} dB, expected < -20 dB'
    # assert np.min(s21_dB[mask]) > -0.5, \
    #     f'FAIL: min(dB(S21)) = {np.min(s21_dB[mask]):.1f} dB, expected > -0.5 dB'
    # assert np.max(s21_dB[mask]) < 0.01, \
    #     f'FAIL: max(dB(S21)) = {np.max(s21_dB[mask]):.2f} dB, expected < +0.01 dB (sign error?)'

    # print('PASS')

    if 1:  # set to 1 for debugging plots
        import matplotlib
        if not args.plot:
            matplotlib.use('Agg')  # headless backend, no display needed
        import matplotlib.pyplot as plt

        fig, axis = plt.subplots(3, 1, num='S-Parameters', tight_layout=True, sharex=True)

        # S11
        axis[0].plot(f/1e9, s11_dB, 'k-', linewidth=2)
        axis[0].set_ylabel('$S_{11}$ (dB)')
        axis[0].grid()
        axis[0].set_xmargin(0)

        # S21
        axis[1].plot(f/1e9, s21_dB, 'r--', linewidth=2)
        axis[1].set_ylabel('$S_{21}$ (dB)')
        axis[1].grid()
        axis[1].set_xmargin(0)

        # S31
        axis[2].plot(f/1e9, s31_dB, 'b-.', linewidth=2)
        axis[2].set_ylabel('$S_{31}$ (dB)')
        axis[2].set_xlabel('Frequency (GHz)')
        axis[2].grid()
        axis[2].set_xmargin(0)

        results_dir = SCRIPT_DIR / "results"
        results_dir.mkdir(parents=True, exist_ok=True)
        plot_path = results_dir / "s_parameters.png"
        print(f'Current working directory before saving plot: {pathlib.Path.cwd()}')
        print(f'Saving S-parameter plot to: {plot_path.resolve()}')
        fig.savefig(plot_path, dpi=200)
        print(f'Plot saved: {plot_path.resolve()}  (exists: {plot_path.exists()})')

        if args.plot:
            plt.show()