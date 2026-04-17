import numpy as np
import subprocess
import os

# -------------------------
# Executable
# -------------------------
repertoire = ''
executable = './main'

# -------------------------
# Base parameters
# -------------------------
input_parameters = {
    'G': 6.674e-11,

    'mT': 5.972e24,
    'mL': 7.3477e22,
    'mA': 8500,

    'rT': 6378.1e3,
    'rL': 1737.4e3,
    'rA': 2.51,

    'rho0': 0.5,
    'lambda': 7238.2,
    'Cx': 0.3,
    'S': 19.7920337176,

    'tfin': 172800,
    'dt': 10,
    'dt_min': 1e-3,
    'dt_max': 300,
    'epsilon': 1e-6,
    'adaptive': 0,
    'nsampling': 10,

    'Nbody': 3,

    'indice_terre': 0,
    'indice_lune': 1,
    'indice_sonde': 2,
    'indice_atmosphere': 0,

    'xT0': 0,
    'yT0': 0,
    'vxT0': 0,
    'vyT0': 0,

    'xL0': 384748e3,
    'yL0': 0,
    'vxL0': 0,
    'vyL0': 1022,

    'xA0': 314159e3,
    'yA0': 0,
    'vxA0': 0,
    'vyA0': 1200,

    'h': 10000,

    'drag_T': 0,
    'drag_L': 0,
    'drag_A': 1,

    'output': 'trajectory.dat'
}

# -------------------------
# Parameter to scan
# -------------------------
paramstr = 'dt'
variable_array = 2**np.arange(3, 15)

# -------------------------
# Output naming
# -------------------------
outstr = f"epsilon_{input_parameters['epsilon']:.2g}"
outdir = f"Scan_{paramstr}_{outstr}"
os.makedirs(outdir, exist_ok=True)

print("Saving results in:", outdir)

# -------------------------
# Scan loop
# -------------------------
for value in variable_array:
    params = input_parameters.copy()
    params[paramstr] = value

    output_file = f"{outstr}_{paramstr}_{value}.txt"
    output_path = os.path.join(outdir, output_file)

    params['output'] = output_path

    param_string = " ".join(f"{k}={v}" for k, v in params.items())

    cmd = f"{repertoire}{executable} {param_string}"

    print(cmd)
    subprocess.run(cmd, shell=True, check=True)
    print("Done.")