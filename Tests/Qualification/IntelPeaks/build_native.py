"""Same native driver and product objects, with a private serial service oracle."""
import argparse
from pathlib import Path
import shlex
import subprocess
import sys

p = argparse.ArgumentParser()
p.add_argument('build', type=Path)
p.add_argument('output', type=Path)
a = p.parse_args()
root = Path(__file__).resolve().parents[3]
base, out = a.build.resolve(), a.output.resolve()
subprocess.run([sys.executable, str(root / 'Tests/Qualification/M18/build_driver.py'),
    str(base / 'compile'), str(Path(__file__).with_name('Native.cpp')),
    str(out / 'app-mode'), '--private'], check=True)
link = shlex.split((base / 'compile/CMakeFiles/LevelScopeStandalone.dir/link.txt').read_text())
i = next(i for i, v in enumerate(link) if v.endswith('/StandaloneApp.cpp.o'))
link[i:i+1] = [str(out / 'app-mode.o'), str(base / 'drivers/peaks/serial-measurement.o')]
link[link.index('-o') + 1] = str(out / 'serial')
subprocess.run(link, cwd=base / 'compile', check=True)
