"""Link the short acceptance driver to the existing development standalone objects."""
import argparse
from pathlib import Path
import shlex
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("base", type=Path)
parser.add_argument("output", type=Path)
args = parser.parse_args()
here = Path(__file__).resolve().parent
base = args.base.resolve()
out = args.output.resolve()
out.mkdir(parents=True, exist_ok=True)
target = base / "CMakeFiles/LevelScopeStandalone.dir"
flags = (target / "flags.make").read_text().splitlines()
options = []
for key in ("CXX_DEFINES", "CXX_INCLUDES", "CXX_FLAGSarm64"):
    options.extend(shlex.split(next(line.split(" = ", 1)[1] for line in flags if line.startswith(key + " = "))))
driver = out / "ShortAcceptance.o"
# Access is relaxed only in the driver to observe paint and supply a file selection
# and injected device. All app objects retain the shipping build settings.
subprocess.run(["/usr/bin/c++", *options, "-x", "objective-c++", "-fno-access-control", "-I" + str(here.parents[1]),
                "-c", str(here / "ShortAcceptance.cpp"), "-o", str(driver)], check=True)
link = shlex.split((target / "link.txt").read_text())
i = link.index("x86_64")
del link[i-1:i+1]
i = next(i for i, item in enumerate(link) if item.endswith("/StandaloneApp.cpp.o"))
link[i] = str(driver)
link[link.index("-o") + 1] = str(out / "short-acceptance")
subprocess.run(link, cwd=base, check=True)
