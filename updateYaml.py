import yaml
import sys
from pathlib import Path

help = f"""
Use {sys.argv[0]} to update/add data to a yaml file
WARNING: this will override your file
WARNING: the "yes" will become "true" and so for the "no"
Usage: 
    python {sys.argv[0]} <yamlfile> <keys>...<any number of keys> <key_to_update> <data>
or 
    python {sys.argv[0]} <yamlfile> <key_to_update> <data>

Examples:
    - `python {sys.argv[0]} "tests/gromacs/info.yml" install_plumed master working`
    will update or add the value "working" to the key "master" under the key "install_plumed" in the 'tests/gromacs/info.yml' file
    It will append 
```
install_plumed:
  master: working
```
    to the file

    - `python {sys.argv[0]} "tests/gromacs/info.yml" install plumed master working`
    will update or add the value "working" to the key "master" under the key  "plumed" unter the key "install" in the 'tests/gromacs/info.yml' file
    It will append 
```
install:
  plumed:
    master: working
```
    to the file       
"""

if len(sys.argv) < 3:
    print("Not enough arguments")
    print(help)
    exit(1)

file = Path(sys.argv[1])
path = sys.argv[2:-2]
key_to_update = sys.argv[-2]
data = sys.argv[-1]

info = {}
if file.is_file():
    with open(file, "r") as f:
        info = yaml.load(f, Loader=yaml.SafeLoader)

# print(yaml.dump(info, sort_keys=False))
t = info
for addr in path:
    if addr in t:
        t = t[addr]
    else:
        t[addr] = {}
        t = t[addr]

t[key_to_update] = data

# print("after")
# print(yaml.dump(info, sort_keys=False,Dumper=yaml.SafeDumper))
# todo backup

with open(file, "w") as f:
    f.write(yaml.dump(info, sort_keys=False))
