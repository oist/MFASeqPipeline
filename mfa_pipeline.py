#!/usr/bin/env python3
import yaml
import argparse
import atexit
import subprocess
import sys
import os
import shutil
from collections import ChainMap

# Exceptions
class ArgumentError(Exception):
    pass
    
class ConfigurationError(Exception):
    pass
    
class ExecutionError(Exception):
    pass

# Argument handling
parser = argparse.ArgumentParser(
    prog="MFASeqPipeline", 
    description="Process MFA-Seq data")

parser.add_argument("config", help="Config file in YAML format")
parser.add_argument("sample", help="Sample name")

# Parse command-line arguments
ARGS, overrides = parser.parse_known_args()
with open(ARGS.config, "r") as cfg:
    CFG = yaml.safe_load(cfg)

CFG["sample"] = ARGS.sample

# Extra command-line arguments override config
while len(overrides) > 0:
    arg = overrides.pop()
    if not arg.startswith("--"):
      raise ArgumentError(f"invalid argument {arg}")
    key = arg[2:-1]
    if len(overrides) == 0:
      raise ArgumentError(f"no value provided for {arg}")
    value = overrides.pop()
    CFG[key] = value


# Function to expand configuration settings
def expand(value, cfg=CFG):
    if isinstance(value, str):
        vp = value.format_map(cfg)
        i = 0
        while value != vp:
            if i > 100:
              raise ConfigurationError("configuration key {key} expands indefinitely")
            value = vp
            vp = value.format_map(cfg)
            i += 1
    return value


# Keep track for each file whether the newest version is in tmp
TEMPFILES = {}
atexit.register(lambda: [ os.remove(f) for f in TEMPFILES])

# Execute steps
for step_name in STEPS:
    title_name = step_name.upper()
    print(f"*** {title_name}")

    # Within each step, we search for configuration keys first within the
    # step's 'cfg' entry and then at the configuration root. Configuration
    # updates are emphemeral, i.e. lost once the step is completed. This
    # is achieved by putting an empty dictionary which received all updates
    # as the first entry of the dictionary chain we use for lookups and updates
    step = CFG[f"{step_name}_step"]
    if not isinstance(step, dict):
        raise ConfigurationError(f"step {step_name} is not properly defined")
    cfg = ChainMap({}, step["override"], CFG)

    # Get upper-case step name, input and output files, and step script
    input_names = step["input"]
    output_names = step["output"]
    script = step["script"]
    
    # Get expanded input and output file names
    # NOTE: Expansion takes place in the *global* context, since the names
    # are supposed to agree between different steps
    input_files = [expand(CFG[f]) for f in input_names]
    output_files = [expand(CFG[f]) for f in output_names]
    
    print("    input:")
    for f in input_files:
        print(f"         - {f}")
    print("    output:")
    for f in output_files:
        print(f"         - {f}")
    
    # Check if the step needs to be executed
    if all([os.path.exists(f) for f in output_files]):
        print(f"    all outputs already exists, skipping step...")
        continue
    
    # Copy input files to tmp and update config accordingly
    for input_name in input_names:
        # Copy file to tmp dir unless that's where the current version is
        # NOTE: Expansion happens in the *global* context here. For the input file the
        # reason is again that the expansion must agree between steps. For the temp file
        # we could in principle allow step-specific overrides, but then the name of the
        # temp file would depend on the step which created it, which seems silly.
        tmp_dir = expand(CFG["tmp_dir"])
        os.makedirs(tmp_dir, exist_ok=True)
        if input_name not in TEMPFILES:
            # Create a temp file for the input
            input_file = expand(CFG[input_name])
            temp_file = os.path.join(tmp_dir, os.path.basename(input_file))
            # Copy input file to temp file
            print(f"    copying {input_file}")
            print(f"    to {temp_file} ...")
            shutil.copyfile(input_file, temp_file)
            TEMPFILES[input_name] = temp_file
        # Update path in cfg
        cfg[input_name] = TEMPFILES[input_name]

    # Generate temp file names for output files
    for output_name in output_names:
        # The output file shouldn't already exists
        if output_name in TEMPFILES:
            temp_file = TEMPFILES[output_name]
            raise ExecutionError(f"step output {temp_file} already exists")
        # Create a temp file for the output
        temp_file = os.path.join(tmp_dir, os.path.basename(expand(CFG[output_name])))
        TEMPFILES[output_name] = temp_file
        # Update path in cfg
        cfg[output_name] = TEMPFILES[output_name]
        
    # Execute step
    try:
        # Create script file
        # NOTE: We use the *local* context in the expansions here, i.e. allow step-specific overrides
        cfg["script_file"] = os.path.join(expand(cfg["tmp_dir"], cfg), f"{step_name}.sh")
        with open(cfg["script_file"], "w+") as f:
            f.write("#!{shell}\n".format_map(cfg))
            f.write(expand(script, cfg))
        os.chmod(cfg["script_file"], 0o755)
        
        # Run script
        print(f"    running step ...")
        execute_cmd = expand(cfg["execute_cmd"], cfg)
        subprocess.run(execute_cmd, shell=True, executable=cfg["shell"], check=True)
    finally:
        os.remove(cfg["script_file"])
    
    # Copy output files from tmp
    for output_name in output_names:
        temp_file = TEMPFILES[output_name]
        output_file = expand(CFG[output_name])
        if not os.path.exists(temp_file):
            raise ExecutionError(f"step failed to produce {temp_file}")
        print(f"    copying {temp_file}")
        print(f"    to {output_file} ...")
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        shutil.copy(temp_file, output_file)

    print(f"    done.")
