#!/usr/bin/env python3
import yaml
import argparse
import atexit
import subprocess
import sys
import os
import io
import shutil
import tempfile
import concurrent.futures
from collections import ChainMap

# Exceptions
class ArgumentError(Exception):
    pass
    
class ConfigurationError(Exception):
    pass
    
class ExecutionError(Exception):
    pass

# Make output unbuffered
sys.stdout = io.TextIOWrapper(open(sys.stdout.fileno(), 'wb', 0), write_through=True)

# Argument handling
parser = argparse.ArgumentParser(
    prog="MFASeqPipeline", 
    description="Process MFA-Seq data")

parser.add_argument("--jobs", default=1, type=int,
                    help="number of parallel jobs to run")
parser.add_argument("--quiet", action=argparse.BooleanOptionalAction,
                    help="do not log to stdout (defaults to true in case of parallel jobs)")
parser.add_argument("--verbose", default=False, action=argparse.BooleanOptionalAction,
                    help="more extensive logging")
parser.add_argument("--steps", type=lambda s: s.split(","),
                    help="comma-separated list of steps to run")
parser.add_argument("config",
                    help="config file in YAML format")
parser.add_argument("targets", nargs="*", metavar="target",
                    help="targets whose jobs to run")
parser.add_argument("settings", nargs="*", metavar="KEY=VALUE",
                    help="override config file settings")

CFG_SEARCHPATH = [ "./", os.path.dirname(os.path.realpath(__file__)) ]

# Set later
CFG = None

# Load config
def load_config(name):
    # Turn config name into a filename
    file_name = f"{name}.yaml" if not name.endswith(".yaml") else name
    # Locate file on search path
    file_path = None
    for dir in CFG_SEARCHPATH:
        file_path = os.path.join(dir, file_name)
        if os.path.exists(file_path):
            break
    if file_path is None:
        raise ConfigurationError(f"unable to location configuration file {file_name}")
    # Load configuration
    with open(file_path, "r") as cfg_file:
        cfg = yaml.safe_load(cfg_file)
    # Recursively load configurations listed in "use"
    if "use" in cfg:
        base_names = cfg["use"]
        if isinstance(base_names, str):
            base_names = [ base_names ]
        if not isinstance(base_names, list) or not all([ isinstance(v, str) for v in base_names ]):
            raise ConfigurationError("'use' must be a string or a list of string of configurations to load")
        for base_name in base_names:
            base = load_config(base_name)
            base.update(cfg)
            cfg = base
    return cfg


# Function to expand configuration settings
def expand(value, cfg=CFG):
    if isinstance(value, str):
        # format_map() uses {{ and }} to represent literal { and }, i.e.
        # double-curly braces are replaced by single braces. This works fine
        # if format_map() is called once, but not if it's called repeatedly
        # until no more replacements are made -- it will eventually reduce
        # any number of curly braces to a single one, and then interpret it,
        # leaving users with no way of producing strings that contain single
        # curly braces. We thus 'quote' {{ and }} by replacing them with {{{{
        # and }}}}, which format_map() then undoes to arrive at the original
        # string. This makes it possible to use {{ and }} to represent { and }
        # despite our iterated calls of format_map()
        vq = value.replace("{{", "{{{{").replace("}}", "}}}}")
        vp = vq.format_map(cfg)
        i = 0
        while value != vp:
            if i > 100:
              raise ConfigurationError("configuration key {key} expands indefinitely")
            value = vp
            # Same quoting of {{ and }} as above
            vq = value.replace("{{", "{{{{").replace("}}", "}}}}")
            vp = vq.format_map(cfg)
            i += 1
    return value.replace("{{", "{").replace("}}", "}")

# Function to execute all requested steps
def execute_job(job):
    # Create job-specific config
    job_cfg = ChainMap(job, CFG)

    # Keep track for each file whether the newest version is in tmp
    tempfiles = {}
    def cleanup_tempfiles():
        for f in tempfiles.values():
            if os.path.exists(f):
                os.remove(f)
    atexit.register(cleanup_tempfiles)
    
    # Execute steps
    try:
        for step_name in job_cfg["steps"]:
            # Within each step, we search for configuration keys first within the
            # step's 'cfg' entry and then at the configuration root. Configuration
            # updates are emphemeral, i.e. lost once the step is completed. This
            # is achieved by putting an empty dictionary which received all updates
            # as the first entry of the dictionary chain we use for lookups and updates
            step = job_cfg[f"{step_name}_step"]
            if not isinstance(step, dict):
                raise ConfigurationError(f"step {step_name} is not properly defined")
            step_cfg = ChainMap({"step": step_name}, step["override"], job_cfg)

            # Get upper-case step name, input and output files, and step script
            input_names = step["input"]
            output_names = step["output"]
            script = step["script"]
            
            # Get expanded input and output file names
            # NOTE: Expansion takes place in the *job* context, since the names
            # are supposed to agree between different steps
            input_files = [expand(job_cfg[f], job_cfg) for f in input_names]
            output_files = [expand(job_cfg[f], job_cfg) for f in output_names]
            
            # Function to write log header listing input and outputs
            def log_header(file):
                title_name = step_name.upper()
                print(f"*** {title_name}", file=file)

                # Report input and output files
                print("    input:", file=file)
                for f in input_files:
                    print(f"         - {f}", file=file)
                print("    output:", file=file)
                for f in output_files:
                    print(f"         - {f}", file=file)

            # Check if the step needs to be executed
            if all([os.path.exists(f) for f in output_files]):
                # Only write to stdout if we"re not running jobs concurrently
                if not ARGS.quiet:
                    log_header(sys.stdout)
                    print(f"    all outputs already exists, skipping step...")
                continue

            # Redirecet stdout to the log file
            if "log_file" in step_cfg:
                log_file = expand(step_cfg["log_file"], step_cfg)
                # Only shadow output to stdout if we're not running jobs concurrently
                if ARGS.quiet:
                    step_stdout = open(log_file, "w+", buffering=1)
                else:
                    step_stdout = subprocess.Popen(["tee", log_file], stdin=subprocess.PIPE,
                                                   text=True, bufsize=1, close_fds=True).stdin
            else:
                step_stdout = sys.stdout

            # Now that we have a log file, write header
            log_header(step_stdout)
            
            # Copy input files to tmp and update config accordingly
            for input_name in input_names:
                # Copy file to tmp dir unless that's where the current version is
                # NOTE: Expansion happens in the *global* context here. For the input file the
                # reason is again that the expansion must agree between steps. For the temp file
                # we could in principle allow step-specific overrides, but then the name of the
                # temp file would depend on the step which created it, which seems silly.
                tmp_dir = expand(job_cfg["tmp_dir"], job_cfg)
                os.makedirs(tmp_dir, exist_ok=True)
                if input_name not in tempfiles:
                    # Create a temp file for the input
                    input_file = expand(job_cfg[input_name], job_cfg)
                    temp_file = os.path.join(tmp_dir, os.path.basename(input_file))
                    # Copy input file to temp file
                    if ARGS.verbose:
                        print(f"    copying {input_file}", file=step_stdout)
                        print(f"    to {temp_file} ...", file=step_stdout)
                    shutil.copyfile(input_file, temp_file)
                    tempfiles[input_name] = temp_file
                # Update path in step-specific cfg
                step_cfg[input_name] = tempfiles[input_name]

            # Generate temp file names for output files
            for output_name in output_names:
                # The output file shouldn't already have a tempfile assigned
                if output_name in tempfiles:
                    temp_file = tempfiles[output_name]
                    raise ExecutionError(f"step output {temp_file} already created by a previous step")
                # Create a temp file for the output
                temp_file = os.path.join(tmp_dir, os.path.basename(expand(job_cfg[output_name], job_cfg)))
                tempfiles[output_name] = temp_file
                # Remove output file should it already exist
                if os.path.exists(tempfiles[output_name]):
                    os.remove(tempfiles[output_name])
                # Update path in step-specific cfg
                step_cfg[output_name] = tempfiles[output_name]
                
            # Execute step
            # NOTE: We use the *local* context in the expansions here, i.e. allow step-specific overrides
            try:
                # Create script file
                (fd, step_cfg["script_file"]) = tempfile.mkstemp(dir=expand(step_cfg["tmp_dir"], step_cfg))
                with os.fdopen(fd, "w") as f:
                    f.write("#!{shell}\n".format_map(step_cfg))
                    f.write(expand(script, step_cfg))
                os.chmod(step_cfg["script_file"], 0o755)
                # Echo script for debugging purposes
                if ARGS.verbose:
                    print(f"    step script:", file=step_stdout)
                    with open(step_cfg["script_file"], "r") as f:
                        for l in f.readlines():
                            l = l.rstrip()
                            print(f"    | {l}", file=step_stdout)
                
                # Run script
                execute_cmd = expand(step_cfg["execute_cmd"], step_cfg)
                if ARGS.verbose:
                    print(f"    executor: {execute_cmd}", file=step_stdout)
                    print(f"    running step ...", file=step_stdout)
                subprocess.run(execute_cmd, shell=True, executable=step_cfg["shell"],
                               stdout=step_stdout, stderr=step_stdout, check=True)
            finally:
                os.remove(step_cfg["script_file"])
            if ARGS.verbose:
                print(f"    ... step completed", file=step_stdout)
            
            # Copy output files from tmp
            for output_name in output_names:
                temp_file = tempfiles[output_name]
                output_file = expand(job_cfg[output_name], job_cfg)
                if not os.path.exists(temp_file):
                    raise ExecutionError(f"step failed to produce {temp_file}")
                if ARGS.verbose:
                    print(f"    copying {temp_file}", file=step_stdout)
                    print(f"    to {output_file} ...", file=step_stdout)
                dir_name = os.path.dirname(output_file)
                if len(dir_name) > 0:
                    os.makedirs(dir_name, exist_ok=True)
                shutil.copy(temp_file, output_file)

            if ARGS.verbose:
                print(f"    done.", file=step_stdout)
    finally:
        # Remove temporary files
        cleanup_tempfiles()
        atexit.unregister(cleanup_tempfiles)

try:
    # Parse arguments
    # NOTE: The parser cannot actually distinguish between target and
    # configuration overrides, and thus puts all extra arguments in targets.
    # We thus must split the targets manually into actual targets and seettings.
    ARGS = parser.parse_args()

    if ARGS.quiet is None:
        ARGS.quiet = (ARGS.jobs > 1)

    # Load config
    CFG = load_config(ARGS.config)

    # Merge command-line arguments into config
    if ARGS.steps is not None:
        CFG["steps"] = ARGS.steps

    # Split ARGS.targets into actual targets and settings, apply setting to config
    targets_or_settings = ARGS.targets
    in_targets = True
    ARGS.targets = []
    ARGS.settings = []
    for i in range(0, len(targets_or_settings)):
        arg = targets_or_settings[i]
        in_targets = in_targets and "=" not in arg
        if in_targets:
            if arg == "all":
                ARGS.targets.extend(CFG["targets"].keys())
            else:
                ARGS.targets.append(arg)
        else:
            key_value = arg.split("=", maxsplit=1)
            if len(key_value) < 2:
                raise ArgumentError(f"invalid argument {arg}, must be of the form KEY=VALUE")
            try:
                key, value = key_value[0], yaml.safe_load(io.StringIO(key_value[1]))
            except yaml.parser.ParserError as e:
                raise ArgumentError(f"invalid YAM syntax in: {key_value[0]}={key_value[1]}\n{e}") from e
            CFG[key] = value

    # Execute all jobs in target
    if len(ARGS.targets) > 0:
        POOL = concurrent.futures.ThreadPoolExecutor(max_workers=ARGS.jobs)
        for target in ARGS.targets:
            if target not in CFG["targets"]:
                raise ArgumentError(f"invalid target {target}")
            if not ARGS.quiet:
                print(f'=== TARGET: {target}')
            for r in POOL.map(execute_job, CFG["targets"][target]):
                pass
            if not ARGS.quiet:
                print(f'=== COMPLETED TARGET: {target}')
    else:
        execute_job({})
except KeyError as e:
    print(f"ERROR: Variable {e} is not defined", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)
