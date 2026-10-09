"""Exercise real space-containing paths and wrapper failure/lock handling."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def run(command, work, env, success=True):
    result = subprocess.run(command, cwd=str(work), env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True, timeout=30)
    if (result.returncode == 0) != success:
        raise AssertionError(result.stdout)
    return result.stdout


def executable(path, script):
    path.write_text("#!/bin/sh\n" + script)
    path.chmod(0o755)


def main():
    wrapper, data, root = (Path(x).resolve() for x in sys.argv[1:])
    root.mkdir(parents=True, exist_ok=True)
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("RABBIT_", "RB_", "OMP_"))
           and k not in ("BADMAP", "PCTID", "MINDEPTH")}
    env["RB_LABEL"] = "space label"
    with tempfile.TemporaryDirectory(prefix="paths ", dir=str(root)) as tmp:
        work = Path(tmp)
        assembly, bam = work / "assembly space.fa", work / "sample one.bam"
        shutil.copyfile(str(data / "contigs.fa"), str(assembly))
        shutil.copyfile(str(data / "contigs-1000.fastq.bam"), str(bam))
        command = [str(wrapper), "--threads", "2", "--seed", "42",
                   "--min-bin-size", "0", str(assembly), str(bam)]
        run(command, work, env)
        members = list(work.glob("*.rabbitbin-*/out.members.tsv"))
        if len(members) != 1 or len(members[0].read_text().splitlines()) != 2:
            raise AssertionError("space-containing paths did not produce members")
        env["RB_LABEL"] = "reuse label"
        if "Using existing depth file" not in run(command, work, env):
            raise AssertionError("existing depth was not reused")

    for case in ("depth_failure", "non_lock_file", "foreign_lock", "option_spaces"):
        with tempfile.TemporaryDirectory(prefix=case, dir=str(root)) as tmp:
            work = Path(tmp)
            shutil.copy2(str(wrapper), str(work / "run_rabbitbin.sh"))
            executable(work / "rabbitbin", 'if [ "$1" = --help ]; then exit 0; fi\n'
                       'printf "%s\\n" "$@" > bin-arguments.txt\n')
            executable(work / "rabbit_depth", 'if [ "$#" = 0 ]; then exit 0; fi\n'
                       'exit 42\n')
            # Fail a wait immediately so the foreign-lock test never sleeps.
            executable(work / "sleep", "exit 1\n")
            (work / "assembly.fa").write_text(">fixture\nACGT\n")
            (work / "sample.bam").touch()
            lock = work / "assembly.fa.depth.txt.BUILDING"
            if case == "non_lock_file":
                lock.write_text("not our lock\n")
            elif case == "foreign_lock":
                lock.symlink_to("another writer")
            elif case == "option_spaces":
                (work / "assembly.fa.depth.txt").write_text("cached fixture\n")
            run([str(work / "run_rabbitbin.sh"), "--test-option", "one two",
                 "assembly.fa", "sample.bam"], work, env,
                success=(case == "option_spaces"))
            if case == "option_spaces":
                args = (work / "bin-arguments.txt").read_text().splitlines()
                if args[:2] != ["--test-option", "one two"]:
                    raise AssertionError("wrapper split a quoted option value")
            elif (work / "bin-arguments.txt").exists():
                raise AssertionError("binning started after a depth/lock failure")
            if case == "foreign_lock" and not lock.is_symlink():
                raise AssertionError("removed another writer's lock")
            if case == "non_lock_file" and lock.read_text() != "not our lock\n":
                raise AssertionError("changed an existing non-lock file")
            if case == "depth_failure" and lock.is_symlink():
                raise AssertionError("failed to release our own lock")
    print("wrapper paths, quoted options, cache reuse and failure handling passed")


if __name__ == "__main__":
    main()
