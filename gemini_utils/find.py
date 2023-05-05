"""
functions for finding files
"""

from __future__ import annotations
from datetime import datetime, timedelta
from pathlib import Path
import shutil
import os
import subprocess

import numpy as np

EXE_PATHS = [
    ".",
    "bin",
    "build",
    "build/bin",
    "build/Debug",
    "build/RelWithDebInfo",
    "build/Release",
]


def filename2datetime(path: Path) -> datetime:
    """
    Gemini3D datafiles use a file naming pattern that we translate to a datetime

    path need not exist.
    """

    name = path.name

    return datetime.strptime(name[:8], "%Y%m%d") + timedelta(seconds=float(name[9:21]))


def config(path: Path) -> Path:
    """given a path or config filename, return the full path to config file"""

    return find_stem(path, stem="config", suffix="nml")


def simsize(path: Path) -> Path:
    """gets path to simsize file"""

    return find_stem(path, stem="simsize")


def gemini_exe(exe: str = None) -> Path:
    """
    find and check that Gemini executable can run on this system
    """

    name = "gemini3d.run"

    if not exe:  # allow for default dict empty
        root = gemini_root()
        for n in EXE_PATHS:
            e = shutil.which(name, path=str(root / n))
            if e:
                break
    if not e:
        raise EnvironmentError(
            f"{name} not found."
            "Set environment variable GEMINI_ROOT to directory above bin/gemini.bin"
        )

    # %% ensure Gemini3D executable is runnable
    gemexe = Path(e).expanduser()
    ret = subprocess.run(
        [str(gemexe)],
        capture_output=True,
        timeout=10,
        text=True,
        cwd=gemexe.parent,
    )
    if ret.returncode == 0:
        pass
    elif ret.returncode == 3221225781 and os.name == "nt":
        # Windows 0xc0000135, missing DLL
        raise RuntimeError(
            "On Windows, it's best to build Gemini3D with static libraries--including all numeric libraries "
            "such as LAPACK.\n"
            "Currently, we are missing a DLL on your system and gemini.bin with shared libs cannot run."
        )
    else:
        raise EnvironmentError(
            f"\n{gemexe} was not runnable on your platform--try rebuilding:\n"
            "gemini3d.setup()\n"
            f"{ret.stderr}"
        )

    return gemexe


def gemini_root() -> Path:
    root = os.environ.get("GEMINI_ROOT")
    if not root:
        raise EnvironmentError(
            "Please set environment variable GEMINI_ROOT to (desired) top-level Gemini3D directory."
            "If Gemini3D is not already there, PyGemini will download and build Gemini3D there."
        )
    return Path(root).expanduser()


def msis_exe(root: Path = None) -> Path | None:
    """
    find MSIS_SETUP executable
    """

    name = "msis_setup"

    paths = [root, os.environ.get("CMAKE_PREFIX_PATH"), os.environ.get("GEMINI_ROOT")]
    paths = [Path(p).expanduser() for p in paths if p]

    if not paths:
        raise EnvironmentError(
            "Specify location of msis_setup executable by environment variable"
            " GEMINI_ROOT or CMAKE_PREFIX_PATH or give gemini_root argument"
        )
    for path in paths:
        for n in EXE_PATHS:
            msis_exe = shutil.which(name, path=str(path / n))
            if msis_exe:
                return Path(msis_exe)

    return None


def frame(simdir: Path, time: datetime) -> Path:
    """
    find frame closest to time
    """

    suffix = ".h5"

    simdir = Path(simdir).expanduser()

    stem = (
        time.strftime("%Y%m%d")
        + f"_{time.hour*3600 + time.minute*60 + time.second:05d}."
        + f"{time.microsecond:06d}"
    )

    fn = simdir / (stem + suffix)
    if fn.is_file():
        return fn

    # %% WORKAROUND for real32 file ticks. This will be removed when datetime-fortran is implemented
    MAX_OFFSET = timedelta(seconds=1)  # 10 ms precision, allow extra accumulated tolerance
    pat = time.strftime("%Y%m%d") + "_*"

    file_times = []
    files = list(simdir.glob(pat + suffix))
    for fn in files:
        file_times.append(filename2datetime(fn))

    if file_times:
        afile_times = np.array(file_times)
        i = abs(afile_times - time).argmin()  # type: ignore

        if abs(afile_times[i] - time) <= MAX_OFFSET:
            return files[i]

    raise FileNotFoundError(f"{stem}{suffix} not found in {simdir}")


def grid(path: Path) -> Path:
    """given a path or filename, return the full path to simgrid file"""

    return find_stem(path, stem="simgrid")


def find_stem(path: Path, stem: str, suffix: str = ".h5") -> Path:
    """find file containing stem"""

    path = Path(path).expanduser()

    if path.is_file():
        if stem in path.stem:
            return path
        else:
            found = find_stem(path.parent, stem, path.suffix)
            if not found:
                raise FileNotFoundError(f"{stem} not found in {path.parent}")
            return found

    if isinstance(suffix, str):
        if not suffix.startswith("."):
            suffix = "." + suffix
        suffixes = [suffix]
    else:
        suffixes = suffix

    if path.is_dir():
        for p in (path, path / "inputs"):
            for suff in suffixes:
                f = p / (stem + suff)
                if f.is_file():
                    return f

    raise FileNotFoundError(f"{stem} not found in {path}")


def inputs(direc: Path, input_dir: Path = None) -> Path:
    """
    find input parameter directory

    direc: pathlib.Path
        top level simulation dir
    input_dir: pathlib.Path
        relative to top level, the parameter dir. E.g. inputs/precip, input/Efield
    """

    direc = Path(direc).expanduser()
    if input_dir:
        input_path = Path(input_dir).expanduser()
        if not input_path.is_absolute():
            input_path = direc / input_path
    else:
        input_path = direc

    if not input_path.is_dir():
        raise NotADirectoryError(direc)

    return input_path
