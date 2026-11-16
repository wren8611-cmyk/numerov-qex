# setup.py

import os
import sys
from pathlib import Path

from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext


# ------------------------------------------------------------
# Platform
# ------------------------------------------------------------

if sys.platform != "darwin" and not sys.platform.startswith("linux"):
    raise RuntimeError(
        "numerov-qex currently supports only macOS and Linux."
    )


# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------

def env_path(name):
    value = os.environ.get(name)
    return Path(value).expanduser() if value else None


def find_directory(candidates, marker):
    for directory in candidates:
        if directory is None:
            continue

        directory = Path(directory).expanduser()

        if (directory / marker).exists():
            return str(directory.resolve())

    return None


def find_library_directory(candidates, library):
    extensions = (".so", ".so.*", ".dylib", ".a")

    for directory in candidates:
        if directory is None:
            continue

        directory = Path(directory).expanduser()

        if not directory.is_dir():
            continue

        for extension in extensions:
            if list(directory.glob(f"lib{library}{extension}")):
                return str(directory.resolve())

    return None


# ------------------------------------------------------------
# Eigen
# ------------------------------------------------------------

def find_eigen():

    root = env_path("EIGEN_ROOT")
    explicit = env_path("EIGEN3_INCLUDE_DIR")

    candidates = []

    if explicit:
        candidates.append(explicit)

    if root:
        candidates.extend([
            root,
            root / "include",
            root / "include" / "eigen3",
        ])

    candidates.extend([
        "/usr/include/eigen3",
        "/usr/local/include/eigen3",
        "/usr/local/include",
        "/opt/homebrew/include/eigen3",
        "/opt/homebrew/include",
        "/opt/local/include/eigen3",
        "/opt/local/include",
    ])

    result = find_directory(candidates, "Eigen/Dense")

    if result is None:
        raise RuntimeError(
            "Eigen3 not found.\n"
            "Set EIGEN_ROOT or EIGEN3_INCLUDE_DIR."
        )

    return result


# ------------------------------------------------------------
# GSL
# ------------------------------------------------------------

def find_gsl():

    root = env_path("GSL_ROOT")
    explicit_include = env_path("GSL_INCLUDE_DIR")
    explicit_lib = env_path("GSL_LIBRARY_DIR")

    include_candidates = []
    lib_candidates = []

    if explicit_include:
        include_candidates.append(explicit_include)

    if explicit_lib:
        lib_candidates.append(explicit_lib)

    if root:
        include_candidates.append(root / "include")
        lib_candidates.extend([
            root / "lib",
            root / "lib64",
        ])

    include_candidates.extend([
        "/usr/include",
        "/usr/local/include",
        "/opt/homebrew/include",
        "/opt/local/include",
    ])

    lib_candidates.extend([
        "/usr/lib",
        "/usr/lib64",
        "/usr/local/lib",
        "/usr/local/lib64",
        "/opt/homebrew/lib",
        "/opt/local/lib",
        "/usr/lib/x86_64-linux-gnu",
        "/usr/lib/aarch64-linux-gnu",
    ])

    include_dir = find_directory(
        include_candidates,
        "gsl/gsl_errno.h",
    )

    lib_dir = find_library_directory(
        lib_candidates,
        "gsl",
    )

    if include_dir is None:
        raise RuntimeError(
            "GSL headers not found.\n"
            "Set GSL_ROOT or GSL_INCLUDE_DIR."
        )

    if lib_dir is None:
        raise RuntimeError(
            "GSL libraries not found.\n"
            "Set GSL_ROOT or GSL_LIBRARY_DIR."
        )

    return include_dir, lib_dir


# ------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------

eigen_include = find_eigen()
gsl_include, gsl_lib = find_gsl()

print(f"Eigen include : {eigen_include}")
print(f"GSL include   : {gsl_include}")
print(f"GSL library   : {gsl_lib}")


# ------------------------------------------------------------
# Extension
# ------------------------------------------------------------

ext_modules = [
    Pybind11Extension(
        "numerov_qex.radialsolver",
        ["src/radialsolver.cxx"],
        include_dirs=[
            "include",
            eigen_include,
            gsl_include,
        ],
        library_dirs=[
            gsl_lib,
        ],
        libraries=[
            "gsl",
            "gslcblas",
        ],
        extra_compile_args=[
            "-O3",
        ],
        cxx_std=17,
    ),
]


# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------

setup(
    packages=["numerov_qex"],
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)