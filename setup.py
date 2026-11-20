# setup.py

import shutil
import subprocess
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


ROOT = Path(__file__).resolve().parent


class MakeExtension(Extension):
    def __init__(self, name):
        super().__init__(name, sources=[])


class MakeBuildExt(build_ext):

    def run(self):
        subprocess.run(
            ["make", f"PYTHON={sys.executable}"],
            cwd=ROOT,
            check=True,
        )

        super().run()

    def build_extension(self, ext):
        source = ROOT / self.get_ext_filename(ext.name)

        if not source.is_file():
            raise RuntimeError(f"Compiled extension not found: {source}")

        destination = Path(self.get_ext_fullpath(ext.name))
        destination.parent.mkdir(parents=True, exist_ok=True)

        if source.resolve() != destination.resolve():
            shutil.copy2(source, destination)


setup(
    packages=["numerov_qex"],
    ext_modules=[MakeExtension("numerov_qex.radialsolver")],
    cmdclass={"build_ext": MakeBuildExt},
)