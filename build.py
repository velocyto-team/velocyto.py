import shutil
import sys
from pathlib import Path

import numpy as np
from Cython.Build import cythonize
from Cython.Compiler.Errors import CompileError
from Cython.Distutils.build_ext import build_ext
from setuptools.dist import Distribution
from setuptools.extension import Extension

# C Extensions

C_COMPILE = True
USE_CYTHON = True

if C_COMPILE:
    package_data: dict = {}
    if USE_CYTHON:
        extensions = [
            Extension(
                "velocyto.speedboosted",
                include_dirs=[np.get_include()],
                sources=["src/velocyto/speedboosted.pyx"],
                extra_compile_args=[
                    "-fopenmp",
                    "-ffast-math",
                    "-O3"
                ],  # NOTE optional flags -O3 -ffast-math -march=native
                extra_link_args=["-fopenmp"],
            )
        ]
        extensions = cythonize(
            extensions,
            include_path=[np.get_include()],
            compiler_directives={"language_level": 3}
            )
    else:
        extensions = [
            Extension(
                "velocyto.speedboosted",
                ["src/velocyto/speedboosted.c"],
                extra_compile_args=["-fopenmp", "-ffast-math"],
                extra_link_args=["-fopenmp"],
            )
        ]
else:
    extensions = []
    if "darwin" in sys.platform:
        compiled = "src/velocyto/speedboosted.cpython-36m-darwin.so"
    elif "win" in sys.platform:
        sys.stdout.write("Sorry, we do not support (or like) Windows OS")
        sys.exit()
    elif "linux" in sys.platform:
        compiled = "src/velocyto/speedboosted.cpython-36m-x86_64-linux-gnu.so"
    package_data = {"velocyto": [compiled]}


class BuildFailed(Exception):
    pass


class ExtBuilder(build_ext):
    # This class allows C extension building to fail.

    built_extensions = []

    def run(self):
        try:
            build_ext.run(self)
        except (FileNotFoundError):
            print(
                "  Unable to build the C extensions, "
                "velocyto will use the pure python code instead."
            )

    def build_extension(self, ext):
        try:
            build_ext.build_extension(self, ext)
        except (CompileError, ValueError):
            print(
                f"  Unable to build the '{ext.name}' C extension, velocyto "
                f"will use the pure python version of the extension."
            )


def build(setup_kwargs):
    """
    This function is mandatory in order to build the extensions.
    """
    distribution = Distribution(
        {
            "name": "src/velocyto",
            "ext_modules": cythonize(
                extensions, compiler_directives={"language_level": "3"}
            ),
        }
    )
    distribution.package_dir = "velocyto"

    cmd = ExtBuilder(distribution)
    cmd.ensure_finalized()
    cmd.run()

    # Copy built extensions back to the project
    for output in cmd.get_outputs():
        output = Path(output)
        relative_extension = Path("src").joinpath(output.relative_to(cmd.build_lib))
        if not output.exists():
            continue

        shutil.copyfile(output, relative_extension)
        mode = relative_extension.stat().st_mode
        mode |= (mode & 0o444) >> 2
        relative_extension.chmod(mode)

    return setup_kwargs


if __name__ == "__main__":
    build({})
