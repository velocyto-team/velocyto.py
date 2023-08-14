# Changelog

** [0.29.0] - 2023-08-07

*** Fixed

- LOOKS LIKE VELOCYTO IS BACK ON THE MENU, BOYS!

*** Changed

- Updated to require Python 3.10 as minimum version
- Updated to use Typer 0.9.0
- Updated dependencies
- Pulled a lot out of `velocyto.commands._run()` into separate functions
- Removed ability to add UMAP and cluster info to the loom file - that will be regenerated in other software, so why?

** [0.28.0] - 2023-02-24

*** Changed

- Code cleanup and improvements thanks to sourcery

** [0.27.0] - 2023-02-23

*** Changed

- Pass the actual "samtools_memory" value to samtools instead of checking the amount of memory available.
  - Also change to passing the argument as a string, allowing for the use of a prefix instead of requireing memory to be
  defined in megabytes.

** [0.26.0] - 2023-02-22

*** Changed

- Various code cleanups
- Add ruff as package linter

## [0.20.0] - [0.25.0]

Various. No longer sure but a lot.

## [0.19.1] - 2022-06-30

### Fixed

- Fixed version issues

## [0.19.0] - 2022-06-29

### Changed

- Removed dumping all imports into __init__.py and then just importing the entire module (which lead to circular imports) and am importing as needed
- removed wildcard importing
- removed unused imports
- commented out a lot of unused assigned variables

### Fixed

- lots of changes to fix things flagged by flake8


## [0.18.0] - 2022-06-29


### Changed

- Black formatting
- Replaced setup.py with pyproject.toml
- remove distutils from build
- lock requirements and versions

[0.27.1]: https://github.com/olivierlacan/keep-a-changelog/releases/tag/0.27.0
[0.19.1]: https://github.com/olivierlacan/keep-a-changelog/releases/tag/0.19.0..0.19.1
[0.19.0]: https://github.com/olivierlacan/keep-a-changelog/releases/tag/0.18.0..0.19.0
[0.18.0]: https://github.com/olivierlacan/keep-a-changelog/releases/tag/0.18.0