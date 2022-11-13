# Changelog

** [0.24.0] - 2023-02-22

*** Changed

- Various code cleanups
- Add ruff as package linter

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


[0.19.1]: https://github.com/olivierlacan/keep-a-changelog/releases/tag/0.19.0..0.19.1
[0.19.0]: https://github.com/olivierlacan/keep-a-changelog/releases/tag/0.18.0..0.19.0
[0.18.0]: https://github.com/olivierlacan/keep-a-changelog/releases/tag/0.18.0