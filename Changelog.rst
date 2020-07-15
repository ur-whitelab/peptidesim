Change Log
==========

v0.5 (2020-07-15)
=======
-----------------------

*New Features*

- Automated creation of cs2backbone data directory
- Reweighting analysis functions

*Enhancements*

- Added tutorials
- Added EDS load/plot utilities
- Moved some methods from scripts to be standard
- No longer generates files in root directory
- More job interrupt resiliency
- Allowed kwargs to be updated when restarting
- New input examples


v0.4 (2020-03-06)
-----------------------

*Enhancements*

- Pickling is handled by peptidesim
- Tests run significantly faster on bluehive
- You can now print peptidesim objects
- Removed no longer used code

v0.3 (2020-02-28)
-----------------------

*New Features*

- Added method assert chemical shift validity for biasing

*Enhancements*

- Updated to python3
- Added documentation examples
- Updated highest supported gromacs version to 2019 and plumed 2.6

*Bug Fixes*

- CPT files are now synched when doing replica-exchange
- Dependencies are now set and pinned
- EDS Bias restarts can work now beyond 2 restarts
