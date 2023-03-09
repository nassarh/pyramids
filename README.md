# pyramids
construct surfaces tessellated with pyramids

![Example - equilateral pyramids tessellating a surface with skew symmetry](pyramidGeneral.png)

## Preliminaries
* A *pyramid* is composed of four triangular lateral faces and one quadrangular base.
* The base has equal opposite edges and is immaterial.
* The lateral triangles are rigid, hinged at edges where they meet.
* Each pyramid can "close" and "open" as a single degree-of-freedom mechanism.
* A *tessellation* is a 2D grid of pyramids that are identical, initially, but that are opened or closed to various degrees.

## Main Observation
* Let one array of pyramids taken in a "diagonal" direction be given. The next array can be constructed by sphere intersection.

## Contents
* ``pyramids.py`` contains functions that allow to construct polyhedral surfaces that are "screw symmetric" deformed tessellations.
* The remaining files are examples.

## How it works
* The input is the deformed state of one pyramid and the rotation matrix that defines the "screw symmetry".
* Letting the rotation act on the pyramid produces the first array of pyramids taken along a "diagonal".
* The next array is found by ensuring that edges maintain their lengths. This involves iterating a function that constructs the intersection of three spheres.
* The construction is pursued as long as said spheres have an intersection. If they don't, the construction stops.

## Dependencies
* ``numpy`` and ``pyvista``

## References
* https://hal.science/hal-01368009v1/file/MainDocument.pdf
* https://arxiv.org/pdf/2207.08752
* https://hal-enpc.archives-ouvertes.fr/hal-01978795v1/file/7OSME-paper-Nassar-Lebee-Monasse-8.pdf
* https://hal-enpc.archives-ouvertes.fr/hal-01691183v1/file/IASS17_revised.pdf

## Acknowledgments
* https://www.nsf.gov/awardsearch/showAward?AWD_ID=2045881
* http://mmcd.univ-paris-est.fr/funded-projects/post-doc-projects/

## Terms
* Author: Hussein Nassar (nassarh@missouri.edu)
* You are free to use and modify for research and education purposes with proper citation and attribution.
* For commercial use, contact the author.
