# AnasenAnalyzer

AnasenAnalyzer is a re-write of the classic (pre-digital-upgrade) analysis toolkit for an ANASEN data set. It aims to solve several problems related to the original system from configuration management to physics implementation and precision control.

Some implementation details are restricted based on legacy code. The input data file structure is assumed to be the classic post-calibration ANASEN scheme, which is not well optimized, and should be modified if a completely new dataset is to be used.

AnasenAnalyzer contains all dependencies in git submodules execpt for the ROOT framework, which is required to be installed.

## Building and Running

To dowload the code base use `git clone --recursive https://github.com/gwm17/AnasenAnalyzer.git`, which will pull in the analyzer code as well as the dependencies. AnasenAnalyzer uses CMake as a build system, so standard operating procedures there. As an example, to build AnasenAnalyzer on a Unix(MacOS/Linux) system use (from the repository directory):

```bash
mkdir build
cd build
cmake ..
make -j 4
```

By default, the binary executable is placed in the repository bin directory, and should be run from the repository as `./bin/AnasenAnalyzer <your_input_file>`, where `<your_input_file>` is a text file containing configuration details (an example is given in input.txt).

## What is different from the old analysis?

### Energy loss

For active target detectors, energy loss through the active target is a critical part of the analysis, but often one where many compromises are made. In the older analysis codes, energy loss was handled by interpolating sanitized SRIM tables, which often had to be "adjusted" to match LISE++ energy loss calcualtions. This lead to several points of confusion (how to adjust to match LISE, too many SRIM files lying around, finding a machine that can actually run SRIM, etc). This code base uses the CAtima library (with some extensions for reverse-integration) to calculate the energy loss, effectively cutting out the middle-man. CAtima is a C++ library based on the Fortan ATIMA library, which is the energy loss library used by LISE++. This means that the analyzer has (relatively) good agreement with LISE and other such calculators out of the box. It also means that specifiying an active target is as simple as defining the target material (chemical composition, density), whereas previously it involved generating a SRIM file for *each* active target-reactant pair. The cost associated with this change is mostly in computation speed. The spline method was cumbersome at configuration time, but very efficient at runtime. CAtima has many tricks to efficiently calculate energy loss *forward*, but reverse integrating energy loss over a path is not so simple and especially for low energy particles (~100s of keV) and dense targets computation times can be longer (a few milliseconds per path). However, the ease of configuration and high accuracy are well worth the trade-off in such a critical area of the analysis.

### Output structure

Previously data was outputed in both a ROOT TTree (using a custom dictionary) and histograms. This can be useful, but in this case it was mostly overkill. Outputing a ROOT tree implies that there will be subsequent analysis stages, which really shouldn't be the case (this is the analyzer afterall). Using trees to prototype analyses can be useful, however they can also lead to results which can not be reproduced or are not clearly defined, as the cuts/analysis sequence is not saved any where. To avoid these pitfalls, AnasenAnalyzer *only* writes histograms and cuts to disk (this also saves some computation time).

### Precision

An often overlooked aspect of analysis is the precision of comparisions of equality between floating point representations of physical values. AnasenAnalyzer has a namespace Precision which allows for well-defined comparisions of floating point values with custom precision specifications. This is particularly important when considering detector geometry and/or gating.

### Configuration

A common issue with analysis tools is how to balance generalization with specific use-cases. AnasenAnalyzer attempts to find a middle ground using configuration files. The configuration files follow a simple script format to specifiy target materials, input data files, 2D graphical cuts (ROOT TCutG), 1D gates (range of values), analysis switches, and output files. This allows for users to save multiple configurations defining various analyses. However, not everything can be specified at runtime. ANASEN has so many use cases that complete generality is impossible. It is left to users to define things like reaction channels, beam energies, etc in the code itself. Examples are included for the analysis of 7Be+d->a+a+p in ANASEN (which should cover a lot of use cases). The Reconstructor class may need to be extended to handle more unqiue analysis methods.

## Requirements

CMake 3.12 or greater
ROOT 6.22 or greater (tested with 6.28), needs to have ROOT CMake extensions for generating dictionaries
