# SPR Kinetics Analysis

A web-based application for analyzing Surface Plasmon Resonance (SPR) binding kinetics data. Features a high-performance C fitting engine with a Python/HTML frontend.

## Features

- **Multi-Cycle Kinetics (MCK)** and **Single-Cycle Kinetics (SCK)** analysis
- **Auto-detection** of data format (MCK vs SCK), association/dissociation timing, and replicate experiments
- **Three binding models**: 1:1 Langmuir, Heterogeneous Ligand, Two-State Conformational Change
- **Global fitting** with Nelder-Mead simplex optimization, RK4 numerical integration, and multi-threaded random restarts
- **Steady-state affinity** analysis with KD determination
- **Biacore export compatibility**: handles variable column formats, fitted overlays, excluded cycles, and duplicate replicates
- **Interactive web UI**: drag-and-drop upload, cycle exclusion toggles, real-time re-fitting, zoomable plots
- **PPTX export**: publication-ready slides (Sensorgram, Steady-State, Residuals, Parameters)
- **Reference subtraction** support

## Requirements

- **C compiler** (GCC or Clang with C11 support)
- **Python 3.6+**
- **python-pptx** (`pip install python-pptx`) — for PPTX export
- **matplotlib** (`pip install matplotlib`) — for PPTX plot generation

## Quick Start

```bash
# Build the C fitting engine
make

# Start the web server
python3 spr_server.py

# Open in browser
open http://localhost:8765
```

## Usage

1. **Upload data**: Drag and drop one or more Biacore export files (.txt) into the upload zone
2. **Optional reference**: Upload a reference channel file for subtraction
3. **Analyze**: The app auto-detects format and timing, then performs global fitting
4. **Refine**: Toggle individual cycles on/off, switch models, and re-fit
5. **Export**: Download results as a PPTX presentation

## Architecture

```
spr_server.py      — Python HTTP server, serves UI and orchestrates fitting
spr_fit            — C binary that performs the actual curve fitting
  spr_fit_main.c   — CLI entry point, JSON output
  spr_io.c/.h      — Data I/O, MCK/SCK parsers, format auto-detection
  spr_models.c/.h  — ODE models (Langmuir, Heterogeneous, Two-State)
  spr_optim.c/.h   — Nelder-Mead optimizer, steady-state analysis
  spr_types.h      — Shared data structures and constants
Makefile           — Build configuration
```

## Input Format

Currently accepts tab-delimited text files exported from **Biacore Evaluation Software**. The parser handles:

- Standard MCK format with `_X` / `_Y` column pairs per cycle
- Instrument-fitted overlay columns (`; Fitted_X` / `; Fitted_Y`)
- Excluded cycles (marked in headers)
- Duplicate/triplicate replicates at the same concentration
- SCK staircase injection traces

## Example Data

The `example_data/` directory contains sample MCK datasets (FUR-CA2 series) that can be used for testing. These include files with and without instrument-fitted overlays, as well as files with duplicate replicates.

## License

MIT
