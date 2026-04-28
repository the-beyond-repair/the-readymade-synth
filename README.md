# The Readymade Synth: Reference Implementation

A minimal SuperCollider implementation of *Prepared Synthesis*, the shape-to-sound method introduced in:

> Antoni Rayzhekov and Martin Murer. 2026. The Readymade Synth: Prepared Synthesis with Everyday Objects. In *Proceedings of the International Conference on New Interfaces for Musical Expression* (NIME '26), June 23–26, 2026, London, UK. ACM, New York, NY, USA, 9 pages. https://doi.org/XX.XXXX/XXXXXXX.XXXXXXX

This repository contains the minimal implementation and algorithmic core of the method, applied to synthetic SVG shapes, so that the analysis pipeline is reproducible and inspectable. It is not the full instrument described in the paper. The computer-vision tracker, MIDI integration, and hardware controls are not included here and will be the subject of a separate publication in the future.

## What it does

A closed 2D contour is treated as a complex periodic signal `z[n] = x[n] + iy[n]`. Its FFT yields positive (CCW) and negative (CW) frequency components, which are mapped to harmonic ratios for an additive bank. Object area drives the envelope; centroid position drives stereo placement.

Pipeline: `SVG polygon → contour resampling → complex FFT → spectral analysis → Klang additive synthesis`.

## Requirements

SuperCollider 3.x, with a working audio server. SVG `<polygon>` files in the `svg/` directory. A small sample set (`circle.svg`, `ellipse.svg`, `gear.svg`, `leaf.svg`, `pebble.svg`, and others) is provided.

## Usage

Open `the_ready_made_synth-minimal_implementation.scd` in the SuperCollider IDE and evaluate the parenthesised boot block. Once the post window prints `Ready.`, call:

```supercollider
~loadAndPlay.("ellipse.svg", 220);
~loadAndPlay.("gear.svg", 220);
~loadAndPlay.("leaf.svg", 220);
.....
```

The second argument is the fundamental frequency in Hz. Loading different shapes at the same frequency makes timbre comparison direct.

## Configuration

Two parameters control analysis behaviour:

```supercollider
~offsetMode = \phase;   // \fixed, \phase, or \symmetry
~normMode   = \joint;   // \joint or \independent
```

`~offsetMode` selects the interleaving offset between positive and negative spectral halves (paper Section 3.3). The phase-derived offset is the default and what the paper adopts.

`~normMode` selects the amplitude normalisation strategy (paper end of Section 3.2). `\joint` divides by the global peak across both halves, preserving the natural pos/neg amplitude ratio and yielding cleaner, more stable timbres. `\independent` divides each half by its own peak, equalising CW and CCW components. Both modes share a single global peak-picking step that retains the strongest `~maxPartials` (default 64) bins from the pooled spectrum.

## License

GPL 3.0. See `LICENSE`.
