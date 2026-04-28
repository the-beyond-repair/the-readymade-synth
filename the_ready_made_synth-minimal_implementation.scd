// ================================================================
// The Readymade Synth: Reference Implementation
//
// A closed 2D contour is treated as a complex periodic signal (x+iy).
// Its DFT yields positive (CCW) and negative (CW) frequency bins —
// the Fourier ellipse decomposition of the shape. These map to
// frequency ratios for additive synthesis, with an interleaving
// offset that ties geometric symmetry to spectral harmonicity.
//
// Pipeline:  SVG polygon → geometry → complex FFT → Klang synthesis
// Requires:  SVG <polygon> files in sc/svg/
// ================================================================

(
s.waitForBoot {

	~svgDir = PathName(thisProcess.nowExecutingPath).pathOnly +/+ "svg";
	~maxPartials = 64;
	~minDbAmp    = -55;
	~offsetMode  = \phase; // \symmetry, \phase, or \fixed
	~normMode    = \joint; // \joint (preserves pos/neg ratio) or \independent (equalises CW/CCW)

	// ── 1. SVG PARSING ──────────────────────────────────────────
	// Extract [[x,y], ...] from <polygon points="..."/>

	~parseSVG = { |path|
		var text = File.readAllString(path);
		var tagStart = text.find("<polygon");
		if(tagStart.notNil) {
			var attrStart;
			var attrEnd;
			attrStart = text.find("points=", false, tagStart) + 8;
			attrEnd = text.find("\"", false, attrStart);
			text[attrStart..attrEnd-1].split($ ).reject{|c| c.size == 0}.collect { |pair|
				var xy = pair.split($,);
				[xy[0].asFloat, xy[1].asFloat]
			};
		};
	};

	// ── 2. CONTOUR PROCESSING ───────────────────────────────────
	// Center on origin, extract geometry (area → envelope, centroid → pan),
	// resample to 2^n uniform arc-length spacing for FFT.

	~processContour = { |contour|
		var segmentLengths, perimeter, fftSize, envX, envY;
		var xCoords = contour.flop[0];
		var yCoords = contour.flop[1];
		var numVertices = xCoords.size;
		var centroidX = xCoords.mean;
		var centroidY = yCoords.mean;
		var area = 0;

		// shoelace area
		numVertices.do { |i|
			var next = (i+1) % numVertices;
			area = area + (xCoords[i]*yCoords[next]) - (xCoords[next]*yCoords[i]);
		};
		area = (area / 2).abs;
		xCoords = xCoords - centroidX;
		yCoords = yCoords - centroidY;

		// arc-length resample via Env (piecewise-linear interpolator)
		segmentLengths = Array.fill(numVertices, { |i|
			var next = (i+1) % numVertices;
			hypot(xCoords[next]-xCoords[i], yCoords[next]-yCoords[i])
		});
		perimeter = segmentLengths.sum;
		fftSize = max(1024, numVertices).nextPowerOfTwo.clip(32, 1024);
		envX = Env(xCoords ++ [xCoords[0]], segmentLengths, \lin);
		envY = Env(yCoords ++ [yCoords[0]], segmentLengths, \lin);

		(
			contour: Array.fill(fftSize, { |i|
				var arcPos = i / fftSize * perimeter;
				[envX.at(arcPos), envY.at(arcPos)]
			}),
			fftSize: fftSize,
			pan: (centroidX / 50 - 1),
			area: (area / 100.squared)
		)
	};

	// ── 3. SPECTRAL ANALYSIS ────────────────────────────────────
	// Complex FFT: x = real, y = imaginary.
	// Positive bins [1..N/2-1] = CCW components, negative [N-1..N/2+1] = CW.
	// Interleaving offset: pos → integer ratios, neg → ratio + offset.

	~analyzeSpectrum = { |data|
		var samples = data[\contour];
		var fftSize = data[\fftSize];
		var spectrum, numBins, posBins, negBins, peakAmp, thresh, offset, posEnergy, negEnergy, posPeak, negPeak;

		spectrum = fft(
			Signal.newFrom(samples.flop[0]),
			Signal.newFrom(samples.flop[1]),
			Signal.fftCosTable(fftSize)
		);
		numBins = (fftSize / 2) - 1;

		// extract pos/neg bins: [harmonic, amplitude, phase]
		posBins = (1..numBins).collect { |bin|
			[bin,
				hypot(spectrum.real[bin], spectrum.imag[bin]),
				atan2(spectrum.imag[bin], spectrum.real[bin])]
		}.reject { |bin| bin[1] == 0 };

		negBins = (1..numBins).collect { |bin|
			var mirrorBin = fftSize - bin;
			[bin,
				hypot(spectrum.real[mirrorBin], spectrum.imag[mirrorBin]),
				atan2(spectrum.imag[mirrorBin], spectrum.real[mirrorBin])]
		}.reject { |bin| bin[1] == 0 };

		// threshold: discard bins below ~minDbAmp relative to peak
		peakAmp = (posBins ++ negBins).collect(_[1]).maxItem ? 0;
		thresh = peakAmp * ~minDbAmp.dbamp;
		posBins = posBins.select { |bin| bin[1] >= thresh };
		negBins = negBins.select { |bin| bin[1] >= thresh };

		// interleaving offset
		offset = case
		{ ~offsetMode == \fixed } { 0.5 }
		{ ~offsetMode == \phase } {
			// from phase difference of fundamental pos/neg pair
			var posF1 = Complex(spectrum.real[1], spectrum.imag[1]);
			var negF1 = Complex(spectrum.real[fftSize-1], spectrum.imag[fftSize-1]).conjugate;
			if((posF1.abs + negF1.abs) > 1e-6) {
				var rawOffset = (posF1.theta - negF1.theta).wrap(-pi, pi).linlin(-pi, pi, 0, 0.5);
				if(rawOffset < 0.05) { 0.5 } { rawOffset }
			} { 0.5 }
		}
		{ ~offsetMode == \symmetry } {
			// symmetric shapes → offset ~0.5 (harmonic), asymmetric → lower
			posEnergy = posBins.collect(_[1]).sum;
			negEnergy = negBins.collect(_[1]).sum;
			(min(posEnergy, negEnergy) / max(posEnergy, negEnergy).max(0.001))
			    .linlin(0, 1, 0.1, 0.5)
		};

		// normalize: \joint divides by global peak (preserves pos/neg ratio);
		// \independent divides each half by its own peak (equalises CW/CCW)
		if(~normMode == \independent) {
			posPeak = posBins.collect(_[1]).maxItem ? 1;
			negPeak = negBins.collect(_[1]).maxItem ? 1;
		} {
			posPeak = negPeak = (posBins ++ negBins).collect(_[1]).maxItem ? 1;
		};
		posBins = posBins.collect{ |bin| [bin[0], bin[1]/posPeak, bin[2]] };
		negBins = negBins.collect{ |bin| [bin[0], bin[1]/negPeak, bin[2]] };

		// peak-pick globally: keep the top ~maxPartials from the pooled spectrum
		// (tag bins with \pos/\neg so the offset can still be applied per half)
		{
			var pooled = (posBins.collect{ |b| b ++ [\pos] } ++ negBins.collect{ |b| b ++ [\neg] })
			    .sort{ |a, b| a[1] > b[1] }
			    .keep(~maxPartials);
			posBins = pooled.select{ |b| b[3] == \pos }.collect{ |b| b[0..2] };
			negBins = pooled.select{ |b| b[3] == \neg }.collect{ |b| b[0..2] };
		}.value;

		// interleave: pos → integer ratios, neg → ratio + offset
		data[\ratios] = posBins.collect{ |bin| bin[0].asFloat }
		             ++ negBins.collect{ |bin| bin[0].asFloat + offset };
		data[\amps]   = posBins.collect(_[1]) ++ negBins.collect(_[1]);
		data[\phases] = posBins.collect(_[2]) ++ negBins.collect(_[2]);
		data[\offset] = offset;

		// sort by ascending ratio
		{ var ord = data[\ratios].order;
			[\ratios, \amps, \phases].do { |k| data[k] = data[k][ord] };
		}.value;
		data
	};

	// ── 4. SYNTHESIS ────────────────────────────────────────────
	// Klang additive bank. Squeeze compresses ratios toward fundamental
	// to fit Nyquist while preserving relative spacing. Envelope from area.

	~play = { |data, fund=220|
		var ratios = data[\ratios].copy;
		var amps = data[\amps].copy;
		var nyquist = s.sampleRate * 0.45;
		var areaSize = data[\area].clip(0.005, 0.6);
		var maxRatio = ratios.maxItem ? 1;
		var squeeze;

		// area → envelope: small shapes = percussive, large = sustained
		var env = Env.perc(
			areaSize.linlin(0.005, 0.6, 0.001, 0.5),
			areaSize.linlin(0.005, 0.6, 0.08,  6.0),
			1.0,
			areaSize.linlin(0.005, 0.6, -8, -2)
		);

		// squeeze all ratios to fit under Nyquist
		squeeze = (((nyquist / fund) - 1) / (maxRatio - 1)).clip(0.01, 1.0);
		ratios = 1 + ((ratios - 1) * squeeze);
		amps = amps.collect{ |a, i| if((ratios[i] * fund) < nyquist) { a } { 0 } };

		if(~synth.notNil and: { ~synth.isPlaying }) { ~synth.free };
		~synth = { Pan2.ar(
			Klang.ar(`[ratios * fund, amps, data[\phases]])
			* EnvGen.kr(env, doneAction: 2) * 0.2,
			data[\pan]
		)}.play;

		"  >> % partials | %Hz | pan:% | area:% | offset:%".format(
			ratios.size, fund, data[\pan].round(0.01),
			data[\area].round(0.001), data[\offset].round(0.001)
		).postln;
	};

	// ── 5. DISPLAY ──────────────────────────────────────────────

	~drawShape = { |contour|
		var pts = contour.collect { |xy| (xy[0] / 100 * 280 + 10) @ (xy[1] / 100 * 280 + 10) };
		if(~shapeWin.isNil or: { ~shapeWin.isClosed }) {
			~shapeWin = Window("Shape", Rect(100, 100, 300, 300)).front;
		};
		~shapeWin.drawFunc = {
			Pen.strokeColor = Color.black;
			Pen.width = 1.5;
			Pen.moveTo(pts[0]);
			pts[1..].do { |pt| Pen.lineTo(pt) };
			Pen.lineTo(pts[0]);
			Pen.stroke;
		};
		~shapeWin.refresh;
	};

	// ── 6. LOAD AND PLAY ────────────────────────────────────────

	~loadAndPlay = { |filename, fund=220|
		var contour = ~parseSVG.(~svgDir +/+ filename);
		if(contour.notNil) {
			var data;
			"── % ──".format(filename).postln;
			~drawShape.(contour);
			data = ~processContour.(contour);
			~play.(~analyzeSpectrum.(data), fund);
		};
	};

	"Ready.".postln;
}; // waitForBoot
)

// ── USAGE (same freq to compare timbres) ────────────────────────
~loadAndPlay.("circle.svg", 220);
~loadAndPlay.("ellipse.svg", 220);
~loadAndPlay.("square.svg", 220);
~loadAndPlay.("triangle.svg", 220);
~loadAndPlay.("gear.svg", 220);
~loadAndPlay.("leaf.svg", 220);
~loadAndPlay.("blob.svg", 220);
~loadAndPlay.("star.svg", 220);
~loadAndPlay.("rosette.svg", 220);
~loadAndPlay.("pebble.svg", 220);

// ── OFFSET MODES ────────────────────────────────────────────────
~offsetMode = \symmetry;  // pos/neg energy ratio → harmonic ↔ inharmonic
~offsetMode = \phase;     // phase difference of fundamental bin pair
~offsetMode = \fixed;     // constant 0.5 offset

// ── NORMALIZATION MODES ─────────────────────────────────────────
~normMode = \joint;        // divide by global peak (preserves pos/neg ratio)
~normMode = \independent;  // divide each half by its own peak (equalises CW/CCW)
