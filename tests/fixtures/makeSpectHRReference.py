"""Reference results of spectHR's own R-peak detection, for SpectHRParityTest.

Heartbeats detects R-peaks exactly as spectHR does, ported to MATLAB: the ECG
polarity check (spectHR.signal.ecg.detect_ecg_polarity), the 0.5 Hz high-pass
that precedes detection (spectHR.signal.filters.butterworth_filter, order 4,
as spectHR.session.Samples.filtered calls it) and the detector itself
(spectHR.signal.rpeak.detect_rpeaks). The plugin does not use spectHR. This
script runs spectHR's own code once, on five minutes of the ECG of spectHR's
VU-AMS example recording and on two made-up signals that reach the parts real
ECG seldom does (the polarity tiebreaker, and scipy's peak finder on plateaus
and close peaks), and writes every input and result to spectHRReference.mat
beside it, which the test holds the port to.

Run it with spectHR's own Python, giving spectHR's repository:

    <spectHR>/.venv/Scripts/python.exe makeSpectHRReference.py <spectHR>
"""
import pathlib
import sys

import numpy as np
import scipy.signal as signal
from scipy.io import savemat

SPECTHR = pathlib.Path(sys.argv[1]).resolve()
sys.path.insert(0, str(SPECTHR / "src"))

from spectHR.dataset import load                                                 # noqa: E402
from spectHR.dataset.preprocessing import apply_canonical_channels, resolve_ecg  # noqa: E402
from spectHR.signal.ecg import detect_ecg_polarity                               # noqa: E402
from spectHR.session import Samples                                              # noqa: E402
from spectHR.signal.rpeak import detect_rpeaks                                   # noqa: E402

out = {}

# ---- five minutes of a real ECG ------------------------------------------
session = apply_canonical_channels(load(SPECTHR / "ExampleData" / "data" / "VU-AMS_5fs_example_data.edf"))
ecg = resolve_ecg(session)
srate = float(ecg.srate)
values = np.asarray(ecg.values, dtype=float)
middle = values.size // 2
half = int(150 * srate)
chunk = values[middle - half:middle + half]
times = np.arange(chunk.size) / srate

# spectHR's own route (apply_ecg_polarity, then apply_beat_detection with its
# default prefilter): the polarity on the whole signal, the upright signal
# high-passed by Samples.filtered (whose rate is 1 / median time step), then
# detect_rpeaks (whose rate is 1 / mean time step).
polarity, debug = detect_ecg_polarity(times, chunk, return_debug=True)
upright = Samples(times, chunk if polarity == "normal" else -chunk, "ecg")
highpassed = upright.filtered(filter_type="highpass", cutoff=0.5)
peaks = detect_rpeaks(highpassed, min_peak_distance_ms=300.0)
inverted, inverted_debug = detect_ecg_polarity(times, -chunk, return_debug=True)

# The times are not stored: (0:n-1) / srate rebuilds them exactly. Of the
# filtered signal, its first, middle and last 5000 samples, where the edges
# and the interior of the filter show.
filtered = np.asarray(highpassed.values)
n = filtered.size
keep = np.r_[0:5000, n // 2 - 2500:n // 2 + 2500, n - 5000:n]
out.update(
    ecg=chunk, srate=srate,
    polarity=polarity, skewness=debug["skewness"],
    highpassedAt=keep + 1, highpassed=filtered[keep], peakTimes=peaks,
    minDistanceSamples=max(1, int((300.0 / 1000.0) * (1.0 / float(np.mean(np.diff(times)))))),
    invertedPolarity=inverted, invertedSkewness=inverted_debug["skewness"],
    # numpy's own sums, which the port must give to the last bit
    numpySum=float(np.sum(chunk)), numpyMeanStep=float(np.mean(np.diff(times))),
    numpyStd=float(np.std(chunk)), numpyMedian=float(np.median(chunk)),
    # the detector alone, on the unfiltered ECG: nothing but numpy's sums and
    # plain arithmetic, so the port must match it to the last bit
    rawPeakTimes=detect_rpeaks(Samples(times, chunk, "ecg"), min_peak_distance_ms=300.0),
    rawThreshold=float(np.median(chunk) + 1.5 * np.std(chunk)),
)

# The same at exactly 500 Hz, as an EEG recording stores its rate: the first
# two minutes, every second sample (ecg(1:2:120000) in MATLAB).
eeg_rate = 500.0
at500 = chunk[: int(120 * srate)][::2]
times500 = np.arange(at500.size) / eeg_rate
polarity500 = detect_ecg_polarity(times500, at500)
upright500 = Samples(times500, at500 if polarity500 == "normal" else -at500, "ecg")
peaks500 = detect_rpeaks(upright500.filtered(filter_type="highpass", cutoff=0.5), min_peak_distance_ms=300.0)
out.update(peakTimes500=peaks500, polarity500=polarity500)

# ---- the polarity tiebreaker ---------------------------------------------
# Spikes up and spikes down on a slow wave, the downward ones taller but
# narrower: searched for the first shape whose skewness is too small to
# decide (|skew| < 0.1), so the peak prominences do. Once with the downward
# spikes dominant, once with the upward ones.
tie_srate = 250.0
t = np.arange(int(40 * tie_srate)) / tie_srate


def spikes(up, down, up_width, down_width):
    signal_ = 0.2 * np.sin(2 * np.pi * 0.3 * t)
    for k, centre in enumerate(np.arange(0.5, 39.5, 0.8)):
        size, width = (up, up_width) if k % 2 == 0 else (-down, down_width)
        signal_ = signal_ + size * np.exp(-0.5 * ((t - centre) / width) ** 2)
    return signal_


for name, flip in (("tie", 1.0), ("tieUp", -1.0)):
    found = None
    for ratio in np.linspace(0.5, 1.0, 501):
        candidate = flip * spikes(1.0, 1.0 / ratio ** 0.5, 0.016, 0.016 * ratio)
        polarity_, debug_ = detect_ecg_polarity(t, candidate, return_debug=True)
        if debug_["decision_source"] == "peak_count":
            found = (candidate, polarity_, debug_)
            break
    if found is None:
        raise SystemExit(f"no signal found that reaches the tiebreaker ({name})")
    candidate, polarity_, debug_ = found
    out.update({name: candidate, name + "Polarity": polarity_, name + "Skewness": debug_["skewness"],
                name + "Score": debug_["peak_score"]})
out["tieTimes"] = t

# ---- scipy's peak finder on what ECG seldom has --------------------------
# Plateaus (flat tops of two and three samples), peaks closer than the
# distance, and two peaks of equal height farther apart than it (numpy's sort
# leaves the order of equal heights unspecified, so equal heights within the
# distance would make the reference itself arbitrary).
rng = np.random.default_rng(3)
raw = rng.normal(size=600)
raw[100:103] = 3.0
raw[200:202] = 2.5
raw[300] = 3.2
raw[320] = 3.2
peaks_hd, _ = signal.find_peaks(raw, height=0.5, distance=7)
peaks_dp, props = signal.find_peaks(raw, distance=5, prominence=1.0)
out.update(findPeaksSignal=raw, findPeaksHeightDistance=peaks_hd + 1,
           findPeaksDistanceProminence=peaks_dp + 1, prominences=props["prominences"])

target = pathlib.Path(__file__).with_name("spectHRReference.mat")
savemat(target, out, do_compression=True)
print(f"wrote {target}: {chunk.size} samples at {srate:g} Hz, {peaks.size} beats, polarity {polarity}, "
      f"min distance {out['minDistanceSamples']} samples, tiebreaks -> {out['tiePolarity']} "
      f"(skew {out['tieSkewness']:.3f}) and {out['tieUpPolarity']} (skew {out['tieUpSkewness']:.3f})")
