# Heartbeats for Alakazam

A plugin for [Alakazam](https://github.com/markspan/Alakazam): every heartbeat as an event, for heartbeat-evoked potentials. The R-peaks of the ECG channel are found exactly as [spectHR](https://github.com/MarkSpan/spectHR) finds them, and each becomes an event with its interval to the next beat and its heart rate. DefineBins then chooses the heartbeats per condition, and Average gives one heartbeat-evoked potential per bin.

## Manual

The manual is on this repository's pages:
<https://markspan.github.io/AlakazamHREvoked/>.

## Installing

In Alakazam, **Install plugins** on the Alakazam tab lists this repository:
select it and press **Install**. It installs the latest release. A zip from
the [releases](https://github.com/markspan/AlakazamHREvoked/releases), or one built with
`build`, installs the same way through **From a link**, and so does the link
to this repository. The tool appears on the Tools tab, under Preprocessing.

It needs Alakazam V0.4.5.3 or newer, the first with **Install plugins**. It
needs no toolbox: spectHR's detection, scipy's filter and peak finder and
numpy's sums are ported in its `+Rpeak` package.

## Building and testing

```matlab
build()                                    % dist/Heartbeats-<version>.zip
runtests('tests')                          % the tool's tests and the install
```

The tests need Alakazam's source: the folder in `ALAKAZAM_ROOT`, or the
Alakazam repository beside this one. `SpectHRParityTest` holds the port to
spectHR's own results, which `tests/fixtures/makeSpectHRReference.py` made
once with spectHR's code; spectHR is not needed to run the tests.

## Licence

GPL-3.0, as Alakazam and spectHR.
