# Heartbeats for Alakazam

A plugin for [Alakazam](https://github.com/markspan/Alakazam): every heartbeat as an event, for heartbeat-evoked potentials. The R-peaks of an ECG are found exactly as [spectHR](https://github.com/MarkSpan/spectHR) finds them, and each becomes an event with its interval to the next beat and its heart rate. The ECG can be a channel of the recording, a stream of the XDF file it was read from (Lab Streaming Layer), or a recording by another device that received the same triggers, its clock lined up with the EEG's on them; a check on the heart's own field in the EEG confirms that the beats are on time. DefineBins then chooses the heartbeats per condition, and Average gives one heartbeat-evoked potential per bin.

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
numpy's sums are ported in its `+Rpeak` package, and the synchronisation is
in its `+Sync` package. An ECG stream of an XDF file is read with EEGLAB's
xdfimport, which Alakazam installs to read XDF files.

## Building and testing

```matlab
build()                                    % dist/Heartbeats-<version>.zip
runtests('tests')                          % the tool's tests and the install
```

The tests need Alakazam's source: the folder in `ALAKAZAM_ROOT`, or the
Alakazam repository beside this one. `SpectHRParityTest` holds the port to
spectHR's own results, which `tests/fixtures/makeSpectHRReference.py` made
once with spectHR's code; spectHR is not needed to run the tests.
`SyncTest` makes up recordings with known heartbeats and clocks (it writes
its own XDF files) and needs EEGLAB with xdfimport; without them its XDF
and second-recording cases are skipped.

## Licence

GPL-3.0, as Alakazam and spectHR.
