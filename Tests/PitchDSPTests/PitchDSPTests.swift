// PitchDSPTests.swift — SPM smoke suite (runs under `swift test`, any
// standard Swift toolchain, no Makefile/cc required).
//
// This is deliberately a thin layer on top of the deep C suite in
// Tests/test_pitchdsp.c (run via `cd Tests && make run`). Each test below
// ports an assertion + tolerance directly from that suite; see the comment
// on each test for the exact C-suite counterpart. Do not loosen these
// tolerances without a corresponding change (and justification) in the C
// suite — they encode known-good behavior, not guesses.

import Testing
import PitchDSP

@Suite("PitchDSP smoke suite")
struct PitchDSPTests {

    // MARK: - Pure sine 440 Hz (two-stage refine accuracy)
    //
    // Ports Suite 5 `sweep_note("A4", 440.0, 48000.0)` → sweep_point's "+0c"
    // case in test_pitchdsp.c: median cents error over the second half of
    // frames must be <= 2.0. This is the two-stage YIN full-rate refine
    // pass being exercised (coarse decimated search + parabolic-interpolated
    // refine at native rate).
    @Test("pure 440 Hz sine refines to within 2 cents (median)")
    func pureSine440RefinesAccurately() {
        let sr: Float = 48000.0
        let n = 8192 * 20
        let samples = genSine(freq: 440.0, sr: sr, n: n)
        let frames = runDetector(samples: samples, sr: sr, config: pitchDetectorDefaultConfig())

        let secondHalf = frames[(frames.count / 2)...]
        var errs: [Double] = []
        for f in secondHalf where f.hz > 0 {
            errs.append(centsDiff(detected: f.hz, expected: 440.0))
        }
        #expect(!errs.isEmpty, "no valid detections in second half")
        errs.sort()
        let median = errs[errs.count / 2]
        #expect(abs(median) <= 2.0, "median cents error \(median) > 2.0")
    }

    // MARK: - A1 55 Hz sine (two-stage YIN's raison d'être: low-frequency accuracy)
    //
    // Ports Suite 2 `synthetic_test("A1_sine", 55.000f, 44100.0f)` in
    // test_pitchdsp.c: production default config, assert_corpus(0.05, 0.05).
    // A1 sits deep in the range the coarse decimated stage exists to serve
    // (5-string-bass-adjacent low end) — this is exactly the case the
    // two-stage design targets over a naive single-rate YIN.
    @Test("A1 55 Hz sine detects with corpus-grade accuracy")
    func a1SineDetectsAccurately() {
        let sr: Float = 44100.0
        let n = 8192 * 40
        let samples = genSine(freq: 55.000, sr: sr, n: n)
        let frames = runDetector(samples: samples, sr: sr, config: pitchDetectorDefaultConfig())
        let m = computeMetrics(frames, expectedHz: 55.000)

        #expect(m.firstFrame >= 0, "no detection within ±100 cents")
        if m.firstFrame >= 0 {
            #expect(m.firstFrame <= 25, "first detection frame \(m.firstFrame) > 25")
        }
        #expect(m.hasStable, "too few frames for stable window")
        guard m.hasStable else { return }
        #expect(m.detectionRate >= 0.80, "detectionRate \(m.detectionRate) < 0.80")
        #expect(m.wrongNoteRate <= 0.05, "wrongNoteRate \(m.wrongNoteRate) > 0.05")
        #expect(m.octaveErrorRate <= 0.05, "octaveErrorRate \(m.octaveErrorRate) > 0.05")
        #expect(m.stdDevCents <= 12.0, "stdDevCents \(m.stdDevCents) > 12.0")
    }

    // MARK: - Silence produces no valid pitch
    //
    // Ports `test_drain_includes_invalid_frames` in test_pitchdsp.c: silence
    // fails the RMS/peak gate, but every analysis cycle still queues a
    // result (invalid, hz = -1) — consumers rely on this for silence
    // debouncing per the header doc on pitchDetectorDrainResults.
    @Test("silence yields only invalid (hz < 0) results")
    func silenceYieldsNoValidPitch() {
        let sr: Float = 48000.0
        let cfg = pitchDetectorDefaultConfig()
        let detector = pitchDetectorCreate(8192, sr, cfg)
        #expect(detector != nil)
        guard let detector else { return }
        defer { pitchDetectorDestroy(detector) }

        let hop = 8192 / Int(cfg.hopDivisor)
        let silence = [Float](repeating: 0, count: hop * 4)
        silence.withUnsafeBufferPointer { buf in
            pitchDetectorProcess(detector, buf.baseAddress, Int32(silence.count))
        }

        var out = [PitchResult](repeating: PitchResult(), count: 16)
        let n = out.withUnsafeMutableBufferPointer { buf in
            pitchDetectorDrainResults(detector, buf.baseAddress, 16)
        }
        #expect(n == 4, "expected 4 invalid results, got \(n)")
        for i in 0..<Int(n) {
            #expect(out[i].hz < 0.0, "silence frame \(i) has hz=\(out[i].hz), want -1")
        }
    }

    // MARK: - maxHz cap raised to 1500 (app's tuner config): E6 detected, not subharmonic
    //
    // Ports the *intent* of `test_maxhz_floor_enforced` /
    // `test_default_config_has_maxhz` in test_pitchdsp.c (which prove maxHz
    // is enforced as a search-range floor) combined with the header's
    // documented behavior: "for a periodic signal whose fundamental is above
    // maxHz, the detector reports the first period multiple inside the
    // range (a subharmonic)". The app overrides the library default maxHz
    // (1000 Hz) to 1500 Hz specifically so E6 (~1318.5 Hz) is detected as
    // itself rather than folding to its subharmonic near 659 Hz (E5) — see
    // AVAudioEngineTunerListeningService.Constants.maxDetectableHz. This
    // test exercises that exact production config/frequency combination and
    // applies the same corpus tolerances as Suite 2 (assert_corpus 0.05/0.05).
    @Test("E6 sine (1318.5 Hz) with maxHz=1500 detects true fundamental, not the subharmonic")
    func maxHzRaisedDetectsE6NotSubharmonic() {
        let sr: Float = 44100.0
        let e6: Float = 1318.510
        var cfg = pitchDetectorDefaultConfig()
        cfg.maxHz = 1500.0

        let n = 8192 * 40
        let samples = genSine(freq: e6, sr: sr, n: n)
        let frames = runDetector(samples: samples, sr: sr, config: cfg)
        let m = computeMetrics(frames, expectedHz: e6)

        #expect(m.firstFrame >= 0, "no detection within ±100 cents of E6")
        #expect(m.hasStable, "too few frames for stable window")
        guard m.hasStable else { return }
        #expect(m.detectionRate >= 0.80, "detectionRate \(m.detectionRate) < 0.80")
        #expect(m.wrongNoteRate <= 0.05, "wrongNoteRate \(m.wrongNoteRate) > 0.05")
        #expect(m.octaveErrorRate <= 0.05, "octaveErrorRate \(m.octaveErrorRate) > 0.05 (subharmonic fold-back)")
        #expect(m.stdDevCents <= 12.0, "stdDevCents \(m.stdDevCents) > 12.0")

        // Direct, unambiguous check on the raw frames: must not collapse to
        // the ~659 Hz subharmonic that the default 1000 Hz maxHz would force.
        let validHz = frames.compactMap { $0.hz > 0 ? $0.hz : nil }
        #expect(!validHz.isEmpty)
        for hz in validHz {
            #expect(hz > 1000.0, "hz=\(hz) folded back below maxHz=1000 default cap — subharmonic lock")
        }
    }
}
