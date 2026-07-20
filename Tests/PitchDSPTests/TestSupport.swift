// TestSupport.swift — shared signal generators and detector-running harness
// for the SPM smoke-test target.
//
// This is a Swift port of the harness in Tests/test_pitchdsp.c (the deep C
// suite run via `cd Tests && make run`). Tolerances asserted by the tests in
// this target are taken directly from that suite — see each test's comment
// for the exact C-suite test it mirrors. This target is the "standard
// toolchain" smoke layer; the C suite (WAV corpus, harmonic corpus, cents
// sweep, white-box internals) remains the deep suite of record.

import PitchDSP
#if canImport(Darwin)
import Darwin
#elseif canImport(Glibc)
import Glibc
#endif

/// One drained analysis result, trimmed to what the tests need.
struct Frame {
    let hz: Float
    let confidence: Float
}

/// Pure sine generator — mirrors `gen_sine` in test_pitchdsp.c.
func genSine(freq: Float, sr: Float, n: Int) -> [Float] {
    (0..<n).map { i in
        sinf(2.0 * Float.pi * freq * Float(i) / sr)
    }
}

/// Runs the detector over `samples` in 512-sample chunks, draining after each
/// chunk — mirrors `run_detector` in test_pitchdsp.c exactly (chunk size,
/// drain cap of 16, windowSize 8192).
func runDetector(samples: [Float], sr: Float, config: PitchDetectorConfig) -> [Frame] {
    guard let detector = pitchDetectorCreate(8192, sr, config) else {
        fatalError("pitchDetectorCreate returned NULL")
    }
    defer { pitchDetectorDestroy(detector) }

    var frames: [Frame] = []
    let chunk = 512
    var i = 0
    while i < samples.count {
        let sz = min(chunk, samples.count - i)
        samples.withUnsafeBufferPointer { buf in
            pitchDetectorProcess(detector, buf.baseAddress! + i, Int32(sz))
        }
        var drained = [PitchResult](repeating: PitchResult(), count: 16)
        let nd = drained.withUnsafeMutableBufferPointer { buf in
            pitchDetectorDrainResults(detector, buf.baseAddress, 16)
        }
        for k in 0..<Int(nd) {
            frames.append(Frame(hz: drained[k].hz, confidence: drained[k].confidence))
        }
        i += sz
    }
    return frames
}

// MARK: - Metrics (port of compute_metrics / assert_corpus in test_pitchdsp.c)

func hzToMidi(_ hz: Float) -> Int {
    guard hz > 0 else { return -1 }
    return Int(69.0 + 12.0 * log2(Double(hz) / 440.0) + 0.5)
}

func centsDiff(detected: Float, expected: Float) -> Double {
    guard detected > 0, expected > 0 else { return 1e9 }
    return 1200.0 * log2(Double(detected) / Double(expected))
}

struct Metrics {
    var firstFrame = -1
    var hasStable = false
    var detectionRate = 0.0
    var wrongNoteRate = 0.0
    var octaveErrorRate = 0.0
    var stdDevCents = 0.0
}

/// Direct port of compute_metrics(). Stable window is frames 14..49, matching
/// the C suite exactly (same windowSize/hopDivisor/chunking → same frame
/// cadence for the same input length).
func computeMetrics(_ frames: [Frame], expectedHz: Float) -> Metrics {
    var m = Metrics()
    let expectedMidi = hzToMidi(expectedHz)

    for (idx, f) in frames.enumerated() {
        if f.hz > 0, abs(centsDiff(detected: f.hz, expected: expectedHz)) <= 100.0 {
            m.firstFrame = idx
            break
        }
    }

    let start = 14, end = 49
    guard frames.count > start else { return m }
    m.hasStable = true
    let last = min(frames.count - 1, end)
    let total = last - start + 1
    var valid = 0, wrong = 0, octave = 0
    var cerr: [Double] = []

    for idx in start...last {
        let hz = frames[idx].hz
        guard hz > 0 else { continue }
        valid += 1
        let dMidi = hzToMidi(hz)
        let c = centsDiff(detected: hz, expected: expectedHz)
        if dMidi != expectedMidi {
            if abs(dMidi - expectedMidi) % 12 == 0 {
                octave += 1
            } else {
                wrong += 1
            }
        } else if abs(c) <= 8.0 {
            cerr.append(abs(c))
        }
    }

    m.detectionRate = Double(valid) / Double(total)
    m.wrongNoteRate = Double(wrong) / Double(total)
    m.octaveErrorRate = Double(octave) / Double(total)

    if !cerr.isEmpty {
        let mean = cerr.reduce(0, +) / Double(cerr.count)
        let variance = cerr.reduce(0.0) { $0 + ($1 - mean) * ($1 - mean) } / Double(cerr.count)
        m.stdDevCents = variance.squareRoot()
    }
    return m
}
