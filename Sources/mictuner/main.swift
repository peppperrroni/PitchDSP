// mictuner — live PitchDSP testing from the Mac microphone.
//
//   swift run mictuner [--yin 0.20] [--fallback 0.30] [--minconf 0.60]
//                      [--minhz 25] [--maxhz 1000] [--window 8192] [--verbose]
//
// Captures the default input device via AVAudioEngine (the same tap mechanism
// the iOS app uses), feeds every buffer to the detector, drains every analysis
// result, and prints one line per frame. Invalid frames print only with
// --verbose. Ctrl-C to stop.
//
// First run: macOS will ask for microphone permission for your terminal app.

#if os(macOS)
import AVFoundation
import PitchDSP

// MARK: - Args

var cfg = pitchDetectorDefaultConfig()
var windowSize: Int32 = 8192
var verbose = false

var it = CommandLine.arguments.dropFirst().makeIterator()
@MainActor
func floatArg(_ flag: String) -> Float {
    guard let v = it.next(), let f = Float(v) else {
        FileHandle.standardError.write("missing/invalid value for \(flag)\n".data(using: .utf8)!)
        exit(2)
    }
    return f
}
while let a = it.next() {
    switch a {
    case "--yin":      cfg.yinThreshold      = floatArg(a)
    case "--fallback": cfg.fallbackThreshold = floatArg(a)
    case "--minconf":  cfg.minConfidence     = floatArg(a)
    case "--minhz":    cfg.minHz             = floatArg(a)
    case "--maxhz":    cfg.maxHz             = floatArg(a)
    case "--window":   windowSize            = Int32(floatArg(a))
    case "--verbose":  verbose = true
    case "-h", "--help":
        print("usage: mictuner [--yin F] [--fallback F] [--minconf F] [--minhz F] [--maxhz F] [--window N] [--verbose]")
        exit(0)
    default:
        FileHandle.standardError.write("unknown argument: \(a)\n".data(using: .utf8)!)
        exit(2)
    }
}

// MARK: - Note formatting

let noteNames = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
func describe(_ hz: Float) -> String {
    let midiExact = 69.0 + 12.0 * log2(Double(hz) / 440.0)
    let midi = Int(midiExact.rounded())
    guard midi >= 0, midi < 128 else { return "?" }
    let cents = Int(((midiExact - Double(midi)) * 100.0).rounded())
    let name = "\(noteNames[midi % 12])\(midi / 12 - 1)"
    return String(format: "%-4@ %+4d¢", name as NSString, cents)
}

// MARK: - Capture

let engine = AVAudioEngine()
let input = engine.inputNode
let format = input.outputFormat(forBus: 0)
let sampleRate = Float(format.sampleRate)

guard let detector = pitchDetectorCreate(windowSize, sampleRate, cfg) else {
    FileHandle.standardError.write("pitchDetectorCreate failed (window must be a multiple of 4)\n".data(using: .utf8)!)
    exit(1)
}

print(String(format: "mictuner: %.0f Hz input, window %d (%.0f ms), hop %d (~%.0f results/s)",
             sampleRate, windowSize, Float(windowSize) / sampleRate * 1000,
             windowSize / cfg.hopDivisor,
             sampleRate / Float(windowSize / cfg.hopDivisor)))
print(String(format: "config: yin %.2f  fallback %.2f  minconf %.2f  range %.0f–%.0f Hz",
             cfg.yinThreshold, cfg.fallbackThreshold, cfg.minConfidence, cfg.minHz, cfg.maxHz))
print("time      hz        note        conf")

// The tap closure runs on AVFAudio's own queue. It must be @Sendable and
// nonisolated: capturing any top-level (MainActor-isolated) var would give the
// closure dynamic MainActor isolation and crash with dispatch_assert_queue
// when the audio queue invokes it. Capture immutable copies instead; the
// detector and scratch buffer are single-thread by API contract (tap thread).
let startTime = Date()
let showInvalid = verbose
nonisolated(unsafe) let det = detector
nonisolated(unsafe) var results = [PitchResult](repeating: PitchResult(), count: 16)

input.installTap(onBus: 0, bufferSize: 1024, format: format) { @Sendable buffer, _ in
    guard let channel = buffer.floatChannelData?[0] else { return }
    pitchDetectorProcess(det, channel, Int32(buffer.frameLength))
    let n = Int(results.withUnsafeMutableBufferPointer {
        pitchDetectorDrainResults(det, $0.baseAddress, Int32($0.count))
    })
    let t = Date().timeIntervalSince(startTime)
    for i in 0..<n {
        let r = results[i]
        if r.hz > 0 {
            print(String(format: "%7.2fs %8.2f Hz  %@  %.3f  rms %.4f", t, r.hz, describe(r.hz), r.confidence, r.rms))
        } else if showInvalid {
            print(String(format: "%7.2fs        —                      rms %.4f", t, r.rms))
        }
    }
}

do {
    try engine.start()
} catch {
    FileHandle.standardError.write("engine start failed: \(error)\n".data(using: .utf8)!)
    exit(1)
}

signal(SIGINT) { _ in
    print("\nbye")
    exit(0)
}
RunLoop.main.run()

#else
fatalError("mictuner is a macOS-only dev tool")
#endif
