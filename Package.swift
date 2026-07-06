// swift-tools-version: 6.2
import PackageDescription

let package = Package(
    name: "PitchDSP",
    platforms: [
        .macOS(.v13)   // for the mictuner dev tool; library itself is platform-agnostic C
    ],
    products: [
        .library(
            name: "PitchDSP",
            targets: ["PitchDSP"]
        ),
    ],
    targets: [
        .target(
            name: "PitchDSP",
            publicHeadersPath: "include"
        ),
        // Dev tool: live pitch detection from the Mac microphone.
        //   swift run mictuner [--yin 0.20] [--fallback 0.30] [--minconf 0.60]
        //                      [--minhz 25] [--maxhz 1000] [--window 8192] [--verbose]
        .executableTarget(
            name: "mictuner",
            dependencies: ["PitchDSP"]
        ),
    ]
)
