// swift-tools-version: 6.2
import PackageDescription

let package = Package(
    name: "PitchDSP",
    // The C core has no OS API dependency (libm only), so these floors are
    // driven by the consumer, not the algorithm: Nocturne pins
    // IPHONEOS_DEPLOYMENT_TARGET = 17.0, and macOS 13 is a conservative
    // modern floor for the standalone `swift test` smoke suite (host-only —
    // no runtime dependency on the app's macOS 15.6 target).
    platforms: [
        .iOS(.v17),
        .macOS(.v13),
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
        .testTarget(
            name: "PitchDSPTests",
            dependencies: ["PitchDSP"]
        ),
    ]
)
