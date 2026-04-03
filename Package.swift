// swift-tools-version: 6.2
import PackageDescription

let package = Package(
    name: "PitchDSP",
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
    ]
)
