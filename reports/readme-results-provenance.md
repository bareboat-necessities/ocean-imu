# README result gallery provenance

The five shared-mode TFG tracking figures show the northbound 1.5 m PM–Stokes record from the v1.4.0 underway dataset. They were rendered with the standard replay plotter from the 200 Hz, default-seed `--regression` replay of `cc133f54807bf50673d0954db55c235a2e781b7e`. Its fixed-record gate passed. This replay uses the CI interval-rate profile. The ordinary gyro sampling configuration is scored separately in the degradation studies. The figure publisher refreshes these files from the following-sea `release-svgs-kalman_tfg_shared_modes` artifact.

The directional spectrum is a property of the input sea. Its existing SVG came from the `release-svgs-kalman_ou_iii` artifact of [build 44101355729](https://github.com/bareboat-necessities/ocean-imu/actions/runs/44101355729).
