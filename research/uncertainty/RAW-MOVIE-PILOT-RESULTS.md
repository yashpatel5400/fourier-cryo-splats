# Complete one-movie acquisition pilot

One preselected EMPIAR-10028 movie contains 16 frames of 4096 x 4096 float32 pixels. The first transfer broke after 922,746,880 saved bytes. A range request recovered the remaining 150,995,968 bytes after checking ETag, Last-Modified, Content-Range and size. The failed record, original partial file and recovery tail remain preserved. The assembled 1,073,742,848-byte movie has SHA-256 `6b31e50e80af4d95073e9eb5fa51d9f58de8353149935bfc0145a236dff973c5`.

All 64 fixed patches and eight disjoint frame pairs were analyzed: 512 patch differences and 1,792 within-patch cross-pair correlations. No patches or frames were selected by outcome.

| Statistic | Minimum | 5th percentile | Median | 95th percentile | Maximum |
| --- | ---: | ---: | ---: | ---: | ---: |
| sd | 341.721 | 351.211 | 357.687 | 370.293 | 402.494 |
| horizontal correlation | 0.645475 | 0.649616 | 0.656851 | 0.664067 | 0.684395 |
| vertical correlation | 0.646312 | 0.652151 | 0.658274 | 0.666227 | 0.686223 |
| disjoint difference correlation | -0.0185616 | -0.00779825 | 0.000657574 | 0.00866096 | 0.0163858 |

Frame means range from 2052.533 to 2060.282; frame SDs range from 365.189 to 366.351, in deposited pixel units. Maximum absolute Parseval discrepancy is 5.82e-11.

Spatial correlation is substantial. Small correlations between disjoint differences are descriptive and do not prove independence, Gaussianity, or covariance transfer. These differences retain dose, motion and detector contributions. The single movie does not validate a homogeneous density, calibrate poses, or resolve the source-group identity limitation. No density intervals or reconstruction are derived from this pilot.

The original movie remains at its recorded EMPIAR URL rather than being duplicated in the release. The release includes every diagnostic, display, protocol, transfer record and analysis source snapshot. The original failed prefix and tail are transport artifacts, not extra scientific observations.
