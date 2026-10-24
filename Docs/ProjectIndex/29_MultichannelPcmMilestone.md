# Subsequent multichannel PCM milestone — staged proposal

Recorded 23 October 2026 with M6 approval; **not implemented or qualified by M6**.
M6 retains 12-channel capacity and ordered height-role interfaces but restricts
Current State to the existing M1 measurement matrix. Unsupported layouts must fail
explicitly, without downmix, fallback to Source or expanded measurement claims.

Stage one should qualify eight-channel 7.1; stage two should qualify twelve-channel
7.1.4. Each is independently bounded and requires approval, source/binary manifests,
actual arm64 Sequoia and Intel High Sierra execution, and coordinated device tests.

Before implementation, specify exact speaker positions and accepted role aliases,
file/channel order and metadata interpretation. A channel count alone does not
identify a layout. Preserve original interleaved order; reject ambiguous or duplicate
roles and require explicit confirmation. Do not silently map side, rear, wide or
height positions onto existing 5.1 roles.

Required work for each stage:

- Review primary BS.1770 provisions for the applicable speaker positions and
  loudness weights, including height channels and LFE exclusion from loudness.
  Implement explicit mappings with independent positional/channel-isolation tests.
- Determine which versioned delivery profiles accept the layout. A measurement
  method supporting it does not establish broadcaster/profile applicability.
- Qualify all-channel sample/true peaks and locations, including LFE and every
  height channel; retain channel identity through reordered input and diagnostics.
- Extend supported measurement requests, per-channel details, M4/M5 diagnostics,
  measured Current evaluation and proposal boundaries deliberately. Preserve
  unavailable metrics and never treat an unqualified layout as compliant.
- Test all-channel processing and pre-monitor sample parity, independent routing
  to 8/12 outputs, mute/duplicate-output rejection and devices with insufficient
  outputs. No implicit downmix, fold-down or channel drop is introduced.
- Extend synthetic/reference fixtures at all qualified rates, permutation/identity
  tests, long-programme memory/performance, stale/cancel/device lifecycle and
  H.264/ProRes transport regressions on both actual Macs. Coordinate audible
  channel-identification and synchronization checks with the user.

Channel-based PCM 7.1/7.1.4 is distinct from Dolby Atmos objects, metadata, ADM/BWF
object import, bed/object rendering, binaural rendering or Atmos delivery approval.
Those need separate format, renderer, licensing and qualification decisions. No
Atmos capability is implied by twelve PCM channels or height-channel labels.

Interleaved rows should be labelled Channel 1, Channel 2, etc. A future separate-file
import milestone should show filenames alongside roles; separate-file import is
not included in M6 or automatically authorized by this follow-up proposal.
