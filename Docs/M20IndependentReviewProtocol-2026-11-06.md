# M20 independent engineering comparison protocol — 6 November 2026

The user supplied these three previously held-out files and explicitly selected the
controlled engineering comparison. Only container/channel metadata and content hashes
were inspected before this protocol was frozen; no holdout loudness/proposal results
informed the defaults. This is an engineering comparison, not a broadcaster delivery
specification or a listening-quality score.

Production is frozen at `f67e2d20bfbf7daa33ebf0152a675b9c8fa04595`, candidate 02.
The only change after candidate 01 is cached immutable native review pages. The 15
selected-control comparisons have passed on Sequoia and High Sierra. Their results
support retaining the approved defaults; independent holdouts cannot retrospectively
become tuning fixtures without explicit disclosure and another held-out review.

## Fixed conditions

Use each complete input at its original 48 kHz/24-bit PCM, with the role order explicitly
recorded in its iXML. Whole-programme applicability is an engineering-test assumption.
Initial Current is identical to Source. Evaluate M19 and explicitly selected M20
independently from that same baseline, with no accumulated processing between methods.

Requirements: IL -23 LUFS with the existing EBU profile's +/-0.2 LU range, true peak
<= -1 dBTP, no M/S limits, and LRA <= floor((baseline LRA - 1 LU) * 10) / 10 LU. This
uses the editable 0.1 LU step and asks for at least 1 LU reduction. If the baseline
has unavailable/unstable LRA or cannot support a positive reduced ceiling, disclose
that limitation instead of inventing a usable target. Full precision decides rules.

Keep policy at 3 dB additional / 12 dB cumulative attenuation, 20% newly affected
programme duration, no count restriction, and 0.5 dB reserve inside those caps.
M20 uses automatic roots, splitting enabled, two splits per seed, refinement enabled,
zero margins, 3 s minimum nominal interval, 1 s unchanged separation, gap fill off,
1 LU need threshold, 1 s persistence, global permission on, comparison intent off.
Keep the existing 64-proxy / three-complete-candidate limits and strict Apply gate.
No default, ranking, search or DSP change is authorized by the holdout result.

Record original content/order, exact requirements/preferences/settings, baseline
metrics, complete explanations, actual private attempts and forecast limitations,
retained plan/activity and observed resource use. A failed search is not proof of
impossibility. If a candidate passes, independently apply/remeasure only in the owned
test model, preserving original inputs. No session from the user's open app is edited.
No hardware opening, audible test, deliverable export or quality conclusion is part
of this numerical pass. Coordinate representative listening separately.

## Immutable input identities

- `/Users/fuzz/N15/Mixdown/DialogAmbiance.wav`
  SHA-256 `602485eef95774826a575b0069d358e6723a31f6928dbd6e04dd0df44c98eb90`; roles L, R.
- `/Users/fuzz/N15/Mixdown/MusicInBackAndForegroundSomeDialog.wav`
  SHA-256 `7ba6404041c94c909df290a6829a0fe92ce594329a1d575490c2bb4a5b6cc8f0`; roles L, R.
- `/Users/fuzz/N15/Mixdown/5.1_Busy.wav`
  SHA-256 `4fa30bfd9a8fbb64e68fc00eaf11f8adf81ead8502d05ce639b22cc07645d7b6`; roles L, R, C, LFE, Ls, Rs.
