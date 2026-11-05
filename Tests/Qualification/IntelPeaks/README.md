# Intel bounded true-peak continuation — 5 November 2026

Base: published M18/main `51e47f300fcc38a386a5a92a7277029175cfab0c`.
Branch: `codex/intel-multicore-true-peak`, in the existing logic checkout.
This is a focused continuation, not the next LRA proposal milestone.

## Impact and required checks

The only product behaviour change removes the arm64 restriction from the
existing macOS measurement-service execution choice. Reuse M18's bounded
per-channel algorithm unchanged. It preserves the serial energy accumulation,
per-channel interpolation order, floating-point environment, earliest peak ties
and EOF drain. At most min(channels, six, logical CPUs minus one) peak lanes run.
Mono, 5.0, 7.1, fewer than three logical CPUs, non-concurrent peak observers,
small blocks and an unavailable floating-point environment retain serial work.
The build switch can disable bounded execution entirely.

Fresh checks on arm64 Sequoia and actual x86_64 High Sierra 10.13.6:

1. Direct core bit-exact serial/bounded parity at all four supported rates,
   channel orders, peak positions, padding, silence, block boundaries, observer
   fallback/errors and floating-point rounding modes. Explicit other-layout
   serial fallback assertions.
2. Existing short measurement and section-proposal smoke checks. No historical
   runner recursion, video load or physical device use.
3. Identical 60-second stereo and 5.1 PCM through the ordinary service objects
   and a private serial-only oracle with the same build identity. All detailed
   snapshot digests, section-peak maps, plans, explanations and seven proposal
   decisions must match exactly. Explicit Apply must invalidate Current; its
   independent measurement must pass and match the serial oracle exactly.
4. Alternating repeated service timings for Source, fresh section Propose,
   reused-attribution Propose and independent Current, with real timeline
   painting and a message-loop heartbeat. Report actual elapsed time and CPU
   utilisation, not extrapolation from the old mono three-hour fixture.
5. Native main-component cancellation during Source, section Propose and Current;
   no partial result adopted, no automatic Apply, restart succeeds, device stays
   closed. UI heartbeat target p95 <=120 ms with a 100 ms timer, maximum <=200 ms;
   cancellation <=1 s. Review a failure explicitly, never hide it by rerunning.
6. On Golem, inspect sampled stacks from the ordinary standalone application
   during eligible analysis to establish execution of dispatched peak workers.
   A forced diagnostic alone cannot satisfy this gate.

The fresh service equality checks include complete windows, numerical statuses,
gates, per-channel peaks and peak positions in snapshot digests. Pinned EBU/ITU
arithmetic/reference evidence from M18 is reused because filter coefficients,
arithmetic and decoding are unchanged; the direct and service exact serial
oracles cover the changed scheduling choice. M18's corrected Intel diagnostic
results are supporting evidence only, not ordinary-app qualification.

Budget: one CPU-heavy job per host, roughly 15–30 minutes per host after build,
under 2 GiB for package/fixtures/results and about 3 GiB for local compilation.
Use fresh owned directories; preserve every M18 package and earlier receipt.
No three-hour, cinema, codec matrix or physical listening rerun is justified by
this scheduling-only change. Expand only for a concrete failure or new impact.

Publication, tags and merge remain for review. Main and every existing tag stay
unchanged. Record candidate source/build/binary identities and the exact serial
oracle overlay with the results; the oracle must not masquerade as the app.
