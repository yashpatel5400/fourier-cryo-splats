# Output-serialization repair and unchanged-design retry

30 September 2026. The v1 screen reached the first EMPIAR-10028 repeat, then
failed to serialize NumPy Boolean rejection flags. Its exception handler also
encountered the same flag. The initial record and generated partial arrays
remain unchanged in `mixture-validation-preflight-v1`; the full failure log
and their hashes are recorded separately. No completed scientific result was
printed. The partial arrays retain the first generated observation, latent
indices and all candidate means.

Convert rejection flags to native Python bool. Write the retry into a new
`mixture-validation-preflight-v2` directory; retain the exact same model,
seed, candidates, 16 repeats, solver settings and 20-minute budget from the
original protocol. This is an implementation repair, not new independent
confirmation: the first simulation has already been generated. Source and
retry note must be committed before execution. Preserve the failed attempt.
