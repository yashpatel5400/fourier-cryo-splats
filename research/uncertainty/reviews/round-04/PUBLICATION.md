# Review publication

The reviewer's final report and response metadata are unchanged. The raw local
provider stream is retained under its recorded hash. `public-events.jsonl`
preserves event order and tool records, replacing internal reasoning blocks
with a redaction marker and block hash. It is not the raw byte stream. This
redaction does not remove or change the public verdict, reviewer findings,
assistant text messages or tool results.

Redacted internal-reasoning blocks: 40.
