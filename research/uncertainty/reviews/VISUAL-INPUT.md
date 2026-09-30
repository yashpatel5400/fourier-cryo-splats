# Visual evidence in subsequent independent reviews

Round 1 received text source and result summaries but no rendered figures. Its
original review remains unmodified, including this limitation. A capability
preflight on the actual requested model, `claude-fable-5-1`, succeeded using
Claude CLI stream-json input with an image content block and stream-json output.
The authenticated provider confirmed the canonical model and first-party
provider. No tools, MCP servers or local customization were enabled.

The preflight asked only for a description of one project figure, not a review
or acceptance judgment. Its response correctly identified all three dataset
panels, axis labels, marker styles and curve behavior. The original response
event stream and manifest are preserved under `vision-preflight-02/`. Attempt 1
failed locally because stream-json input requires stream-json output; its error
is retained under `vision-preflight-01/`.

The review runner now sends the source packet followed by every manuscript page
rendered at 150 dpi. Each page, source file and complete multimodal payload has
a recorded SHA-256 digest. The exact input can be reconstituted from prompt.txt
and the saved page PNGs. Raw response events and the final result are preserved.
New theory, tables and completed or explicitly in-progress development records
are included. The next full review has not been invoked: visual capability is
not evidence that the scientific objections have been resolved.

A local-only full packet preflight subsequently passed with 31 rendered pages,
1,731,068 text bytes and 19,517,883 multimodal-input bytes. The first two local
packet builds exceeded a conservative text-context guard and did not invoke a
model. The passing projection retains every non-history outcome, includes
counts and endpoints of long iteration histories, and factors repeated snapshot
paths through their retained SHA-256 hashes. These omissions are explicit in
the reviewer instructions and each projected file's manifest. Intermediate
optimization traces remain available in their original public JSON files.
The preflight used the explicitly incomplete response plan only to test local
packet construction; it did not submit that plan as a completed rebuttal.
