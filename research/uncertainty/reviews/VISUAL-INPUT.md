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
rendered at an explicitly selected 125 or 150 dpi (150 by default). Each page, source file and complete multimodal payload has
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

At the later 41-page checkpoint, two local builds exceeded the unchanged
1,800,000-byte text guard (4,524,133 and 2,979,244 bytes). No model was invoked.
The next preview uses the documented read-only evidence mode: a 1,482,400-byte
initial text contains the manuscript, methods and summary-level evidence, while
exact per-case originals and archived sources are accessible through an indexed
copy. All 41 rendered pages are still included. The source hashes distinguish
the full originals from displayed eight-significant-digit numerical summaries;
the latter are not offered as numerical-error proofs. Columnar records preserve
all cases, including failures, and verified duplicate baseline rows are omitted
only because their original case records remain available. This is a successful
packet construction, not a new scientific review or acceptance assessment.

At the 53-page development checkpoint, an explicit 125-dpi preview preserved
all pages without cropping. Its text and multimodal inputs were 1,430,791 and
27,667,013 bytes, with 1,210 indexed evidence files. Pages 8, 32, 41 and 42
were visually inspected for legibility. This preview predates later manuscript
and completed-result changes and is not the final round-2 packet. Final input
sizing and page inspection must be repeated after the remaining outcomes.
