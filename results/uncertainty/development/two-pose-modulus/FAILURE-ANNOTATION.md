# Postmortem metadata disclosure

The initial run completed its first numerical fit but crashed while serializing
a NumPy boolean. Its exception handler hit the same serialization problem.
The on-disk case JSON therefore still held the initial incomplete record.
After the process exited, the agent added `error`, `failure_log` and
`failure_log_sha256` to that record to document the actual failure. Those fields
were not successfully written by the crashed script itself. No completed fit,
convergence or successful scientific outcome is claimed for this attempt.
The corrected run uses a distinct directory and preserves this record and its
original source snapshot.
