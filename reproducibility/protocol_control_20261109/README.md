# Protocol-control implementation

This directory contains the post-fit scoring and fail-closed validation code used
for the registered control. The public core does not redistribute the registered
prediction fixture, participant-level prediction table or multiplicity file.
The complete evidence package supplies the fixture under its access terms. The
scoring specification is public. The checker validates required columns,
duplicate keys, participant labels, probabilities, run coverage, multiplicities
and incompatible declarations. The executable public worked example is under
`../examples/custom_panel/` and uses synthetic inputs.
