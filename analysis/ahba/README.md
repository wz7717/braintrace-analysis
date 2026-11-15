# AHBA frozen post-transform target

This directory supplies the Methods S32 post-transform target for the Editorial Concise manuscript. The SHA-256, axes, data scope and upstream terms appear in the binding. [Frozen target binding](FROZEN_TARGET_BINDING.json)

- Load query_matrix.npy as little-endian float32, with 231 rows and 5,000 columns. It is the exact query member of the original mixed replay NPZ. Rows follow the index in file order; columns follow the frozen floor+clip genes.json.
- The index has profile_id, donor_id, allowed_truth_json and prediction_failure. It verifies the 231 groups derived from the two official Allen SampleAnnot files (242 original assays). The original truth sets retain 223 eligible groups, 56 singleton and 167 allowed-set, and eight out-of-scope groups. No technical failures occurred.
- The full-library total belongs to the archived specimen reporting output. It is omitted from this minimal index. To regenerate the historical CSV bytes, obtain the exact full report inputs and original replay contract.
- The fixed 5,000-gene axis, state, classifier coefficients, intercepts and ordered classes can be taken from the identified public software release. The exact donor-equal Bo 10 by 5,000 float64 centroid member used by the Pearson comparator requires access through the project custodian under the applicable source and contribution conditions. Its identity is in the binding. An authorized reconstruction must match the frozen shape, axes and hash.
- This is a start from the transformed target. The full upstream pathway starts at the official Allen counts and needs the supplied transform, bound configuration, mappings, frozen state and additional authorized source inputs. The saved evidence describes the earlier local replay.

[Official data and source credit](../../sources/allen_ahba/SOURCE_CREDITS.csv) · [Allen Institute data-use and attribution conditions](https://alleninstitute.org/legal/terms-of-use) · [Historical local replay summary](../evidence/S32_AHBA_POSTTRANSFORM_CHECK.json)
