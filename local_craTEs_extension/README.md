# Local craTEs Extension

Source these files after `library(craTEs)` to override selected package
functions without reinstalling the package.

- `getActivities.R`: keeps the package API, adds `lambda_max_by_sample`, and
  returns `activities` as a base matrix.
- `getSignifActivities.R`: keeps the package API for paired groups. In
  `contrast_mode = "auto"`, equal-sized groups use the original paired
  implementation; unequal-sized groups use an explicit group-mean contrast
  instead of R's accidental vector recycling.
- `preprocess_E_N_for_activities_hg38.R`: hg38-compatible preprocessing helper.
