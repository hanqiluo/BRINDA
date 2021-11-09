## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Resubmission
This is a resubmission. In this version I have:
- Reduced the title to 65 characters
- Added description
- Streamlined the notification
- use exist("agp/crp", dataset) instead of !all(is.na(dataset$agp/crp)) to supress warning messages
