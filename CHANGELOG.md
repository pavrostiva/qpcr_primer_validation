# Changelog

The version is shown at the bottom of the sidebar and in the footer of every PDF page.

## 1.4.0 (2026-11-25)
- New buttons to drop the most dilute / the most concentrated dilution point, for one primer or for all primers at once, and to reset to the advisor's recommendation. At least 3 dilution points are always kept.
- The Smart LOQ advisor now explains itself when dropping the most dilute point(s) cannot fix a curve, instead of staying silent.

## 1.3.1 (2026-11-25)
- Fix: a layout cell that contains just `NTC` (or water / H2O / NFW) is now a negative control of the primer in the same row, not a primer named "NTC". The NTC check now works for such layouts.
- Summary table: status and efficiency come first; "< 3 points" is shown instead of "No Amp" when too few dilutions are left.
- Shorter Excel button label.

## 1.3.0 (2026-11-18)
- Show the app version in the sidebar and in the PDF footer.
- Well manager: tick or untick a whole dilution at once; wells that are not part of the curve are listed.
- Smart LOQ advisor can drop up to two of the most dilute points.

## 1.2.0
- PDF report with one page per primer; every primer is shown at once; the Excel report follows the ticked wells.

## 1.1.0
- Fixes: ZIP from macOS, primer names with underscores, primers without any Cq.
