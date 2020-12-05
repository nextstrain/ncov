localrules: download, upload, clean

ruleorder: align > download_aligned
ruleorder: filter > download_filtered
ruleorder: refilter > download_refiltered
ruleorder: mask > download_masked
ruleorder: diagnostic > download_diagnostic
