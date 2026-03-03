mkdir -p rev/lists

# memfix: mem_fix OR rerun_missing OR inspect_parse OR parse_error (if any)
xsv search -s action mem_fix meta/revision_plan.tsv | xsv select species | tail -n +2 > rev/lists/to_resubmit_memfix.txt
xsv search -s action rerun_missing meta/revision_plan.tsv | xsv select species | tail -n +2 >> rev/lists/to_resubmit_memfix.txt
xsv search -s action inspect_parse meta/revision_plan.tsv | xsv select species | tail -n +2 >> rev/lists/to_resubmit_memfix.txt
xsv search -s action inspect_parse meta/revision_plan.tsv | xsv select species | tail -n +2 >> rev/lists/to_resubmit_memfix.txt

sort -u rev/lists/to_resubmit_memfix.txt -o rev/lists/to_resubmit_memfix.txt

# deadend: dead_end_gt0==1 AND not in memfix list
xsv search -s dead_end_gt0 1 meta/revision_plan.tsv \
| xsv select species \
| tail -n +2 \
| sort -u \
| comm -23 - rev/lists/to_resubmit_memfix.txt > rev/lists/to_resubmit_deadend.txt

wc -l rev/lists/to_resubmit_memfix.txt rev/lists/to_resubmit_deadend.txt
