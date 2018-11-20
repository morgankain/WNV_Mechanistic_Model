PUSH = perl -wf $(filter %.pl, $^) $(filter-out %.pl, $^) > $@

ebird_stripped.tsv: select.csv ebird_data.tsv select.pl
	$(PUSH)
	