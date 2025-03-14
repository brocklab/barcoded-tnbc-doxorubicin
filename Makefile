help: ## prints help for targets with comments
	@cat $(MAKEFILE_LIST) | grep -E '^[a-zA-Z_-]+:.*?## .*$$' | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'

bootstrap: ## generate pixi environment and run post-install
	pixi install
	pixi run post-install

data/stat.txt: .FORCE
	find . \
		\! -path ".git" \
		\! -path ".jj" \
		\! -path "*/.pixi/*" \
		\! -path "*/env/*" \
		\! -path "*cellranger*" \
		\! -path "./historical/*" \
		-a -path '*/data/*' \
		-o -path '*/outs/*' \
		-exec stat --terse {} \; > $@ \
	&& jj ci -m 'update stat' $@


DAG_TARGETS = \
	$(foreach num, 1 2 3 4 5 6 ,figure-$(num)) \
	$(foreach num, 1 2 3 4 5 6 7 8,supp-figure-$(num))

dags/figuregraph.svg: .FORCE
	snakemake --filegraph $(DAG_TARGETS) | sed "s|$(CURDIR)|\.|g" | dot -Tsvg > $@

dags/dag.svg: .FORCE
	snakemake -F --dag $(DAG_TARGETS) | dot -Tsvg > $@

.FORCE:
.PHONY: .FORCE bootstrap
