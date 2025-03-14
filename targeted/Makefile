PYCASHIER = \
	docker run \
	--rm -it \
	-v "$(shell pwd)":/data \
	-u "$(shell id -u):$(shell id -g)" \
	ghcr.io/brocklab/pycashier:v2024.1004

merge:
	$(PYCASHIER) merge

extract:
	$(PYCASHIER) extract

receipt:
	$(PYCASHIER) receipt -i data/outs -o data/combined.tsv

extract-fc0:
	$(PYCASHIER) extract -fc 0 -o data/outs-fc0

receipt-fc0:
	$(PYCASHIER) receipt -i data/outs-fc0 -o data/combined-fc0.tsv

raw:
	./scripts/fetch-raw.py

.PHONY: merge extract raw env
