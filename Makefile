.DEFAULT_GOAL := help
.PHONY: help papers paper1 paper2 paper3 checks check-full verify verify-full

help:
	@printf '%s\n' \
	  'make papers       Build all three manuscripts' \
	  'make paper1       Build the general-constraints paper' \
	  'make paper2       Build the quantum-transition paper' \
	  'make paper3       Build the magnetic-sheaf technical note' \
	  'make checks       Run ten Python checks' \
	  'make check-full   Run all 38 Python/Sage scripts' \
	  'make verify       Run Python checks and build all papers' \
	  'make verify-full  Run all checks and build all papers'

papers:
	@./scripts/build.sh

paper1 paper2 paper3:
	@./scripts/build.sh $(patsubst paper%,%,$@)

checks:
	@./scripts/check.sh

check-full:
	@./scripts/check.sh --full

verify: checks papers

verify-full: check-full papers
