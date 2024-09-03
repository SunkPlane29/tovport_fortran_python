.PHONY: enter-repl
enter-repl: init-out
	julia --threads=auto --project=. -i -L "main.jl"

.PHONY: init-out
init-out:
	mkdir -p out

.PHONY: run
run: init-out
	time -p julia --threads=auto --project=. -L "main.jl" -e "main()"
