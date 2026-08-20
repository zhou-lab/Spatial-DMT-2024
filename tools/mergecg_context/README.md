# mergecg_context

`mergecg_context` merges adjacent C- and G-strand CpG calls from the
nine-column output of `biscuit vcf2bed -e -t c`. It reproduces the
single-sample text emitted by `biscuit mergecg`, but uses the reference base
already present in column four instead of loading each reference chromosome.
Its merge and rounding behavior follows BISCUIT's MIT-licensed `mergecg`.

Input must be coordinate sorted. Non-CpG rows are ignored. Output is written
to stdout.
