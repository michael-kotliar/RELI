## RELI (refactored)

Params:
```
    --snp            SNP file
    --ld             LD file
    --target         Target interval file
    --chr            Chromosome length file. Should be sorted lexicographicaly (1, 10, 11 ... 19, 2, 20, 21 ... X, Y)
    --null           Null model file
    --dbsnp          SNP database file. Shouldn't have any header line
    --out            Output directory, default "./results"
    --prefix         Output file prefix, default "reli"
    --rep            Permutation number, default 2000
    --corr           Correction multiplier, default 1
```


Does the chromosome order in `--chr` file influence the results?

`--chr`, `--dbsnp`, ` --null`, `--target`, `--snp`, `--ld` are not supposed to have any header 

from `--ld` all `:` should be removed