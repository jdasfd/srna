# Classified bacteria in their habitat preference with land plants

I classified bacteria into four different types according to bacteria-plant(host) relations. Four categories are: endo/epiphyte, environment, gut and marine.

All the habitat information of bacteria selected were acquired from different references. Please go and check `ASSEMBLY.xlsx` for more details.

## Extract bacteria name with its accession numbers

In RefSeq genome, there were different assembly species. In order to match the species to its correspoding accession number, a few step should be finished.

```bash
```

## Split bacterial species to different group

After we specified 161 (only include category 1-4) species habitat (details in `ASSEMBLY.xlsx`), we grouped them according to the bacteria spatial distance with the landplants. I numbered every group preventing too much characters in a tsv file. After that I saved them into `name.tsv`.

|     group     | number | species number |
| :-----------: | :----: | :------------: |
| endo/epiphyte |   1    |       38       |
|  environment  |   2    |       50       |
|      gut      |   3    |       68       |
|    marine     |   4    |       5        |
