# RIsearch2

RIsearch 1 and RIresearch 2 were all collected to one package RIsearch-1.2.tar.gz.

```bash
cd ~/share
wget -P /tmp https://rth.dk/resources/risearch/download/RIsearch-1.2.tar.gz
# combined RIsearch1/2 pacakge
tar -zxvf /tmp/RIsearch-1.2.tar.gz
cd RIsearch-1.2
cd RIsearch2
./rebuild.sh
echo '# Risearch2' >> ~/.bashrc
echo 'export PATH="/home/j/share/RIsearch-1.2/RIsearch2/bin:$PATH"' >> ~/.bashrc
```

```bash
risearch2.x -h

================================ RIsearch v2.1 ===============================
================ Energy based RNA-RNA interaction predictions ================

Usage: risearch2.x [options]

  -h,         --help
                 show this message
------------------------------------------------------------------------------
--------------------------- SUFFIX ARRAY CREATION ----------------------------
  -c <FILE>,  --create=FILE (.fa or .fa.gz)
                 create suffix array for target sequence(s) together with
                 their reverse complements, FASTA format, use '-' for stdin
  -o <FILE>,  --output=FILE
                 save created suffix array to given index file path
------------------------------------------------------------------------------
--------------------------- INTERACTION PREDICTION ---------------------------
  -q <FILE>,  --query=FILE (.fa or .fa.gz)
                 FASTA file for query sequence(s), use '-' for stdin
  -i <FILE>,  --index=FILE
                 pregenerated suffix array file for target sequence(s)
  -s n:m/l,   --seed=n:m/l
                 set seed length (-s l = length only; -s n:m = full interval;
                 -s n:m/l = length in interval; default -s 6)
  -l <int>,   --extension=L
                 max extension length(L) on the seed (do DP for max this length
                 up- and downstream of seed) (default L=20)
  -e <float>, --energy=dG
                 set deltaG energy threshold (in kcal/mol) to filter predictions
                 (default=-20)
  -z mat,     --matrix=mat
                 set energy matrix to t99 or t04 (default) for RNA-RNA duplexes
  -d <int>,   --penalty=dP
                 per-nucleotide extension penalty given in dacal/mol
                 (recommended: 30, default: 0)
  -t <int>,   --threads=N
                 set maximum number of threads to use (default=1)
  -p,         --report_alignment
                 report predictions in detailed format
  -p2,        --report_alignment=2
                 report predictions in a simple format together with CIGAR-like
                 string for interaction structure
  -p3,        --report_alignment=3
                 report predictions in a simple format together with
                 binding site (3'->5'), flanking 5'end (3'->5') and
                 flanking 3'end (5'->3') sequences of the target
                 (required for post-processing of CRISPR off-target predictions)
  --noGUseed     consider G-U wobble pairs as mismatch within the seed
                 (only for locating seeds, energy model is not affected)
  --verbose      verbose output
------------------------------------------------------------------------------
---------------------------- EXPERIMENTAL OPTIONS ----------------------------
  -m c:p,     --mismatch=c:p
                 introduce mismatched seeds
                 Set the max num of mismatches (c) allowed in the seed and
                 min num of consecutive matches required at seed start/end (p)
                 ! These seeds will not overlap with perfect complementary seeds
                 (default -m 0:0  (no mismatch);
                 if you set c>0, please also set p>0 to avoid overlaps)
  -x <float>, --seed_energy=F
                 set energy per length threshold that filters seeds (default=0)
                 
```



