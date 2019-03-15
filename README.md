# tRNA gene classifier

This program uses DNA data alone to predict tRNA gene expression, using binary (active/inactive) classifications. The premise of this project was based heavily on the observations made in this study: https://www.pnas.org/content/pnas/early/2018/08/17/1801240115.full.pdf. A separate manuscript focused on the development of this classifier will be posted on bioRxiv very soon.



To run, simply type ./tRNA

## List of command line flags:

##### --path: the path to the directory where the allPenaltiesPct.txt file is located (necessary for reading in mutation effects)
##### --print: generation printing interval (--print 500 will print every 500 generations, etc.)
##### -b: burn-in - output will start printing after a set number of generations pass (-b 0 prints results starting from beginning, -b 50000 prints only after 50000 generations have passed)
##### -n: population size - number of individuals, will be held constant throughout simulation, and each individual is diploid (default: 1,000)
##### -g: number of generations for simulation to run (default: 1,000,000)
##### --ug: germline mutation rate (default: 1e-6)
##### --us: somatic mutation rate (default: 1e-5)
##### --dup: duplication rate
##### --del: deletion rate
##### -s: seed, to ensure non-unique results on many simulations begun at the same time
##### --start: number of tRNA genes in the genome initially (default: 1)
##### --pseudo: begin with a tRNA pseudogene in addition to the number of functional tRNA genes
##### -c: fraction of duplications that are local (default: 1)
##### -m: length of chromosome in morgans
##### --run: name of a run, to name folders for optional output of each run
##### --output-frequencies: outputs frequency log file, which contains frequencies over time of each tRNA present at the end of the simulation
##### --output-lifespans: outputs lifespan log file, which contains frequency and lifespan data for all tRNAs created during the simulation

The following are based on the Nowak paper linked above:

##### --model-1:
Two tRNAs initially, with exactly the same mutation rate and function. All mutations are completely inactivating. No somatic mutations, duplications or deletions are possible.

##### --model-2:
Two tRNAs, one with function = 1 and mutation rate = --ug; one with function = 0.8 and mutation rate = --ug / 100. All mutations are completely inactivating. No somatic mutations, duplications or deletions are possible.

##### --model-4:
Any number of tRNAs, each with function = 1 and mutation rate == --ug, but somatic mutation rates are now possible. When invoking this model, individual fitness = 1 - ((developmental error rate)^(number of functional tRNA genes)). See Nowak paper for more detailed explanation.
###### --model-4-count: number of tRNA genes in genome initially in addition to --start value (only used with --model-4; default = 1).
###### --model-4-deverr: developmental error rate (only used with --model-4; default = 1e-4).

All questions should be directed to bthornlo@ucsc.edu.
