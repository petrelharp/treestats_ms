import msprime
import libsequence
import allel
import timeit
import numpy as np

NE = 1e4
#RECRATE = 1e-8
RHO=1e4
RECRATE=1e3/(4*NE)
NU = [RECRATE]

def make_tree_sequence(nsam):
    """
    Make a tree seq w/no mutations
    """
    return msprime.simulate(nsam, Ne=NE, recombination_rate=RECRATE)


def mutate_tree_sequence(ts, Nu):
    """ 
    Copy a tree seq and mutate it, return
    mutated copy
    """
    tscopy = ts.dump_tables().tree_sequence()
    return msprime.mutate(tscopy, Nu)


def tskit_tajd(ts):
    return ts.Tajimas_D([[i for i in ts.samples()]])


def libseq_tajd_via_genotype_matrix(ts, pos):
    gm = ts.genotype_matrix()
    vm = libsequence.VariantMatrix(gm, pos)
    ac = vm.count_alleles()
    return libsequence.tajd(ac)


def allel_tajd(ts):
    gm = ts.genotype_matrix()
    ha = allel.HaplotypeArray(gm)
    ac = ha.count_alleles()
    return allel.tajima_d(ac)


print("toolkit nsam nmutations Nu nbytes seconds")
for nsam in np.logspace(np.log10(25), np.log10(1e6), 25).astype(np.int32):
    ts = make_tree_sequence(nsam)
    for Nu in NU:  # Skip no mutation as that is obviously dumb
        tscopy = mutate_tree_sequence(ts, Nu)
        # Conveniently, we can get the exact size
        # of the genotype matrix w/o ever calculating it
        nbytes = tscopy.num_sites * tscopy.num_samples
        timer = timeit.Timer("tskit_tajd(tscopy)", globals=globals())
        res = min(timer.repeat(repeat=5, number=1))
        p = tscopy.tables.sites.position
        print(f"tskit {nsam} {tscopy.num_mutations} {Nu} {nbytes} {res}")
        timer = timeit.Timer(
            "libseq_tajd_via_genotype_matrix(tscopy,p)", globals=globals())
        res = min(timer.repeat(repeat=5, number=1))
        print(f"libseq {nsam} {tscopy.num_mutations} {Nu} {nbytes} {res}")
        timer = timeit.Timer("allel_tajd(tscopy)", globals=globals())
        res = min(timer.repeat(repeat=5, number=1))
        print(f"allel {nsam} {tscopy.num_mutations} {Nu} {nbytes} {res}")
