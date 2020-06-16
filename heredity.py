import csv
import itertools
import sys
import numpy

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    # initialize probability datas
    joint_pb = 0
    mutation = PROBS["mutation"]

    # create empty dictionaries to store probability distributions for each gene type
    pb_one_gene = {}
    pb_two_genes = {}
    pb_no_gene = {}

    # default probabilities all to 1 for later multiplication purposes
    for person in people:
        pb_one_gene[person] = 1
        pb_two_genes[person] = 1
        pb_no_gene[person] = 1

    # loop over each individual and calculate for their corresponding probability
    for person in people:
        # remember the individual in a convenient manner
        x = people[person]
        
        # remember the gene numbers of mother and father
        # then calculate for the probabilities of genes inheritance based on the 
        # genes of mother and father, considering chance for mutation
        if not x["mother"] == None:
            if x["mother"] in one_gene:
                mother_genes = 1
            elif x["mother"] in two_genes:
                mother_genes = 2
            else:
                mother_genes = 0
            pb_from_mother = ((1 - mother_genes / 2) * mutation +
                              mother_genes / 2 * (1 - mutation))
        if not x["father"] == None:
            if x["father"] in one_gene:
                father_genes = 1
            elif x["father"] in two_genes:
                father_genes = 2
            else:
                father_genes = 0
            pb_from_father = ((1 - father_genes / 2) * mutation +
                              father_genes / 2 * (1 - mutation))

        # if the person belongs to the one_gene group
        if person in one_gene:
            # if the data for their parents is missing, set the probability of getting one gene 
            # to the unconditional proability times the probability of whether they have trait or not
            if x["mother"] == None and x["father"] == None:
                pb_one_gene[person] = PROBS["gene"][1]
            
            # calculate for the probability of the individual obtaining one gene based on the genes
            # inheritance probabilities from their parent
            else:
                pb_one_gene[person] = (pb_from_mother * (1 - pb_from_father) +
                                       (1 - pb_from_mother) * pb_from_father)

            # multiply the probailities of maintaining the genes while manifesting the traits
            if person in have_trait:
                pb_one_gene[person] *= PROBS["trait"][1][True]
            else:
                pb_one_gene[person] *= PROBS["trait"][1][False]

        # if the person belongs to the two_genes group
        elif person in two_genes:
            # if the data for their parents is missing, set the probability of getting two genes 
            # to the unconditional proability times the probability of whether they have trait or not
            if x["mother"] == None and x["father"] == None:
                pb_two_genes[person] = PROBS["gene"][2]

            # calculate for the probability of the individual obtaining two genes based on the genes
            # inheritance probabilities from their parent
            else:
                pb_two_genes[person] = pb_from_mother * pb_from_father
                
            if person in have_trait:
                pb_two_genes[person] *= PROBS["trait"][2][True]
            else:
                pb_two_genes[person] *= PROBS["trait"][2][False]

        # the person then must belong to the no_gene group         
        else:
            # if the data for their parents is missing, set the probability of getting no gene 
            # to the unconditional proability times the probability of whether they have trait or not
            if x["mother"] == None and x["father"] == None: 
                pb_no_gene[person] = PROBS["gene"][0]

            # calculate for the probability of the individual obtaining no gene at all based on the genes
            # inheritance probabilities from their parent
            else:
                pb_no_gene[person] = (1 - pb_from_mother) * (1 - pb_from_father)
                
            if person in have_trait:
                pb_no_gene[person] *= PROBS["trait"][0][True]
            else:
                pb_no_gene[person] *= PROBS["trait"][0][False]

    # calculate for the joint probability by multiplying all individual probabilities together
    joint_pb = (numpy.prod(list(pb_one_gene.values())) * 
                numpy.prod(list(pb_two_genes.values())) * 
                numpy.prod(list(pb_no_gene.values())))
    return joint_pb


def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    for person in probabilities:
        if person in one_gene:
            probabilities[person]["gene"][1] += p
        elif person in two_genes:
            probabilities[person]["gene"][2] += p
        else:
            probabilities[person]["gene"][0] += p
        if person in have_trait:
            probabilities[person]["trait"][True] += p
        else:
            probabilities[person]["trait"][False] += p

def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    for person in probabilities:
        # normalize probability distribution of genes
        sum_pb_gene = sum(list(probabilities[person]["gene"].values()))
        ratio = 1 / sum_pb_gene
        for i in range(0, 3):
            probabilities[person]["gene"][i] *= ratio

        # normalize probability distribution of trait
        sum_pb_trait = sum(list(probabilities[person]["trait"].values()))
        ratio2 = 1 / sum_pb_trait
        for j in range(0, 2):
            probabilities[person]["trait"][j] *= ratio2

if __name__ == "__main__":
    main()
