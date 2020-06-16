## Description
Heredity is an algorithm that calculates the probability distribution of a person having a particular number of specified genes and genetic trait given their family's genetic inheritance and their traits using a Bayesian network structure
- a modifiable dictionary, PROBS, is created to store the numbers of genes, traits, and their corresponding unconditional probability that will be used in calculation
- chance of genetic mutation is also considered in the calculation and is defaulted to be 0.01

## Usage
```
$ python heredity.py [file]
```
- file should be a formatted csv file
  - column should follow name, mother, father, trait
  - see example

```
$ python heredity example
Harry:
  Gene:
    2: 0.0092    
    1: 0.4557    
    0: 0.5351    
  Trait:
    True: 0.2665 
    False: 0.7335
James:
  Gene:
    2: 0.1976
    1: 0.5106
    0: 0.2918
  Trait:
    True: 1.0000
    False: 0.0000
Lily:
  Gene:
    2: 0.0036
    1: 0.0136
    0: 0.9827
  Trait:
    True: 0.0000
    False: 1.0000
```
