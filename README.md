# pro2. calculate the Sum-Of-Pair score of the MSA
* 曹柏泓
* 113753108

## Description

* Write a Python program to calculate the multiple sequence alignment's sum-of-pair score (SoP).
* Creating your own program, i.e., hw2.py.
* Packages you can use: numpy, pandas
* You write a program with a function named calculate_SoP, ie.
```
def calculate_SoP(input_path, score_path, gopen, gextend):
    .
    .
    .
    .
```

## File

* hw2_ref.py: You can start from this reference code and try to write your own comment in English
* pam100.txt
* pam250.txt
* test1.fasta

## Parameters

* input_path: fasta file (ex. test1.fasta)
* score_path: score file (ex. pam250.txt)
* gopen: gap open penalty
* gextend: gap extend penalty

## Command

Please go ahead and execute your code using the following command.


```Python
calculate_SoP("examples/test1.fasta", "pam250.txt", -10, -2) #score=1047
calculate_SoP("examples/test2.fasta", "pam100.txt", -8, -2) #score=606
```
 

## Evaluation

10 testing data(5 public, 5 private)

The correct answer gets 10 points for each testing data.

### Penalty

* High code similarity to others: YOUR SCORE = 0

## References
Please provide the code along with its reference. For example, you can cite it as: ```# ChatGPT, respond to “your prompt,” on February 16, 2023```. Below is an example of a reference format summarizing the use of ChatGPT for R programming

>You are the R Language expert.
>Please help me to write a function called “k_fold”.
>Using a given dataset to train the random forest model, and using the k-fold cross-validation to evaluate the best model parameters. Here is the instruction for the function requirements:\
>Function name: k_fold\
>Function parameters:
