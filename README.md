# 🧬 Subset Sum Optimization using Genetic Algorithm (GA)

This project implements a **Genetic Algorithm (GA)** to solve the classical **Subset Sum Problem (SSP)**. Instead of using brute-force techniques, it uses evolutionary principles like selection, crossover, and mutation to evolve a population of candidate solutions toward a target sum.

---

## 🧠 Problem Overview

> **Goal**: From a given set of integers, find a subset whose sum is exactly or closely equal to a target value.

The Subset Sum Problem is an NP-complete problem, making it an ideal candidate for heuristic approaches like Genetic Algorithms. This project demonstrates the effectiveness of GAs in approximating optimal solutions when exhaustive search is computationally expensive.

---

## 🧰 Key Features

- ✅ **Binary-encoded individuals** representing subset inclusion
- ✅ **Fitness evaluation** based on how close the subset sum is to the target
- ✅ **Tournament selection** to pick parents for reproduction
- ✅ **Single-point crossover** for combining parent chromosomes
- ✅ **Bit-flip mutation** for exploration
- ✅ **Real-time console logging**
- ✅ **Results exported to `GA_results.csv`**

---

## ⚙️ Genetic Algorithm Parameters

```java
int populationSize = 10;         // Number of individuals in the population
double crossoverRate = 0.7;      // Probability of crossover between pairs
double mutationRate = 0.05;      // Mutation rate per gene (bit-flip)
int maxGenerations = 50;         // Number of generations to run


Main.java                    # Entry point of the application
 ├── class GeneticAlgorithm  # Contains the core logic for GA operations
 ├── class Individual        # Represents a chromosome and evaluates fitness
GA_results.csv               # Output file with per-generation accuracy data

```

Sample Output

Generation 1: Best Fitness = 0.0192
Generation 2: Best Fitness = 0.0258
...
Perfect solution found in generation 27: 
Subset: [12, 28, 80, 50, 44, 56, 30]
Sum: 300
Fitness: 1.0
Average Accuracy: 96.34%
