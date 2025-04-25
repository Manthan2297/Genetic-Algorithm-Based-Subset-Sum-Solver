
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Main {

    public static void main(String[] args) {
        // Define the set and target sum
        List<Integer> numberList = List.of(3, 34, 4, 12, 5, 2, 25, 31, 60, 91, 47, 73, 17, 53, 28, 39, 67, 80, 36, 50, 15, 95, 44, 78, 20, 10, 13, 56, 89, 14, 38, 70, 9, 40, 22, 7, 76, 58, 49, 85);
        int targetSum = 300;

        // GA Parameters
        int populationSize = 10;
        double crossoverRate = 0.7;
        double mutationRate = 0.05;
        int maxGenerations = 50;

        // Initialize GA and run
        GeneticAlgorithm ga = new GeneticAlgorithm(numberList, targetSum, populationSize, crossoverRate, mutationRate, maxGenerations);
        Individual bestSolution = ga.run();

        // Print the best solution
        System.out.println("Best Solution: " + bestSolution);
    }
}

class Individual {

    private List<Integer> chromosome;
    private double fitness;
    private List<Integer> numberList;
    private int targetSum;

    // Constructor
    public Individual(List<Integer> numberList, int targetSum, List<Integer> chromosome) {
        this.numberList = numberList;
        this.targetSum = targetSum;
        this.chromosome = chromosome;
        calculateFitness();
    }

    // Calculate fitness based on the absolute difference between subset sum and target
    public void calculateFitness() {
        int subsetSum = getSubsetSum();
        int diff = Math.abs(targetSum - subsetSum);
        this.fitness = 1.0 / (diff + 1);
    }

    // Get the fitness value
    public double getFitness() {
        return fitness;
    }

    // Get the subset sum represented by the chromosome
    public int getSubsetSum() {
        int sum = 0;
        for (int i = 0; i < chromosome.size(); i++) {
            if (chromosome.get(i) == 1) {
                sum += numberList.get(i);
            }
        }
        return sum;
    }

    // Return chromosome
    public List<Integer> getChromosome() {
        return chromosome;
    }

    // Return subset represented by the chromosome
    public List<Integer> getSubset() {
        List<Integer> subset = new ArrayList<>();
        for (int i = 0; i < chromosome.size(); i++) {
            if (chromosome.get(i) == 1) {
                subset.add(numberList.get(i));
            }
        }
        return subset;
    }

    @Override
    public String toString() {
        return "Subset: " + getSubset() + ", Sum: " + getSubsetSum() + ", Fitness: " + fitness;
    }

    // Create a copy of the individual
    public Individual copy() {
        return new Individual(numberList, targetSum, new ArrayList<>(chromosome));
    }
}

class GeneticAlgorithm {

    private List<Integer> numberList;
    private int targetSum;
    private int populationSize;
    private double crossoverRate;
    private double mutationRate;
    private int maxGenerations;
    private List<Individual> population;
    private Individual bestIndividual;

    // Constructor
    public GeneticAlgorithm(List<Integer> numberList, int targetSum, int populationSize, double crossoverRate, double mutationRate, int maxGenerations) {
        this.numberList = numberList;
        this.targetSum = targetSum;
        this.populationSize = populationSize;
        this.crossoverRate = crossoverRate;
        this.mutationRate = mutationRate;
        this.maxGenerations = maxGenerations;
        this.population = new ArrayList<>();
    }

    // Run the genetic algorithm
    public Individual run() {
        double totalAccuracyAcrossGenerations = 0;
        int iterationCount = 0;

        try (FileWriter csvWriter = new FileWriter("GA_results.csv")) {
            csvWriter.append("Generation,Accuracy\n");

            initializePopulation();
            evaluatePopulation();

            for (int generation = 1; generation <= maxGenerations; generation++) {
                List<Individual> selectedIndividuals = selection();
                List<Individual> offspring = crossover(selectedIndividuals);
                mutation(offspring);
                population = offspring;
                evaluatePopulation();

                // Calculate and track accuracy
                double accuracy = calculateAccuracy(bestIndividual.getSubsetSum(), targetSum);
                totalAccuracyAcrossGenerations += accuracy;
                iterationCount++;

                // Write the generation and accuracy to CSV
                csvWriter.append(generation + "," + accuracy + "\n");

                // Display the best fitness for the current generation
                System.out.println("Generation " + generation + ": Best Fitness = " + bestIndividual.getFitness());

                // Print when a perfect solution is found, but don't stop the algorithm
                if (bestIndividual.getSubsetSum() == targetSum) {
                    System.out.println("Perfect solution found in generation " + generation + ": " + bestIndividual);
                }
            }

            // Calculate the average accuracy across generations
            double averageAccuracy = totalAccuracyAcrossGenerations / iterationCount;
            System.out.println("Average Accuracy: " + averageAccuracy);
            csvWriter.append("Average Accuracy," + averageAccuracy + "\n");

        } catch (IOException e) {
            System.out.println("Error writing to CSV file: " + e.getMessage());
        }

        return bestIndividual;
    }

    // Initialize the population with random individuals
    private void initializePopulation() {
        Random random = new Random();
        for (int i = 0; i < populationSize; i++) {
            List<Integer> chromosome = new ArrayList<>();
            for (int j = 0; j < numberList.size(); j++) {
                chromosome.add(random.nextBoolean() ? 1 : 0);
            }
            population.add(new Individual(numberList, targetSum, chromosome));
        }
    }

    // Evaluate the fitness of each individual in the population
    private void evaluatePopulation() {
        for (Individual individual : population) {
            individual.calculateFitness();
            updateBestIndividual(individual);
        }
    }

    // Update the best individual found so far
    private void updateBestIndividual(Individual individual) {
        if (bestIndividual == null || individual.getFitness() > bestIndividual.getFitness()) {
            bestIndividual = individual.copy();
        }
    }

    // Perform selection using tournament selection
    private List<Individual> selection() {
        List<Individual> selected = new ArrayList<>();
        Random random = new Random();
        for (int i = 0; i < populationSize; i++) {
            Individual best = null;
            for (int j = 0; j < 3; j++) {  // Tournament size of 3
                Individual contender = population.get(random.nextInt(population.size()));
                if (best == null || contender.getFitness() > best.getFitness()) {
                    best = contender;
                }
            }
            selected.add(best);
        }
        return selected;
    }

    // Perform single-point crossover
    private List<Individual> crossover(List<Individual> selectedIndividuals) {
        Random random = new Random();
        List<Individual> offspring = new ArrayList<>();

        for (int i = 0; i < selectedIndividuals.size(); i += 2) {
            Individual parent1 = selectedIndividuals.get(i);
            Individual parent2 = selectedIndividuals.get(Math.min(i + 1, selectedIndividuals.size() - 1));

            List<Integer> childChromosome1 = new ArrayList<>(parent1.getChromosome());
            List<Integer> childChromosome2 = new ArrayList<>(parent2.getChromosome());

            if (random.nextDouble() < crossoverRate) {
                int crossoverPoint = random.nextInt(numberList.size());
                for (int j = 0; j < crossoverPoint; j++) {
                    childChromosome1.set(j, parent2.getChromosome().get(j));
                    childChromosome2.set(j, parent1.getChromosome().get(j));
                }
            }

            offspring.add(new Individual(numberList, targetSum, childChromosome1));
            offspring.add(new Individual(numberList, targetSum, childChromosome2));
        }

        return offspring;
    }

    // Perform bit-flip mutation
    private void mutation(List<Individual> offspring) {
        Random random = new Random();
        for (Individual individual : offspring) {
            for (int i = 0; i < individual.getChromosome().size(); i++) {
                if (random.nextDouble() < mutationRate) {
                    int currentValue = individual.getChromosome().get(i);
                    individual.getChromosome().set(i, currentValue == 0 ? 1 : 0);
                }
            }
        }
    }

    // Calculate accuracy as the percentage difference from the target sum
    private double calculateAccuracy(int subsetSum, int targetSum) {
        if (subsetSum <= targetSum) {
            return (subsetSum / (double) targetSum) * 100;
        } else {
            return (targetSum / (double) subsetSum) * 100;
        }
    }
}
