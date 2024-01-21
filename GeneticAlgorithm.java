import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class GeneticAlgorithm {


    private static final int populationSize = 50;
    private static final double crossoverRate = 0.7;
    private static final double mutationRate = 0.001;
    private static final int genomeLength = 10;
    private static int totalRuns = 30;
	private static boolean bestFitnessFlag = false;
	
	public static void main(String args[]) {	   
	   runGA(populationSize, crossoverRate, mutationRate);
    }
	
	public static void runGA(int populationSize, double crossoverRate, double mutationRate) {
        
        int generation = 0; 		
	    // Step 1 - Make Population 
        List<String> population = makePopulation(populationSize, genomeLength);
		 
		 System.out.println("Population Size : "+populationSize);
		 System.out.println("Genome Length : "+genomeLength);
		 
		while (generation < totalRuns && !bestFitnessFlag) { 	// End loop when 30 runs are complete or the best fitness reaches 10. 
		// Step 2 - Evaluate Fitness 
			List<Double> fitness = evaluateFitness(population); 
		 
		// Step 3 - Find Average and Best Fitness; 
			List<Double> avgAndBest = AvgandBestFitness(fitness);
			if (avgAndBest.get(1) == 10.0) 
				bestFitnessFlag = true;
		 System.out.println("Generation "+generation+" : average fitness "+avgAndBest.get(0)+", best fitness "+ avgAndBest.get(1));
		 
		// Step 4 - Select Parents (Roulette Wheel selection)  		 
		// Step 5 - Crossover 
		    //List<String> crossoverGenomes = crossover(parents.get(0), parents.get(1));
		    
		    List<String> offspring = new ArrayList<>();
			while (offspring.size() < populationSize) {
				List<String> parents = selectPair(population);
				offspring.addAll(crossover(parents.get(0), parents.get(1)));
			}
		 
		// Step 6 - Mutation
	 
		 List<String> mutatedOffspring = new ArrayList<>();
			for (String child : offspring) {
				mutatedOffspring.add(mutate(child, mutationRate));
			}

			// Replace the current population with the combined population of offspring and
			// parents
			population = new ArrayList<>(mutatedOffspring);		 
		 
		 generation++;
		}
	}
			
	// Returns a randomly generated population of a given size, represented as a list of genomes of the specified length
	public static List<String> makePopulation(int size, int length) {
		ArrayList<String> populationArray = new ArrayList<>();
		for(int i=0; i<size; i++) {
			String generatedGenome = null;
			try {
				generatedGenome = randomGenome(length);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			populationArray.add(i, generatedGenome);	
		}
        return populationArray;		
	}
	
	// Returns the fitness for each genome of a given population
	public static List<Double> evaluateFitness(List<String> population) {
		// Evaluate fitness for each genome 
		List<Double> populationFitness = new ArrayList<>();
		
		for (int i = 0; i < population.size(); i++) {
		   	double fitnessGenome = fitness(population.get(i));
			populationFitness.add(i, fitnessGenome);
		}	
		
		return populationFitness;
	}
	
	// Returns Average and Best fitness of a given population
	public static List<Double> AvgandBestFitness(List<Double> populationFitness) {
		
		List<Double> AvgBestFitness = new ArrayList<>();
		double averageFitness; 
		double bestFitness = 0.0;
		double sum = 0.0;
		
		// Evaluate Average & Best fitness  
		for (int j = 0; j < populationFitness.size(); j++) {
			sum += populationFitness.get(j);		
			
			// best fitness 
			bestFitness = Math.max(bestFitness, populationFitness.get(j));
		}	

        // calculate average fitness 
        averageFitness = sum / populationFitness.size();

        AvgBestFitness.add(0, averageFitness);
        AvgBestFitness.add(1, bestFitness); 

        return AvgBestFitness;		
	}
	
	// Selects two genomes using fitness proportionate selection (Roulette-Wheel selection)
	private static List<String> selectPair(List<String> population) {
		List<String> selectedPair = new ArrayList<>();
		List<Double> fitnessValues = evaluateFitness(population);

		// Select the first parent
		int parent1Index = rouletteWheelSelection(population, fitnessValues);
		selectedPair.add(population.get(parent1Index));

		// Select the second parent (ensure it's different from the first)
		int parent2Index;
		do {
			parent2Index = rouletteWheelSelection(population, fitnessValues);
		} while (parent2Index == parent1Index);

		selectedPair.add(population.get(parent2Index));

		return selectedPair;
	}
	
	
	public static double calculateTotalFitness(List<Double> fitnessValues) {
        double totalFitness = 0.0;
        for (double d : fitnessValues) {
            totalFitness += d; 
        }
        return totalFitness;
    }
	

	public static int rouletteWheelSelection(List<String> population, List<Double> fitnessValues) {
		double totalFitness = calculateTotalFitness(fitnessValues);
		double randomValue = Math.random() * totalFitness;

		double cumulativeFitness = 0;
		for (int i = 0; i < population.size(); i++) {
			cumulativeFitness += fitnessValues.get(i);
			if (cumulativeFitness >= randomValue) {
				return i; // Return the index of the selected person
			}
		}

		// If nothing, return the last person 
		return population.size() - 1;
    }
	
	// Crossover function with a random crossover point and a specified crossover rate of 0.7
	public static List<String> crossover(String parent1, String parent2) {
		Random random = new Random();

		// Check if crossover should occur based on the crossover rate
		if (random.nextDouble() <= crossoverRate) {
			int crossoverPoint = random.nextInt(parent1.length());

			String child1 = parent1.substring(0, crossoverPoint) + parent2.substring(crossoverPoint);
			String child2 = parent2.substring(0, crossoverPoint) + parent1.substring(crossoverPoint);

			return Arrays.asList(child1, child2);
		} else {
			// If crossover doesn't occur, return the parents unchanged
			return Arrays.asList(parent1, parent2);
		}
	}
	

	
	public static String mutate(String genome, double mutationRate) {
		Random random = new Random();
		StringBuilder mutatedGenome = new StringBuilder();

		for (int i = 0; i < genome.length(); i++) {
			char currentChar = genome.charAt(i);
			if (random.nextDouble() < mutationRate) {
				// Flip the bit with probability mutationRate
				mutatedGenome.append((currentChar == '0') ? '1' : '0');
			} else {
				mutatedGenome.append(currentChar);
			}
		}

		return mutatedGenome.toString();
	}
		
	
	// Creates a random genome of given length
	public static String randomGenome(int genomeLength){
		Random random = new Random();
		StringBuilder genome = new StringBuilder();
        for (int i = 0; i < genomeLength; i++) {            
			genome.append(random.nextInt(2));
        }
        return genome.toString();
    
	}	
	
	// Returns the fitness of a genome 
	public static double fitness(String genome) {
        double fitness = 0.0;
		for(char c : genome.toCharArray()) {
			if(c == '1')
				fitness++;
		}
		return fitness;
	}	
}

	