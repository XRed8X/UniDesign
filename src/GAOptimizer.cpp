/*******************************************************************************************************************************
Copyright (c) Xiaoqiang Huang

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy,
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************************************************/

#pragma warning(disable:4244)
#pragma warning(disable:4305)
#include "GAOptimizer.h"
#include "ProteinDesign.h"
#include "Sequence.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <random>

// External flags from Main.cpp
extern BOOL FLAG_EVOLUTION;
extern BOOL FLAG_PPI;
extern BOOL FLAG_PROT_LIG;
extern BOOL FLAG_ENZYME;

extern char FILE_BESTSEQS[MAX_LEN_FILE_NAME + 1];
extern char FILE_BESTSTRUCT[MAX_LEN_FILE_NAME + 1];
extern char FILE_BEST_ALL_SITES[MAX_LEN_FILE_NAME + 1];
extern char FILE_BEST_MUT_SITES[MAX_LEN_FILE_NAME + 1];
extern char FILE_DESSEQS[MAX_LEN_FILE_NAME + 1];
extern char FILE_DESROT_NDX[MAX_LEN_FILE_NAME + 1];

// Random number generator for GA
static std::random_device rd;
static std::mt19937 gen(rd());

// ================================
// Context Management Functions
// ================================

int GAContextInitialize(GAContext* context, Structure* pStructure, RotamerList* pRotamerList)
{
    context->pStructure = pStructure;
    context->pRotamerList = pRotamerList;
    context->design_site_count = StructureGetDesignSiteCount(pStructure);
    context->use_evolution_score = FLAG_EVOLUTION;
    context->use_binding_energy = (FLAG_PPI || FLAG_PROT_LIG || FLAG_ENZYME);
    
    // Initialize rotamer type and count arrays
    context->ppRotamerTypes = (StringArray**)malloc(sizeof(StringArray*) * context->design_site_count);
    context->ppRotamerCounts = (IntArray**)malloc(sizeof(IntArray*) * context->design_site_count);
    
    for (int i = 0; i < context->design_site_count; i++)
    {
        context->ppRotamerTypes[i] = (StringArray*)malloc(sizeof(StringArray));
        context->ppRotamerCounts[i] = (IntArray*)malloc(sizeof(IntArray));
    }
    
    DesignSiteShowRotamerTypeAndCount(pRotamerList, pStructure, 
                                      context->ppRotamerTypes, 
                                      context->ppRotamerCounts);
    
    return Success;
}

void GAContextDestroy(GAContext* context)
{
    if (context->ppRotamerTypes != nullptr)
    {
        for (int i = 0; i < context->design_site_count; i++)
        {
            if (context->ppRotamerTypes[i] != nullptr)
            {
                StringArrayDestroy(context->ppRotamerTypes[i]);
                free(context->ppRotamerTypes[i]);
            }
        }
        free(context->ppRotamerTypes);
        context->ppRotamerTypes = nullptr;
    }
    
    if (context->ppRotamerCounts != nullptr)
    {
        for (int i = 0; i < context->design_site_count; i++)
        {
            if (context->ppRotamerCounts[i] != nullptr)
            {
                IntArrayDestroy(context->ppRotamerCounts[i]);
                free(context->ppRotamerCounts[i]);
            }
        }
        free(context->ppRotamerCounts);
        context->ppRotamerCounts = nullptr;
    }
}

// ================================
// Sequence Conversion Functions
// ================================

int GASolutionToSequence(const GASolution& solution, Sequence* sequence, GAContext* context)
{
    SequenceCreate(sequence);
    sequence->desSiteCount = context->design_site_count;
    
    // Copy rotamer indices from solution
    IntArrayResize(&sequence->rotNdxs, context->design_site_count);
    for (int i = 0; i < context->design_site_count; i++)
    {
        IntArraySet(&sequence->rotNdxs, i, solution.rotamer_indices[i]);
    }
    
    // Energy values will be calculated by EvaluateFitness
    sequence->etot = solution.total_energy;
    sequence->ephy = solution.physical_energy;
    sequence->ebin = solution.binding_energy;
    sequence->eevo = solution.evolution_energy;
    
    return Success;
}

int SequenceToGASolution(const Sequence* sequence, GASolution& solution, GAContext* context)
{
    solution.rotamer_indices.resize(context->design_site_count);
    
    for (int i = 0; i < context->design_site_count; i++)
    {
        solution.rotamer_indices[i] = IntArrayGet((IntArray*)&sequence->rotNdxs, i);
    }
    
    solution.total_energy = sequence->etot;
    solution.physical_energy = sequence->ephy;
    solution.binding_energy = sequence->ebin;
    solution.evolution_energy = sequence->eevo;
    
    return Success;
}

// ================================
// Fitness Evaluation
// ================================

int EvaluateFitness(GASolution& solution, GAContext* context)
{
    Sequence sequence;
    GASolutionToSequence(solution, &sequence, context);
    
    // Use existing UniDesign energy calculation
    int result = SequenceEnergy(context->pStructure, &sequence);
    
    if (result == Success)
    {
        // Update solution with calculated energies
        solution.total_energy = sequence.etot;
        solution.physical_energy = sequence.ephy;
        solution.binding_energy = sequence.ebin;
        solution.evolution_energy = sequence.eevo;
    }
    else
    {
        // If energy calculation fails, assign a penalty
        solution.total_energy = 1e10;
        solution.physical_energy = 1e10;
        solution.binding_energy = 0.0;
        solution.evolution_energy = 0.0;
    }
    
    SequenceDestroy(&sequence);
    return result;
}

// ================================
// Population Initialization
// ================================

int GenerateInitialPopulation(std::vector<GASolution>& population, 
                              GAContext* context, 
                              int population_size)
{
    population.resize(population_size);
    
    for (int i = 0; i < population_size; i++)
    {
        GASolution& individual = population[i];
        individual.rotamer_indices.resize(context->design_site_count);
        
        // Generate random rotamer indices for each design site
        for (int j = 0; j < context->design_site_count; j++)
        {
            int rotamer_count = context->pRotamerList->rotamerCount[j];
            
            // Count available rotamers
            int available_count = 0;
            for (int k = 0; k < rotamer_count; k++)
            {
                if (context->pRotamerList->remainFlag[j][k])
                {
                    available_count++;
                }
            }
            
            if (available_count > 0)
            {
                // Select a random available rotamer
                std::uniform_int_distribution<> dis(0, available_count - 1);
                int selected = dis(gen);
                
                // Find the actual rotamer index
                int count = 0;
                for (int k = 0; k < rotamer_count; k++)
                {
                    if (context->pRotamerList->remainFlag[j][k])
                    {
                        if (count == selected)
                        {
                            individual.rotamer_indices[j] = k;
                            break;
                        }
                        count++;
                    }
                }
            }
            else
            {
                individual.rotamer_indices[j] = 0;
            }
        }
        
        // Evaluate fitness for this individual
        EvaluateFitness(individual, context);
    }
    
    printf("Generated initial population of %d individuals\n", population_size);
    return Success;
}

// ================================
// Genetic Operators
// ================================

void Crossover(const GASolution& parent1, const GASolution& parent2,
               GASolution& offspring1, GASolution& offspring2,
               GAContext* context)
{
    int gene_count = context->design_site_count;
    offspring1.rotamer_indices.resize(gene_count);
    offspring2.rotamer_indices.resize(gene_count);
    
    // Single-point crossover
    std::uniform_int_distribution<> dis(1, gene_count - 1);
    int crossover_point = dis(gen);
    
    for (int i = 0; i < gene_count; i++)
    {
        if (i < crossover_point)
        {
            offspring1.rotamer_indices[i] = parent1.rotamer_indices[i];
            offspring2.rotamer_indices[i] = parent2.rotamer_indices[i];
        }
        else
        {
            offspring1.rotamer_indices[i] = parent2.rotamer_indices[i];
            offspring2.rotamer_indices[i] = parent1.rotamer_indices[i];
        }
    }
}

void Mutate(GASolution& solution, double mutation_rate, GAContext* context)
{
    std::uniform_real_distribution<> prob_dis(0.0, 1.0);
    
    for (int i = 0; i < context->design_site_count; i++)
    {
        if (prob_dis(gen) < mutation_rate)
        {
            // Mutate this position
            int rotamer_count = context->pRotamerList->rotamerCount[i];
            
            // Count available rotamers
            int available_count = 0;
            for (int k = 0; k < rotamer_count; k++)
            {
                if (context->pRotamerList->remainFlag[i][k])
                {
                    available_count++;
                }
            }
            
            if (available_count > 1)
            {
                // Select a different rotamer
                std::uniform_int_distribution<> dis(0, available_count - 1);
                int selected = dis(gen);
                
                int count = 0;
                for (int k = 0; k < rotamer_count; k++)
                {
                    if (context->pRotamerList->remainFlag[i][k])
                    {
                        if (count == selected)
                        {
                            solution.rotamer_indices[i] = k;
                            break;
                        }
                        count++;
                    }
                }
            }
        }
    }
}

// ================================
// Selection Operators
// ================================

int TournamentSelection(const std::vector<GASolution>& population, int tournament_size)
{
    std::uniform_int_distribution<> dis(0, population.size() - 1);
    
    int best_index = dis(gen);
    double best_fitness = population[best_index].total_energy;
    
    for (int i = 1; i < tournament_size; i++)
    {
        int candidate_index = dis(gen);
        double candidate_fitness = population[candidate_index].total_energy;
        
        // Lower energy is better
        if (candidate_fitness < best_fitness)
        {
            best_index = candidate_index;
            best_fitness = candidate_fitness;
        }
    }
    
    return best_index;
}

// ================================
// Output Functions
// ================================

void PrintGAProgress(int generation, double best_fitness, double avg_fitness)
{
    printf("Generation %4d: Best Energy = %12.4f, Avg Energy = %12.4f\n", 
           generation, best_fitness, avg_fitness);
}

int SaveBestSolution(const GASolution& best_solution, GAContext* context)
{
    Sequence best_sequence;
    GASolutionToSequence(best_solution, &best_sequence, context);
    
    // Save best sequence
    FILE* pFileSeq = fopen(FILE_BESTSEQS, "w");
    if (pFileSeq != nullptr)
    {
        SequenceWriteDesignFasta(&best_sequence, context->pStructure, 0, pFileSeq);
        fclose(pFileSeq);
        printf("Best sequence saved to %s\n", FILE_BESTSEQS);
    }
    
    // Save best structure
    DesignShowMinEnergyDesignStructure(context->pStructure, &best_sequence, FILE_BESTSTRUCT);
    printf("Best structure saved to %s\n", FILE_BESTSTRUCT);
    
    // Save design sites
    DesignShowMinEnergyDesignSites(context->pStructure, &best_sequence, FILE_BEST_ALL_SITES);
    
    // Save mutable sites
    DesignShowMinEnergyDesignMutableSites(context->pStructure, &best_sequence, FILE_BEST_MUT_SITES);
    
    // Save rotamer indices
    FILE* pFileRot = fopen(FILE_DESROT_NDX, "w");
    if (pFileRot != nullptr)
    {
        SequenceWriteDesignRotamer(&best_sequence, context->pStructure, 0, pFileRot);
        fclose(pFileRot);
    }
    
    SequenceDestroy(&best_sequence);
    return Success;
}

// ================================
// Main GA Optimization Function
// ================================

int RunGAOptimization(Structure* pStructure, RotamerList* pRotamerList, GAParameters params)
{
    printf("\n");
    printf("========================================\n");
    printf("  Genetic Algorithm Optimizer\n");
    printf("========================================\n");
    printf("Population Size:    %d\n", params.population_size);
    printf("Max Generations:    %d\n", params.max_generations);
    printf("Crossover Rate:     %.2f\n", params.crossover_rate);
    printf("Mutation Rate:      %.2f\n", params.mutation_rate);
    printf("Elite Count:        %d\n", params.elite_count);
    printf("Selection:          %s\n", params.use_tournament_selection ? "Tournament" : "Roulette");
    if (params.use_tournament_selection)
    {
        printf("Tournament Size:    %d\n", params.tournament_size);
    }
    printf("========================================\n\n");
    
    // Initialize GA context
    GAContext context;
    GAContextInitialize(&context, pStructure, pRotamerList);
    
    // Generate initial population
    std::vector<GASolution> population;
    GenerateInitialPopulation(population, &context, params.population_size);
    
    // Track best solution across all generations
    GASolution global_best = population[0];
    for (const auto& ind : population)
    {
        if (ind.total_energy < global_best.total_energy)
        {
            global_best = ind;
        }
    }
    
    printf("\nInitial best energy: %.4f\n\n", global_best.total_energy);
    
    // Main GA loop
    for (int generation = 0; generation < params.max_generations; generation++)
    {
        // Create new population
        std::vector<GASolution> new_population;
        new_population.reserve(params.population_size);
        
        // Elitism: preserve best individuals
        std::vector<GASolution> sorted_pop = population;
        std::sort(sorted_pop.begin(), sorted_pop.end(), 
                 [](const GASolution& a, const GASolution& b) {
                     return a.total_energy < b.total_energy;
                 });
        
        for (int i = 0; i < params.elite_count && i < sorted_pop.size(); i++)
        {
            new_population.push_back(sorted_pop[i]);
        }
        
        // Generate offspring through selection, crossover, and mutation
        while (new_population.size() < params.population_size)
        {
            // Selection
            int parent1_idx, parent2_idx;
            if (params.use_tournament_selection)
            {
                parent1_idx = TournamentSelection(population, params.tournament_size);
                parent2_idx = TournamentSelection(population, params.tournament_size);
            }
            else
            {
                // Simple random selection (can be improved with roulette wheel)
                std::uniform_int_distribution<> dis(0, population.size() - 1);
                parent1_idx = dis(gen);
                parent2_idx = dis(gen);
            }
            
            const GASolution& parent1 = population[parent1_idx];
            const GASolution& parent2 = population[parent2_idx];
            
            // Crossover
            std::uniform_real_distribution<> prob_dis(0.0, 1.0);
            GASolution offspring1, offspring2;
            
            if (prob_dis(gen) < params.crossover_rate)
            {
                Crossover(parent1, parent2, offspring1, offspring2, &context);
            }
            else
            {
                offspring1.rotamer_indices = parent1.rotamer_indices;
                offspring2.rotamer_indices = parent2.rotamer_indices;
            }
            
            // Mutation
            Mutate(offspring1, params.mutation_rate, &context);
            Mutate(offspring2, params.mutation_rate, &context);
            
            // Evaluate offspring
            EvaluateFitness(offspring1, &context);
            EvaluateFitness(offspring2, &context);
            
            new_population.push_back(offspring1);
            if (new_population.size() < params.population_size)
            {
                new_population.push_back(offspring2);
            }
        }
        
        // Replace old population
        population = new_population;
        
        // Update global best
        for (const auto& ind : population)
        {
            if (ind.total_energy < global_best.total_energy)
            {
                global_best = ind;
            }
        }
        
        // Calculate statistics
        double sum = 0.0;
        for (const auto& ind : population)
        {
            sum += ind.total_energy;
        }
        double avg_fitness = sum / population.size();
        
        // Print progress
        if (generation % 5 == 0 || generation == params.max_generations - 1)
        {
            PrintGAProgress(generation, global_best.total_energy, avg_fitness);
        }
    }
    
    printf("\n");
    printf("========================================\n");
    printf("  GA Optimization Complete\n");
    printf("========================================\n");
    printf("Final Best Energy:  %.4f\n", global_best.total_energy);
    printf("  Physical:         %.4f\n", global_best.physical_energy);
    printf("  Binding:          %.4f\n", global_best.binding_energy);
    printf("  Evolution:        %.4f\n", global_best.evolution_energy);
    printf("========================================\n\n");
    
    // Save best solution
    SaveBestSolution(global_best, &context);
    
    // Cleanup
    GAContextDestroy(&context);
    
    return Success;
}
