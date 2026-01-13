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

#ifndef GA_OPTIMIZER_H
#define GA_OPTIMIZER_H

#include "Structure.h"
#include "Sequence.h"
#include "ProteinDesign.h"
#include "openGA.hpp"
#include <vector>

// ================================
// GA Solution Representation
// ================================
// Each individual in the GA represents a sequence design solution
// Genes encode rotamer indices for each design site
struct GASolution
{
    std::vector<int> rotamer_indices;  // Rotamer index for each design site
    
    // Fitness components (calculated by UniDesign)
    double total_energy;               // Total fitness (optimized by GA)
    double physical_energy;            // Physics-based energy
    double binding_energy;             // Protein-protein or protein-ligand binding
    double evolution_energy;           // PSSM-based evolutionary score
    
    // Constructor
    GASolution() : total_energy(0.0), physical_energy(0.0), 
                   binding_energy(0.0), evolution_energy(0.0) {}
};

// ================================
// GA Parameters Configuration
// ================================
struct GAParameters
{
    int population_size;               // Number of individuals in population
    int max_generations;               // Maximum number of generations
    double crossover_rate;             // Probability of crossover (0.0 - 1.0)
    double mutation_rate;              // Probability of mutation (0.0 - 1.0)
    int elite_count;                   // Number of elite individuals to preserve
    bool use_tournament_selection;     // Use tournament vs roulette selection
    int tournament_size;               // Size of tournament if used
    
    // Default constructor with recommended values
    GAParameters() : 
        population_size(100), 
        max_generations(50),
        crossover_rate(0.8),
        mutation_rate(0.1),
        elite_count(2),
        use_tournament_selection(true),
        tournament_size(3) {}
};

// ================================
// GA-UniDesign Bridge Context
// ================================
// This structure holds all UniDesign context needed for fitness evaluation
struct GAContext
{
    Structure* pStructure;             // Protein structure
    RotamerList* pRotamerList;         // Available rotamers
    StringArray** ppRotamerTypes;      // Rotamer type names
    IntArray** ppRotamerCounts;        // Rotamer counts per site
    int design_site_count;             // Number of design sites
    
    // Energy calculation flags
    bool use_evolution_score;          // Include PSSM/evolutionary terms
    bool use_binding_energy;           // Calculate binding energy (PPI/PLI)
    
    GAContext() : pStructure(nullptr), pRotamerList(nullptr),
                  ppRotamerTypes(nullptr), ppRotamerCounts(nullptr),
                  design_site_count(0), use_evolution_score(false),
                  use_binding_energy(false) {}
};

// ================================
// Main GA Optimization Functions
// ================================

/**
 * @brief Main entry point for GA-based protein design optimization
 * 
 * This function replaces SimulatedAnnealing in the ProteinDesign workflow
 * when --optimizer=GA is specified.
 * 
 * @param pStructure Protein structure with design sites
 * @param pRotamerList List of available rotamers for design sites
 * @param params GA algorithm parameters
 * @return Success code
 */
int RunGAOptimization(Structure* pStructure, RotamerList* pRotamerList, GAParameters params);

/**
 * @brief Initialize GA context from UniDesign structures
 * 
 * @param context Output GA context
 * @param pStructure Input protein structure
 * @param pRotamerList Input rotamer list
 * @return Success code
 */
int GAContextInitialize(GAContext* context, Structure* pStructure, RotamerList* pRotamerList);

/**
 * @brief Free GA context resources
 * 
 * @param context Context to destroy
 */
void GAContextDestroy(GAContext* context);

/**
 * @brief Generate a random initial population of solutions
 * 
 * @param population Output population vector
 * @param context GA context with design information
 * @param population_size Number of individuals to create
 * @return Success code
 */
int GenerateInitialPopulation(std::vector<GASolution>& population, 
                              GAContext* context, 
                              int population_size);

/**
 * @brief Evaluate fitness of a single solution using UniDesign energy functions
 * 
 * This is the critical bridge function that converts GA solutions to
 * Sequence objects and evaluates them using SequenceEnergy and related functions.
 * 
 * @param solution Solution to evaluate (will be modified with fitness values)
 * @param context GA context with structure and rotamer information
 * @return Success code
 */
int EvaluateFitness(GASolution& solution, GAContext* context);

/**
 * @brief Convert GASolution to UniDesign Sequence for energy evaluation
 * 
 * @param solution Input GA solution
 * @param sequence Output UniDesign sequence
 * @param context GA context
 * @return Success code
 */
int GASolutionToSequence(const GASolution& solution, Sequence* sequence, GAContext* context);

/**
 * @brief Convert UniDesign Sequence back to GASolution
 * 
 * @param sequence Input UniDesign sequence
 * @param solution Output GA solution
 * @param context GA context
 * @return Success code
 */
int SequenceToGASolution(const Sequence* sequence, GASolution& solution, GAContext* context);

/**
 * @brief Perform crossover operation between two parent solutions
 * 
 * Uses single-point or uniform crossover to create offspring
 * 
 * @param parent1 First parent
 * @param parent2 Second parent
 * @param offspring1 First offspring (output)
 * @param offspring2 Second offspring (output)
 * @param context GA context
 */
void Crossover(const GASolution& parent1, const GASolution& parent2,
               GASolution& offspring1, GASolution& offspring2,
               GAContext* context);

/**
 * @brief Perform mutation operation on a solution
 * 
 * Randomly changes rotamer indices at one or more design sites
 * 
 * @param solution Solution to mutate (modified in place)
 * @param mutation_rate Probability of mutating each gene
 * @param context GA context
 */
void Mutate(GASolution& solution, double mutation_rate, GAContext* context);

/**
 * @brief Select parent using tournament selection
 * 
 * @param population Current population
 * @param tournament_size Number of individuals in tournament
 * @return Index of selected parent
 */
int TournamentSelection(const std::vector<GASolution>& population, int tournament_size);

/**
 * @brief Save best solution from GA optimization
 * 
 * Writes the best sequence and structure using UniDesign output functions
 * 
 * @param best_solution Best solution found
 * @param context GA context
 * @return Success code
 */
int SaveBestSolution(const GASolution& best_solution, GAContext* context);

/**
 * @brief Print GA optimization progress
 * 
 * @param generation Current generation number
 * @param best_fitness Best fitness in current generation
 * @param avg_fitness Average fitness in population
 */
void PrintGAProgress(int generation, double best_fitness, double avg_fitness);

#endif // GA_OPTIMIZER_H
