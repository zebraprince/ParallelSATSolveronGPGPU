#ifndef __STRUCT_DEF__
#define __STRUCT_DEF__

//Constants
#define SAT 1
#define UNSAT 0
#define UNKNOWN -1
#define UNSET 0
#define SET 1
#define batchSize 192
#define thread_per_block 128
//constant probability of a random walk
const int P = 6;//probability 0.6
//random walk after trapped by trap times of flipping with local best
const int Trapped = 4;

const int RestartTimes = 4;

typedef struct {
	int clsSize;
} Clause;

typedef struct {
	int formulaState;
	int numofCls;
	int numofVar;
	int numofTotalLit;
	Clause *cls;
	int *literal;
	int *numofLitBeforeCls;
	int MAX_CLS_SIZE;
}Formula;

typedef struct {
	int var;
	int val;
	int weight;
}GuideVariable;

#endif
