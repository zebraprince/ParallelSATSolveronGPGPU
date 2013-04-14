#ifndef __STRUCT_DEF__
#define __STRUCT_DEF__

//Constants
#define Sat 1
#define Unsat 0
#define Unknown -1
#define LitTrue 1
#define LitFalse 0
#define Unassigned -1

#define batchSize 12
#define Thread_per_Block 256



typedef struct {
	int countUnassignment;
	int clsState;
	int clsSize;
} Clause;

typedef struct {
	int formulaState;
	int numofCls;
	int numofLit;
	int numofTotalLit;
	Clause *cls;
	int *expPath;
	int *literal;
	int *litState;
	int *numofLitBeforeCls;
}Formula;

typedef struct {
	int *notifyUnit;	
	int numofFormula;
	Clause *cls_b;
	int *litState_b;
	int *formulaState_b;
	int *expPath_b;
	int *inferPath_b;
	int *notifyUnknown;
	int *notifyUnsat;
}FormulaBatch;


#endif
