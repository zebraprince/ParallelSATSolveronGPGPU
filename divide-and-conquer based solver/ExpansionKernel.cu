#ifndef _EXPANSION_KERNEL_H_
#define _EXPANSION_KERNEL_H_

#include "struct_def.h"


__global__ void Expansion_Kernel(Formula origin, FormulaBatch s, FormulaBatch r, int expPara){

	int formulaIdx = blockIdx.y;
	int blockCol = blockIdx.x;
	int clauseCol = threadIdx.x;

	if(blockCol*Thread_per_Block + clauseCol < origin.numofCls){
		
		__shared__ int d_expPath[batchSize];

		for(int expCounter=0; expCounter<expPara; expCounter++){
		
			int clauseIdx_r = formulaIdx*origin.numofCls*expPara + expCounter*origin.numofCls + blockCol*Thread_per_Block+clauseCol;
			int clauseIdx_s = formulaIdx*origin.numofCls + blockCol*Thread_per_Block+clauseCol;
			r.cls_b[clauseIdx_r].clsSize = s.cls_b[clauseIdx_s].clsSize;
			r.cls_b[clauseIdx_r].countUnassignment = s.cls_b[clauseIdx_s].countUnassignment;
			r.cls_b[clauseIdx_r].clsState = s.cls_b[clauseIdx_s].clsState;
			
			d_expPath[expCounter] = s.expPath_b[formulaIdx*batchSize + expCounter];
			int currentClsSize = s.cls_b[clauseIdx_s].clsSize;
			r.formulaState_b[formulaIdx*expPara + expCounter] = s.formulaState_b[formulaIdx];
			__syncthreads();

			for(int litCounterinCls = 0; litCounterinCls < currentClsSize; litCounterinCls++){		
			
				int litIdx_s = formulaIdx*origin.numofTotalLit + origin.numofLitBeforeCls[blockCol*Thread_per_Block+clauseCol] + litCounterinCls;
				int litIdx_r = formulaIdx*origin.numofTotalLit*expPara + expCounter*origin.numofTotalLit + origin.numofLitBeforeCls[blockCol*Thread_per_Block+clauseCol] + litCounterinCls;
				int litIdx_origin = origin.numofLitBeforeCls[blockCol*Thread_per_Block+clauseCol] + litCounterinCls;

				r.litState_b[litIdx_r] = s.litState_b[litIdx_s];
				
				if(r.cls_b[clauseIdx_r].clsState == Unknown)	
					for(int pathCounter=0; pathCounter<=expCounter; pathCounter++){

						if( (r.litState_b[litIdx_r] == Unassigned) ){
							if( ((origin.literal[litIdx_origin]+d_expPath[pathCounter]==0) && (pathCounter==expCounter)) || ((origin.literal[litIdx_origin]-d_expPath[pathCounter]==0) && (pathCounter!=expCounter)) ){
								r.litState_b[litIdx_r] = LitFalse;
								r.cls_b[clauseIdx_r].countUnassignment--;
								if((r.cls_b[clauseIdx_r].countUnassignment == 0) && (r.cls_b[clauseIdx_r].clsState == Unknown))
									r.cls_b[clauseIdx_r].clsState = Unsat;
							}
							else if( ((origin.literal[litIdx_origin]+d_expPath[pathCounter]==0) && (pathCounter!=expCounter)) || ((origin.literal[litIdx_origin]-d_expPath[pathCounter]==0) && (pathCounter==expCounter)) ){
								r.litState_b[litIdx_r] = LitTrue;	
								r.cls_b[clauseIdx_r].countUnassignment--;
								if(r.cls_b[clauseIdx_r].clsState == Unknown)
									r.cls_b[clauseIdx_r].clsState = Sat;
							}
						}
					}
			}	

		}
	}
}
#endif

