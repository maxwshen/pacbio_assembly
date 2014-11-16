#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define NEITHER       0
#define UP            1
#define LEFT          2
#define UP_AND_LEFT   3
#define LENGTH        1000
#define ARTIFCAL      50 

int COV;

int main(int argc, char* argv[]) {
	FILE *fr_reads, *fr_genome;
	char ReadsFile[100];
	char GenomeFile[100];
	FILE *fw_sample4;
	FILE *fw_sample3;
	char CorrectedGenomeFile[100];
	char Top3GenomeFile[100];
	char Seq_header[100];
   	char line[LENGTH];
	char Read[100][LENGTH];
	char Sample4[100][LENGTH];
	//char Sample3[100][LENGTH];
	int file_i;
	int file_j;
	int iii,jjj;
	int ii;
	int i1,i2,i3,i4;
	int file_end_tag;
	int header_tag;
	int exit_tag;
	int Sample4_start[100];
	int Sample4_index[100][4];
	int Sample4_end[100];
	int Score[100];
	int M1_Score, M2_Score;
	int ScoreCov;
	char Top3Read[LENGTH];
	int Max_i1, Max_i2, Max_i3;
	int temp_Read_start;
	int temp_Read_end;
	int temp_score;
	char temp_Read[LENGTH];
	//srand(time(NULL));
	srand(1129);	
	//int r = rand();

	for (file_i = 0;file_i < 1;file_i++){
			sprintf(ReadsFile,"%s",argv[1]);
			//sprintf(GenomeFile,"genome.%d.fasta",file_i);
			fr_reads = fopen(ReadsFile, "r");  
			//fr_genome = fopen(GenomeFile, "r");  
			if (fr_reads != NULL){	
				printf("Reads file %s \n", ReadsFile);		
				iii = 0;	
				file_end_tag = 1;
				header_tag = 1;
				for (ii=0; ii++; ii<100)
					Score[ii] = 0;
				while(file_end_tag){ // read the reads file
					line[0] = 'X';
					line[1] = 'Y';
					fgets(line, LENGTH, fr_reads);
					if (header_tag == 1){
						//sprintf(Seq_header, "%s", line);
						header_tag = 0;
					}
					if ((line[0] == 'X')||(line[1] == 'Y'))
						file_end_tag = 0;
					else{
						fgets(Read[iii], LENGTH, fr_reads);
						if ((strlen(Read[iii])<550) && (strlen(Read[iii])>450) && (iii<99)){
								iii++;
						}
					}
				}
			   	fclose(fr_reads);
				COV = iii;
				ScoreCov = 0;
				if (COV > 3){
					M1_Score = 0;
					for (ii=0; ii<COV; ii++){
						rand3number(ii, &i2, &i3);
						Score[ii] = GetLOCAL(Read[ii], Read[i2], 0, 500, &temp_Read_start, &temp_Read_end);
						if (Score[ii] > 150)
							Score[ii] = Score[ii] + GetLOCAL(Read[ii], Read[i3], 0, 500, &temp_Read_start, &temp_Read_end);
						if (Score[ii] > M1_Score){
							Max_i1 = ii;
							Max_i2 = i2;
							Max_i3 = i3;
							M1_Score = Score[ii];
						}
					}
					SampleAny3(Read[Max_i1],Read[Max_i2],Read[Max_i3],Top3Read);
					//printf("\nNo[%d],[%d],[%d]  MaxScore = %d\n",Max_i1,Max_i2,Max_i3,M1_Score); 
					/*sprintf(Top3GenomeFile,"Top3_COV%d_N150_read%d_i%d.fasta",COV,file_i,file_j);
					fw_sample3 = fopen(Top3GenomeFile, "w");
					if (fw_sample3 == NULL){
					    		printf("write Genome file %d error!\n", file_i);
					      		getchar();
					}
					fprintf (fw_sample3,">read%d_i%d_%d\n", file_i,file_j,COV);
					fprintf (fw_sample3,"%s", Top3Read);
					fclose(fw_sample3);*/

					ScoreCov = 0;
					for (ii=0; ii<COV; ii++){
						Score[ii] = GetLOCAL(Top3Read, Read[ii], 0, 500, &temp_Read_start, &temp_Read_end);
						//printf("No[%d],%d\n",ii,Score[ii]); 
						if (Score[ii] > 300)
							ScoreCov++;
					}
					iii = 0;
					if (ScoreCov > 3){
						while(exit_tag < 16){
							rand4number(&i1, &i2, &i3, &i4,Score);
							//printf("Sample Readsfile[%d]: %d %d %d %d\n", iii, i1, i2, i3, i4);
							SampleAny4(Read[i1], Read[i2], Read[i3], Read[i4], Sample4[iii]);
							iii++;
							rand4number(&i1, &i2, &i3, &i4,Score);
							//printf("Sample Readsfile[%d]: %d %d %d %d\n", iii, i1, i2, i3, i4);
							SampleAny4(Read[i1], Read[i2], Read[i3], Read[i4], Sample4[iii]);
							iii++;
							rand4number(&i1, &i2, &i3, &i4,Score);
							//printf("Sample Readsfile[%d]: %d %d %d %d\n", iii, i1, i2, i3, i4);
							SampleAny4(Read[i1], Read[i2], Read[i3], Read[i4], Sample4[iii]);
							iii++;
							rand4number(&i1, &i2, &i3, &i4,Score);
							//printf("Sample Readsfile[%d]: %d %d %d %d\n", iii, i1, i2, i3, i4);
							SampleAny4(Read[i1], Read[i2], Read[i3], Read[i4], Sample4[iii]);
							iii++;
							SampleAny4(Sample4[iii-4], Sample4[iii-3], Sample4[iii-2], Sample4[iii-1], Sample4[iii]);
 							//printf("%d: %s\n", iii,Sample4[iii]);getchar();
							iii++;
							exit_tag++;
						}
						SampleAny4(Sample4[4], Sample4[9], Sample4[14], Sample4[19], Sample4[90]);
						SampleAny4(Sample4[24], Sample4[29], Sample4[34], Sample4[39], Sample4[91]);
						SampleAny4(Sample4[44], Sample4[49], Sample4[54], Sample4[59], Sample4[92]);
						SampleAny4(Sample4[64], Sample4[69], Sample4[74], Sample4[79], Sample4[93]);
						SampleAny4(Sample4[90], Sample4[91], Sample4[92], Sample4[93], Sample4[94]);

						if ((GetLOCAL(Sample4[94], Sample4[90], 0, 500, &temp_Read_start, &temp_Read_end) > strlen(Sample4[94])*0.9) && (GetLOCAL(Sample4[94], Sample4[91], 0, 500, &temp_Read_start, &temp_Read_end) > strlen(Sample4[94])*0.9) && (GetLOCAL(Sample4[94], Sample4[92], 0, 500, &temp_Read_start, &temp_Read_end) > strlen(Sample4[94])*0.9) && (GetLOCAL(Sample4[94], Sample4[93], 0, 500, &temp_Read_start, &temp_Read_end) > strlen(Sample4[94])*0.9)){				
							sprintf(CorrectedGenomeFile,"%s",argv[2]);
							fw_sample4 = fopen(CorrectedGenomeFile, "w");
							if (fw_sample4 == NULL){
						     		printf("write Genome file %d error!\n", file_i);
						       		getchar();
							}
							fprintf (fw_sample4,">corrected_read_Cov%d\n", ScoreCov);
							fprintf (fw_sample4,"%s", Sample4[94]);
							fclose(fw_sample4);
							printf("Corrected! (Effective Cov = %d) \n", ScoreCov);
						}else{
							printf("Not Corrected. (Effective Cov = %d)\n", ScoreCov);
						}
					}else{
						printf("Not Corrected. (Effective Cov = %d)\n", ScoreCov);
					}
				}else{
					printf("Not Corrected. (Effective Cov = %d)\n", ScoreCov);
				}
			}
		
	}
}


int rand4number(int* i1, int* i2, int* i3, int* i4, int* Score){
	(*i1) = (int)(rand()) % COV;
	(*i2) = (int)(rand()) % COV;
	(*i3) = (int)(rand()) % COV;
	(*i4) = (int)(rand()) % COV;
	while (Score[(*i1)] < 300){
		(*i1) = rand() % COV;
	}
	while ((*i2) == (*i1)||(Score[(*i2)] < 300)){
		(*i2) = rand() % COV;
	}
	while ( ((*i3) == (*i1)) || ((*i3) == (*i2))|| (Score[(*i3)] < 300) ){
		(*i3) = (int)(rand()) % COV;
	}
	while ( ((*i4) == (*i1)) || ((*i4) == (*i2)) || ((*i4) == (*i3)) || (Score[(*i4)] < 300) ){
		(*i4) = (int)(rand()) % COV;
	}
}

int rand3number(int i1, int* i2, int* i3){
	//(*i1) = (int)(rand()) % COV;
	(*i2) = (int)(rand()) % COV;
	(*i3) = (int)(rand()) % COV;
	while ((*i2) == i1){
		(*i2) = rand() % COV;
	}
	while ( ((*i3) == i1) || ((*i3) == (*i2)) ){
		(*i3) = (int)(rand()) % COV;
	}
}

int SampleAny4(char *a, char *b, char *c, char *d, char* s1){
	char temp_ab[LENGTH];
	char temp_cd[LENGTH];
	GetSCS(a,b, temp_ab);
	GetSCS(c,d, temp_cd);
	GetLCS(temp_ab, temp_cd, s1);
}

int SampleAny3(char *a, char *b, char *c, char *abc){
	char temp_ab[LENGTH];
	char temp_bc[LENGTH];
	char temp_ac[LENGTH];
	char temp1[LENGTH];
	GetSCS(a,b, temp_ab);
	GetSCS(b,c, temp_bc);
	GetSCS(c,a, temp_ac);
	GetLCS(temp_ab, temp_bc, temp1);
	GetLCS(temp1, temp_ac, abc);
}




int** Make2DIntArray(int arraySizeX, int arraySizeY) {
	int** theArray;
	int i; 
	theArray = (int**) malloc(arraySizeX*sizeof(int*));
	for (i = 0; i < arraySizeX; i++)    
		theArray[i] = (int*)malloc(arraySizeY*sizeof(int));
	return theArray; 
} 


/* 
 * Note that this will allocate space for the result, and you will need 
 * to free it yourself to avoid a memory leak. 
 */
int GetLCS(char* a, char* b, char* lcs) {
	int n = strlen(a);
	int m = strlen(b);

	//int S = Make2DIntArray(n+1,m+1);
	//int R = Make2DIntArray(n+1,m+1);

	int S[999][999];
	int R[999][999];

	int ii;
	int jj;

	int i;

	int pos;


	/* It is important to use <=, not <.  The next two for-loops are initialization */
	for(ii = 0; ii <= n; ++ii) {
		S[ii][0] = 0;
		R[ii][0] = UP;
	}
	for(jj = 0; jj <= m; ++jj) {
		S[0][jj] = 0;
		R[0][jj] = LEFT;
	}



	/* This is the main dynamic programming loop that computes the score and */
	/* backtracking arrays. */
	for(ii = 1; ii <= n; ++ii) {
		for(jj = 1; jj <= m; ++jj) { 

			if( a[ii-1] == b[jj-1] ) {
				S[ii][jj] = S[ii-1][jj-1] + 1;
				R[ii][jj] = UP_AND_LEFT;
			}

			else {
				S[ii][jj] = S[ii-1][jj-1] + 0;
				R[ii][jj] = NEITHER;
			}

			if( S[ii-1][jj] >= S[ii][jj] ) {	
				S[ii][jj] = S[ii-1][jj];
				R[ii][jj] = UP;
			}

			if( S[ii][jj-1] >= S[ii][jj] ) {
				S[ii][jj] = S[ii][jj-1];
				R[ii][jj] = LEFT;
			}
		}
	}

	/* The length of the longest substring is S[n][m] */
	ii = n; 
	jj = m;
	pos = S[ii][jj];
	lcs[pos--] = '\0';


	/* Trace the backtracking matrix. */
	while( ii > 0 || jj > 0 ) {
		if( R[ii][jj] == UP_AND_LEFT ) {
			ii--;
			jj--;
			lcs[pos--] = a[ii];
		}
	
		else if( R[ii][jj] == UP ) {
			ii--;
		}

		else if( R[ii][jj] == LEFT ) {
			jj--;
		}
	}


	/*for (i = 0; i < n+1; i++){    
		free(S[i]);
		free(R[i]); 
	}
	free(S);
	free(R);*/
}

int GetSCS(char* a, char* b, char* scs) {
	int n = strlen(a);
	int m = strlen(b);

	//int** S = Make2DIntArray(n+1,m+1);
	//int** R = Make2DIntArray(n+1,m+1);
	int S[999][999];
	int R[999][999];

	int ii;
	int jj;

	int i;

	int pos;
	int max_length;
	int last_pos;


	/* It is important to use <=, not <.  The next two for-loops are initialization */
	for(ii = 0; ii <= n; ++ii) {
		S[ii][0] = 0;
		R[ii][0] = UP;
	}
	for(jj = 0; jj <= m; ++jj) {
		S[0][jj] = 0;
		R[0][jj] = LEFT;
	}

	//printf("scs[%d] = %c\n", 0, scs[0]); getchar();

	/* This is the main dynamic programming loop that computes the score and */
	/* backtracking arrays. */
	pos = 0;
	for(ii = 1; ii <= n; ++ii) {
		for(jj = 1; jj <= m; ++jj) { 

			if( a[ii-1] == b[jj-1] ) {
				S[ii][jj] = S[ii-1][jj-1] + 1;
				R[ii][jj] = UP_AND_LEFT;
			}

			else {
				S[ii][jj] = S[ii-1][jj-1] + 0;
				R[ii][jj] = NEITHER;
			}

			if( S[ii-1][jj] >= S[ii][jj] ) {	
				S[ii][jj] = S[ii-1][jj];
				R[ii][jj] = UP;
			}

			if( S[ii][jj-1] >= S[ii][jj] ) {
				S[ii][jj] = S[ii][jj-1];
				R[ii][jj] = LEFT;
			}
		}
	}

	/* The length of the longest substring is S[n][m] */
	ii = n; 
	jj = m;
	//pos = S[ii][jj];
	pos = n + m - S[ii][jj];
	max_length = pos;

	scs[pos--] = '\0';
	for(i= 0; i < pos-1; i++){
		scs[i] = 'N';
	}
	//printf("scs[%d] = %c\n", pos, scs[pos]); getchar();

	int prev = 1;
	int test_ii;
	int test_jj;
	/* Trace the backtracking matrix. */
	while( ii > 0 || jj > 0 ) {
		if( R[ii][jj] == UP_AND_LEFT ) {
			ii--;
			jj--;
			scs[pos--] = a[ii];
			prev = prev + 1;
 			last_pos = pos+1;
		} else if( R[ii][jj] == UP ) {
			ii--;
			scs[pos--] = a[ii];
		} else if( R[ii][jj] == LEFT ) {
			jj--;
			scs[pos--] = b[jj];
		}
		//printf("scs[%d] = %c\n", pos+1, scs[pos+1]); getchar();
	}

	/*for (i = 0; i < n+1; i++){    
		free(S[i]);
		free(R[i]); 
	}
	free(S);
	free(R);*/
}


int GetLOCAL(char* sa, char* b, int a1, int a2, int *b1, int *b2) {
	int n = strlen(sa);

	int n1 = strlen(sa);
	int m = strlen(b);

	//printf("len(a) = %d len(b) = %d\n", n1, m); getchar(); 

	(*b1) = -1;
	(*b2) = -1;

	if ((m < 0.8*n)||(m>1.4*n)){
		//printf("Sequence B is too short!!!\n"); 
		return(0);
	}
	//int** S = Make2DIntArray(n+1,m+1);
	//int** R = Make2DIntArray(n+1,m+1);

	int S[999][999];
	int R[999][999];

	char a[1000];

	int ii;
	int jj;

	int i;

	int pos2;
	int pos1;

	int score = 0;


	/* It is important to use <=, not <.  The next two for-loops are initialization */
	for(ii = 0; ii <= n; ++ii) {
		a[ii] = sa[a1+ii];
		S[ii][0] = 0;
		R[ii][0] = NEITHER;
	}
	for(jj = 0; jj <= m; ++jj) {
		S[0][jj] = 0;
		R[0][jj] = NEITHER;
	}



	/* This is the main dynamic programming loop that computes the score and */
	/* backtracking arrays. */
	for(ii = 1; ii <= n; ++ii) {
		for(jj = 1; jj <= m; ++jj) { 
			S[ii][jj] = 0;
			R[ii][jj] = NEITHER;
			if( a[ii-1] == b[jj-1] ) {
				if (S[ii-1][jj-1] + 1 > S[ii][jj]){
					S[ii][jj] = S[ii-1][jj-1] + 1;
					R[ii][jj] = UP_AND_LEFT;
				}
			}else{
				if (S[ii-1][jj-1] - 1 > S[ii][jj]){
					S[ii][jj] = S[ii-1][jj-1] - 1;
					R[ii][jj] = UP_AND_LEFT;
				}
			}

			if( S[ii-1][jj] - 1 > S[ii][jj] ) {	
				S[ii][jj] = S[ii-1][jj] - 1;
				R[ii][jj] = UP;
			}

			if( S[ii][jj-1] -1 > S[ii][jj] ) {
				S[ii][jj] = S[ii][jj-1] - 1;
				R[ii][jj] = LEFT;
			}
		}
	}

	/* The length of the longest substring is S[n][m] */

	score = 0;
	for(jj = 0; jj <= m; ++jj) {
		//printf("Score = %d j = %d \n", S[n][jj], jj ); getchar(); 
		if (S[n][jj] > score){
			score = S[n][jj];
			pos2 = jj;	
		}
	}

	//printf("Score = %d j = %d \n", score, pos2); getchar(); 
	ii = n-1;
	jj = pos2-1;

	if (score > 0.89*n){
		/* Trace the backtracking matrix. */
		while( (ii > 0) &&  (jj > 0) && (R[ii][jj] != NEITHER) ) {
			if( R[ii][jj] == UP_AND_LEFT ) {
				ii--;
				jj--;
				//lcs[pos--] = a[ii];
			}	
			else if( R[ii][jj] == UP ) {
				ii--;
			}
			else if( R[ii][jj] == LEFT ) {
				jj--;
			}
		}
		(*b2) = pos2-1;
		(*b1) = jj;
	}

	/*for (i = 0; i < n+1; i++){    
		free(S[i]);
		free(R[i]); 
	}
	free(S);
	free(R);*/

	return(score);
}

