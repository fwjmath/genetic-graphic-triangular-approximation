#include <stdio.h>
#include <stdlib.h>
#include <io.h>

const int GeneNum = 100;
const int CommutateNum = 10;
const int Population = 30;
const int Mutation_Rate_Mask = 2047;

struct pixel{
	int R;
	int G;
	int B;
	int m;
};

struct triangle{
	int point[3][2];
	int R;
	int G;
	int B;
};

struct Amoeba{
	triangle gene[GeneNum];
	pixel appearance[200][200];
	int evaluation;
};




////////////////////////////////////////
pixel pic[200][200];
pixel pixel_average;
int height, width, iternum;
Amoeba pool[Population];
Amoeba best_now;

void init_pic(){
	FILE* fio;
	int i,j,k,r,g,b;
	fio=fopen("data.txt","r");
	fscanf(fio,"%d%d", &height, &width);
	for(i=0;i<height;i++){
		for(j=0;j<width;j++){
			fscanf(fio,"%d%d%d",&r,&g,&b);
			pic[i][j].R=r;
			pic[i][j].G=g;
			pic[i][j].B=b;
			pic[i][j].m=1;
		}
	}
	fclose(fio);
	pixel_average=pic[0][0];
	return;
}

int sgn_int(int a){
	if(a>0) return 1; else if(a<0) return -1; else return 0;
}

void cover_triangle(Amoeba* obj, int gen){
	int i, j, k, l, p, maxx, minx, maxy, miny, sgn1, sgn2;
	triangle tri;
	tri=obj->gene[gen];
	maxx=0;
	minx=height;
	maxy=0;
	miny=width;
	for(i=0;i<3;i++){
		if(tri.point[i][0]>maxx) maxx=tri.point[i][0];
		if(tri.point[i][0]<minx) minx=tri.point[i][0];
		if(tri.point[i][1]>maxy) maxy=tri.point[i][1];
		if(tri.point[i][1]<miny) miny=tri.point[i][1];
	}
	for(p=0;p<3;p++){
		l=p+1;
		if(l==3) l=0;
		sgn1=tri.point[l][0]-tri.point[p][0];
		sgn2=tri.point[l][1]-tri.point[p][1];
		if(sgn_int(sgn1)==0){
			i=tri.point[p][0];
			j=tri.point[p][1];
			sgn2=sgn_int(sgn2);
			do{
				obj->appearance[i][j].m=-abs(obj->appearance[i][j].m);
				j+=sgn2;
			}while(j!=tri.point[l][1]);
		}
		else if(sgn_int(sgn2)==0){
			i=tri.point[p][0];
			j=tri.point[p][1];
			sgn1=sgn_int(sgn1);
			do{
				obj->appearance[i][j].m=-abs(obj->appearance[i][j].m);
				i+=sgn1;
			}while(i!=tri.point[l][0]);
		}
		else if(abs(sgn1)>abs(sgn2)){
			i=tri.point[p][0];
			j=tri.point[p][1];
			k=0;
			do{
				obj->appearance[i][j].m=-abs(obj->appearance[i][j].m);
				i+=sgn_int(sgn1);
				k+=abs(sgn2);
				if((abs(k)<<1)>abs(sgn1)){
					j+=sgn_int(sgn2);
					k-=abs(sgn1);
				}
			}while((i!=tri.point[l][0])||(j!=tri.point[l][1]));
			obj->appearance[i][j].m=-abs(obj->appearance[i][j].m);
		}
		else{
			i=tri.point[p][0];
			j=tri.point[p][1];
			k=0;
			do{
				obj->appearance[i][j].m=-abs(obj->appearance[i][j].m);
				j+=sgn_int(sgn2);
				k+=abs(sgn1);
				if((abs(k)<<1)>abs(sgn2)){
					i+=sgn_int(sgn1);
					k-=abs(sgn2);
				}
			}while((i!=tri.point[l][0])||(j!=tri.point[l][1]));
			obj->appearance[i][j].m=-abs(obj->appearance[i][j].m);
		}
	}
	for(i=minx;i<=maxx;i++){
		k=miny;
		while(obj->appearance[i][k].m>=0) k++;
		l=maxy;
		while(obj->appearance[i][l].m>=0) l--;
		for(j=k;j<=l;j++){
			obj->appearance[i][j].m=1+abs(obj->appearance[i][j].m); 
			obj->appearance[i][j].R+=tri.R;
			obj->appearance[i][j].G+=tri.G;
			obj->appearance[i][j].B+=tri.B;
		}
	}
	return;
}

void delete_triangle(Amoeba* obj, int gen){
	int i, j, k, l, p, maxx, minx, maxy, miny, sgn1, sgn2;
	triangle tri;
	tri=obj->gene[gen];
	maxx=0;
	minx=height;
	maxy=0;
	miny=width;
	for(i=0;i<3;i++){
		if(tri.point[i][0]>maxx) maxx=tri.point[i][0];
		if(tri.point[i][0]<minx) minx=tri.point[i][0];
		if(tri.point[i][1]>maxy) maxy=tri.point[i][1];
		if(tri.point[i][1]<miny) miny=tri.point[i][1];
	}
	for(p=0;p<3;p++){
		l=p+1;
		if(l==3) l=0;
		sgn1=tri.point[l][0]-tri.point[p][0];
		sgn2=tri.point[l][1]-tri.point[p][1];
		if(sgn_int(sgn1)==0){
			i=tri.point[p][0];
			j=tri.point[p][1];
			sgn2=sgn_int(sgn2);
			do{
				obj->appearance[i][j].m=-abs(obj->appearance[i][j].m);
				j+=sgn2;
			}while(j!=tri.point[l][1]);
		}
		else if(sgn_int(sgn2)==0){
			i=tri.point[p][0];
			j=tri.point[p][1];
			sgn1=sgn_int(sgn1);
			do{
				obj->appearance[i][j].m=-abs(obj->appearance[i][j].m);
				i+=sgn1;
			}while(i!=tri.point[l][0]);
		}
		else if(abs(sgn1)>abs(sgn2)){
			i=tri.point[p][0];
			j=tri.point[p][1];
			k=0;
			do{
				obj->appearance[i][j].m=-abs(obj->appearance[i][j].m);
				i+=sgn_int(sgn1);
				k+=abs(sgn2);
				if((abs(k)<<1)>abs(sgn1)){
					j+=sgn_int(sgn2);
					k-=abs(sgn1);
				}
			}while((i!=tri.point[l][0])||(j!=tri.point[l][1]));
			obj->appearance[i][j].m=-abs(obj->appearance[i][j].m);
		}
		else{
			i=tri.point[p][0];
			j=tri.point[p][1];
			k=0;
			do{
				obj->appearance[i][j].m=-abs(obj->appearance[i][j].m);
				j+=sgn_int(sgn2);
				k+=abs(sgn1);
				if((abs(k)<<1)>abs(sgn2)){
					i+=sgn_int(sgn1);
					k-=abs(sgn2);
				}
			}while((i!=tri.point[l][0])||(j!=tri.point[l][1]));
			obj->appearance[i][j].m=-abs(obj->appearance[i][j].m);
		}
	}
	for(i=minx;i<=maxx;i++){
		k=miny;
		while(obj->appearance[i][k].m>=0) k++;
		l=maxy;
		while(obj->appearance[i][l].m>=0) l--;
		for(j=k;j<=l;j++){
			obj->appearance[i][j].m=abs(obj->appearance[i][j].m)-1; 
			obj->appearance[i][j].R-=tri.R;
			obj->appearance[i][j].G-=tri.G;
			obj->appearance[i][j].B-=tri.B;
		}
	}
	return;
}
void evaluate(Amoeba* obj){
	int i, j, s, m1;
	s=0;
	for(i=0;i<height;i++){
		for(j=0;j<width;j++){
			m1=obj->appearance[i][j].m;
			s+=abs((obj->appearance[i][j].R)/m1-pic[i][j].R)+abs((obj->appearance[i][j].G)/m1-pic[i][j].G)+abs((obj->appearance[i][j].B)/m1-pic[i][j].B);
		}
	}
	obj->evaluation=s;
	return;
}

void init_pool(){
	int i, j, k, l;
	__asm RDTSC;
	__asm xor eax, edx;
	__asm mov i, eax;
	srand((unsigned int)i);
	for(i=0;i<Population;i++){
		for(j=0;j<height;j++){
			for(k=0;k<width;k++){
				pool[i].appearance[j][k]=pixel_average;
			}
		}
		for(j=0;j<GeneNum;j++){
			for(k=0;k<3;k++){
				pool[i].gene[j].point[k][0]=rand()%height;
				pool[i].gene[j].point[k][1]=rand()%width;
			}
			pool[i].gene[j].R=rand()&255;
			pool[i].gene[j].G=rand()&255;
			pool[i].gene[j].B=rand()&255;
			cover_triangle(&pool[i],j);
		}
		evaluate(&pool[i]);
	}
	return;
}

void print_best(){
	int i, j, k, l, r, g, b, x;
	FILE* fio;
	fio=fopen("result.txt","a");
	fprintf(fio,"----------------------------------\n");
	fprintf(fio,"%d %d %d %d %d\n",height,width,GeneNum,iternum,best_now.evaluation);
	fprintf(fio,"%d %d %d\n",pixel_average.R,pixel_average.G,pixel_average.B);
	for(i=0;i<GeneNum;i++){
		fprintf(fio,"%d %d %d %d %d %d %d %d %d\n", best_now.gene[i].point[0][0],best_now.gene[i].point[0][1],best_now.gene[i].point[1][0],best_now.gene[i].point[1][1],best_now.gene[i].point[2][0],best_now.gene[i].point[2][1],best_now.gene[i].R,best_now.gene[i].G,best_now.gene[i].B);
	}
	/*for(i=0;i<height;i++){
		for(j=0;j<width;j++){
			r=best_now.appearance[i][j].R;
			x=(best_now.appearance[i][j].m);
			r/=x;
			g=best_now.appearance[i][j].G;
			g/=x;
			b=best_now.appearance[i][j].B;
			b/=x;
			fprintf(fio,"%d %d %d\n",r,g,b);
		}
	}*/
	fclose(fio);
	return;
}

void print_info(){
	printf("Generation: %d. Best evaluation: %d\n", iternum, best_now.evaluation);
	return;
}

void iterate_generation(){
	int i, j, k, l;
	int tri[Population][2];
	//sort by evaluation
	for(i=0;i<Population;i++){
		tri[i][0]=pool[i].evaluation;
		tri[i][1]=i;
	}
	for(i=0;i<Population;i++){
		for(j=0;j<Population-1;j++){
			if(tri[j][0]>tri[j+1][0]){
				l=tri[j][0];
				tri[j][0]=tri[j+1][0];
				tri[j+1][0]=l;
				l=tri[j][1];
				tri[j][1]=tri[j+1][1];
				tri[j+1][1]=l;
			}
		}
	}
	//Update Best Amoeba
	if(best_now.evaluation>pool[tri[0][1]].evaluation) best_now=pool[tri[0][1]];
	//Replace "failure" with "success", with commutation
	for(i=0;i<Population/3;i++){
		pool[tri[Population-1-i][1]]=pool[tri[i][1]];
		j=rand()%(Population/3);
		j+=Population/3;
		for(k=0;k<CommutateNum;k++){
			l=rand()%GeneNum;
			delete_triangle(&pool[tri[Population-1-i][1]],l);
			pool[tri[Population-1-i][1]].gene[l]=pool[j].gene[l];
			cover_triangle(&pool[tri[Population-1-i][1]],l);
		}
		evaluate(&pool[tri[Population-1-i][1]]);
	}
	//Mutation
	for(i=0;i<Population;i++){
		l=0;
		for(j=0;j<GeneNum;j++){
			for(k=0;k<9;k++){
				if((rand()&Mutation_Rate_Mask)==14){
					delete_triangle(&pool[i],j);
					if(k<6){
						if((k&1)==0){
							pool[i].gene[j].point[(k>>1)%3][0]+=(rand()%(height>>1))-(height>>2);
							if(pool[i].gene[j].point[(k>>1)%3][0]<0) pool[i].gene[j].point[(k>>1)%3][0]+=height;
							if(pool[i].gene[j].point[(k>>1)%3][0]>=height) pool[i].gene[j].point[(k>>1)%3][0]-=height;
						}
						else{
							pool[i].gene[j].point[(k>>1)%3][1]+=(rand()%(width>>1))-(width>>2);
							if(pool[i].gene[j].point[(k>>1)%3][1]<0) pool[i].gene[j].point[(k>>1)%3][1]+=width;
							if(pool[i].gene[j].point[(k>>1)%3][1]>=width) pool[i].gene[j].point[(k>>1)%3][1]-=width;
						}
					}
					else{
						if(k==6){
							pool[i].gene[j].R+=rand()%255;
							pool[i].gene[j].R&=255;
						}
						else if(k==7){
							pool[i].gene[j].G+=rand()%255;
							pool[i].gene[j].G&=255;
						}
						else{
							pool[i].gene[j].B+=rand()%255;
							pool[i].gene[j].B&=255;
						}
					}
					cover_triangle(&pool[i],j);
					l=1;
				}
			}
		}
		if(l==1) evaluate(&pool[i]);
	}
	iternum++;
	return;
}

bool CheckStop(){
	return access("stop.txt",0)==0;
}

int main(){
	int i;
	char c;
	init_pic();
	init_pool();
	best_now=pool[0];
	iternum=0;
	while(true){
		/*for(i=0;i<2000;i++){
			iterate_generation();
		}*/
		scanf("%d",&i);
		while(i>=iternum) iterate_generation();
		print_info();
		print_best();
	}
	return 0;
}