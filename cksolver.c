#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<complex.h>
#include<math.h>
#define SUCCESS 1
#define EMPTY_LINE -1
#define INVALID_NAME -2
#define INVALID_NODE -3
#define INVALID_NUMBER -4
#define INVALID_ELEMENT -5
#define EOL_JUNK -6 // For unexpected tokens
#define MAX_LENGTH 512
typedef struct element  
{       struct element *left,*right;
        char depname[15],name[15],loc[4][15];
        double value;
        int n1,n2,n3,n4,eqnid;
} element;
element* createNode ( char* name , double value , char * loc1 ,
 char * loc2 , char * loc3 , char * loc4 , char* depname) 
{
        element* node = malloc(sizeof(element));
        node->left=node->right=NULL;
        node->value=value;
        if ( name != NULL )
                strcpy(node->name,name);
        if ( depname != NULL )
                strcpy(node->depname,depname);
        if ( loc1 != NULL )
                strcpy(node->loc[0],loc1);
        if ( loc2 != NULL )
                strcpy(node->loc[1],loc2);
        if ( loc3 != NULL )
                strcpy(node->loc[2],loc3);
        if ( loc4 != NULL )
                strcpy(node->loc[3],loc4);
        return node;
}
typedef struct elemList 
{
        element *left , *right;
        int count;
}elemList;
elemList* initList () 
{
        elemList * list=malloc(sizeof(elemList));
        list->left=list->right=NULL;
        list->count=0;
        return list;
}
//Add to right of list
void addRight(elemList* list , element* node ) 
{
        node->left=list->right;
        node->right=NULL;
        if ( list -> right )
        {
                list->right->right=node;
                list->right=node;
        }
        else
        {
                list->left=node;
                list->right=node;
        }
        list->count++;
}
//Delete from left of list : 
void deleteLeft(elemList* list) 
{
        if(list->count)
        {
                if ( list->left->right )
                {
                        element * temp = list->left;
                        list->left->right->left=NULL;
                        list->left=list->left->right;
                        free(temp);
                }
                else
                {
                        free(list->left);
                        list->right=list->left=NULL;
                }
                list->count--;
        }
}
//Function to traverse a name list and find a key.
element* searcheList ( elemList *l , char * key )
{
    element * t = l->left;
    do
        {
                if ( strcmp(t->name,key)==0 )
                    return t;
                t=t->right;
        }while(t!=NULL);
    return NULL;
}
//Parser and other string Function definitions:
int checkName ( char * name )
{       
        int i=0;
        for(;name[i]!='\0';i++)
                if(!(isalnum(name[i])))
                        return 0;
        return 1;
}
int getNode ( char * name )
{
        char * temp = name+strlen(name)-1;
        if(!(isdigit(*temp)))
                return INVALID_NODE;
        while ( temp!=name && isdigit(*(temp-1)))
                temp--;
        int r=atoi(temp);
        if(r<0)
                return INVALID_NODE;
        return r;
} 
int getValue( char* data , double * val )
{
        //Find base number.
        if(sscanf(data,"%lf",val)!=1)
                return 0;
        int i=strlen(data)-1;
        while ( (isalpha(data[i-1])) && isalpha(data[i]))
                i--;
        // Find any multipliers and change value accordingly
        if (i<strlen(data) && isalpha(data[i]))
                if(strcmp((data+i),"meg")==0)
                        *val=(*val ) * 1.0e6;
                else if(strcmp((data+i),"k")==0)
                        *val=(*val ) * 1.0e3;
                else if(strcmp((data+i),"m")==0)
                        *val=(*val ) * 1.0e-3;
                else if(strcmp((data+i),"u")==0)
                        *val=(*val ) * 1.0e-6;
                else if(strcmp((data+i),"n")==0)
                        *val=(*val ) * 1.0e-9;
                else
                        return 0;
        return 1;
}
int parseElement(char * desc , element ** node)
{
        //Get name and location of element
        char * name = strtok(desc," \t\n");
        if(name==NULL)
                return EMPTY_LINE;
        char* data[5];
        if(!(checkName(name)))
                return INVALID_NAME;
        data[0]=strtok(NULL," \t\n");
        data[1]=strtok(NULL," \t\n");
        if(data[0]==NULL||data[1]==NULL||!
(checkName(data[0]) && checkName(data[1])))
                return INVALID_NODE;
        //Get type specific data of element
        switch (name[0])
        {
                case 'R': case 'r' :
                case 'L' : case 'l' :
                case 'C' : case 'c' :
                case 'I' : case 'i' :
                case 'V' : case 'v' :
                {
                        data[2]=strtok(NULL," \t\n");
                                double val;
                        if(data[2]==NULL)
                                return INVALID_NUMBER;
                        if(getValue(data[2],&val)==0)
                                return INVALID_NUMBER;
                        data[4]=strtok(NULL," \t\n");
                        if(data[4]!=NULL && data[4][0]!='#')
                                return EOL_JUNK;
                     *node=createNode(name,val,data[0],data[1],NULL,NULL,NULL);
                        return SUCCESS;
                }
                break;
                case 'E' : case 'e' :
                case 'G' : case 'g' :
                {
                        data[2]=strtok(NULL," \t\n");
                        data[3]=strtok(NULL," \t\n");
                        if(data[2]==NULL||data[3]==NULL||!
(checkName(data[2]) && checkName(data[3])))
                                return INVALID_NODE;
                        data[4]=strtok(NULL," \t\n");
                        double val;
                        if(data[4]==NULL)
                                return INVALID_NUMBER;
                        if(getValue(data[4],&val)==0)
                                return INVALID_NUMBER;
                        data[4]=strtok(NULL," \t\n");
                        if(data[4]!=NULL && data[4][0]!='#')
                                return EOL_JUNK;
                        *node=createNode(name,val,data[0],data[1],data[2],data[3],NULL);
                        return SUCCESS;
                }
                break;
                case 'H' : case 'h' :
                case 'F' : case 'f' :
                {
                        data[2]=strtok(NULL," \t\n");
                        if(data[2]==NULL || !(checkName(data[2])))
                                return INVALID_NAME;
                        if(data[2][0]!='V' && data[2][0]!='v')
                                return INVALID_NAME;
                        data[3]=strtok(NULL," \t\n");
                        double val;
                        if(data[3]==NULL)
                                return INVALID_NUMBER;
                        if(getValue(data[3],&val)==0)
                                return INVALID_NUMBER;
                        data[4]=strtok(NULL," \t\n");
                        if(data[4]!=NULL && data[4][0]!='#')
                                return EOL_JUNK;
                        *node=createNode(name,val,data[0],data[1],NULL,NULL,data[2]);
                        return SUCCESS;
                }
                break;                
                default :
                        return INVALID_ELEMENT;
        }
}
//List structure to contain node names :
typedef struct nname  
{       struct nname *left,*right;
        char name[15];
        int num;
} nname;
nname* createName ( char* name , int num ) 
{
        nname* node = malloc(sizeof(nname));
        node->left=node->right=NULL;
        node->num=num;
        if ( name != NULL )
                strcpy(node->name,name);
        return node;
}
typedef struct nameList 
{
        nname *left , *right;
        int count;
}nameList;
nameList* initnList () 
{
        nameList * list=malloc(sizeof(nameList));
        list->left=list->right=NULL;
        list->count=0;
        return list;
}
//Add to right of list
void addRightn(nameList* list , nname* node ) 
{
        node->left=list->right;
        node->right=NULL;
        if ( list -> right )
        {
                list->right->right=node;
                list->right=node;
        }
        else
        {
                list->left=node;
                list->right=node;
        }
        list->count++;
}
//Delete from left of list : 
void deleteLeftn(nameList* list) 
{
        if(list->count)
        {
                if ( list->left->right )
                {
                        nname * temp = list->left;
                        list->left->right->left=NULL;
                        list->left=list->left->right;
                        free(temp);
                }
                else
                {
                        free(list->left);
                        list->right=list->left=NULL;
                }
                list->count--;
        }
}
//Function to traverse a name list and find a key.
nname* searchnList ( nameList *l , char * key )
{
    nname * t = l->left;
    do
        {
                if ( strcmp(t->name,key)==0 )
                    return t;
                t=t->right;
        }while(t!=NULL);
        return NULL;
}
//Equation solving modules
//Row operation function
void rowop( double complex ** a , int n , int R1 , int R2 , double complex factor )
{
        int i=0;
        for(;i<n;i++)
                a[R1][i]+=a[R2][i]*factor;
}
#define EPSILON 1.0e-10
/*
Linear equation solver that uses Gaussian elimination
to solve a system of linear equations , with n equations
and variables. The (complex) coefficients in each equation are
provided in coeff[n][n] , and the constants on the right hand
side of each equations are provided in the column vector rhs[n].
*/
int gaussianSolver ( double complex ** coeff , double complex * rhs , int n )
{
        int i=0,j=0;
        //Outer iteration over first n-1 columns
        for(;i<n-1;i++)
        {
                j=i;
                //See if coeff[j][j] is non zero
                if(cabs(coeff[j][i])<EPSILON)
                {
                        //if not iterate over rows to find a non zero value
                        int ind=-1;
                        for(;j<n;j++)
                                if(cabs(coeff[j][i])>EPSILON)
                                {
                                        ind = j;
                                        break;
                                }
                        //If not found , terminate
                        if(ind==-1)
                                return 0;
                        //If found , add row ind to row i
                        rowop(coeff,n,i,ind,1.0);
                        rhs[i]+=rhs[ind];
                }
                //Reduce lower elements of the column to zero
                for(j=i+1;j<n;j++)
                        if(cabs(coeff[j][i])>EPSILON)
                        {
                                rhs[j]-=rhs[i]*(coeff[j][i]/coeff[i][i]);
                                rowop(coeff,n,j,i,-(coeff[j][i]/coeff[i][i]));
                        }
        }
        if(cabs(coeff[n-1][n-1])<EPSILON)
                return 0;
        //Solve by reverse substitution
        rhs[n-1]/=coeff[n-1][n-1];
        for(i=n-2;i>=0;i--)
        {
                for(j=n-1;j>i;j--)
                        rhs[i]-=rhs[j]*coeff[i][j];
                rhs[i]/=coeff[i][i];
        }
        return 1;
}
//Circuit solver module
int circuitSolve (nameList *nlist , elemList * list ,
int nnodes , int nvs , double freq , double complex ** result )
{
        //Initialize the two matrices
        double complex ** A = malloc((nnodes+nvs)*sizeof(double complex*));
        int i=0,j=0;
        for(;j<(nnodes+nvs);j++)
        {
                A[j]=malloc(nnodes*sizeof(double complex));
                int i=0;
                for(;i<nnodes;i++)
                        A[j][i]=0.0;
        }
        double complex ** B = malloc((nnodes+nvs+1)*sizeof(double complex*));
        j=0;
        for(;j<(nnodes+nvs);j++)
        {
                B[j]=malloc(nvs*sizeof(double complex));
                int i=0;
                for(;i<nvs;i++)
                        B[j][i]=0.0;
        }
        double complex * b = malloc((nnodes+nvs)*sizeof(double complex));
        j=0;
        for(;j<(nnodes+nvs);j++)
                b[j]=0.0;
        //Fill in matrices while traversing the element list
        element* t=list->left;
        int flag =1;
        do
        {
                switch(t->name[0])
                {
                        case 'R' : case 'r' :
                        {
                                double complex z = t->value;
                                //Add to voltage matrix
                                A[t->n1][t->n1]+=1.0/z;
                                A[t->n2][t->n2]+=1.0/z;
                                A[t->n1][t->n2]+=-1.0/z;
                                A[t->n2][t->n1]+=-1.0/z;
                        }
                        break;
                        case 'L' : case 'l' :
                        {
                                double complex z = t->value * freq * I;
                                //Add to voltage matrix
                                A[t->n1][t->n1]+=1.0/z;
                                A[t->n2][t->n2]+=1.0/z;
                                A[t->n1][t->n2]+=-1.0/z;
                                A[t->n2][t->n1]+=-1.0/z;
                        }
                        break;
                        case 'C' : case 'c' :
                        {
                                double complex z = -I/ ( freq *t->value );
                                //Add to voltage matrix
                                A[t->n1][t->n1]+=1.0/z;
                                A[t->n2][t->n2]+=1.0/z;
                                A[t->n1][t->n2]+=-1.0/z;
                                A[t->n2][t->n1]+=-1.0/z;
                        }
                        break;
                        case 'I' : case 'i' :
                        {
                                //Add to source vector
                                double complex z = t->value;
                                b[t->n1]+=z; b[t->n2]+=-z;
                        }
                        break;
                        case 'V' : case 'v' :
                        {
                                double complex z = t->value;
                                //Add to current matrix
                                B[t->n1][t->eqnid]=1.0;
                                B[t->n2][t->eqnid]=-1.0;
                                //Add to voltage matrix
                                A[nnodes+t->eqnid][t->n1]+=1.0;
                                A[nnodes+t->eqnid][t->n2]+=-1.0;
                                //Add to source vector
                                b[nnodes+t->eqnid]=z;
                        }
                        break;
                        case 'E' : case 'e' :
                        {
                                nname * c,*d;
                                c=searchnList(nlist,t->loc[2]);
                                d=searchnList(nlist,t->loc[3]);
                                //Check if dependent nodes exist
                                if((c==NULL && getNode(t->loc[2])!=0)
||(d==NULL && getNode(t->loc[3])!=0))
                                {
                                        printf("Invalid node name for VCVS %s\n"
,t->name);
                                        flag=0;
                                }
                                else
                                {
                                        if(c==NULL)
                                                c=nlist->left;
                                        if(d==NULL)
                                                d=nlist->left;
                                        double complex z = t->value;
                                        //Add to current matrix
                                        B[t->n1][t->eqnid]+=1.0;
                                        B[t->n2][t->eqnid]+=-1.0;
                                        //Add to voltage matrix
                                        A[nnodes+t->eqnid][t->n1]+=1.0;
                                        A[nnodes+t->eqnid][t->n2]+=-1.0;
                                        A[nnodes+t->eqnid][c->num]+=-z;
                                        A[nnodes+t->eqnid][d->num]+=z;
                                }
                        }
                        break;
                        case 'G' : case 'g' :
                        {
                                nname * c,*d;
                                c=searchnList(nlist,t->loc[2]);
                                d=searchnList(nlist,t->loc[3]);
                                //Check if dependent nodes exist
                                if((c==NULL && getNode(t->loc[2])!=0)||
(d==NULL && getNode(t->loc[3])!=0))
                                {
                                        printf("Invalid node name for VCCS %s\n",
t->name);
                                        flag=0;
                                }
                                else
                                {
                                        if(c==NULL)
                                                c=nlist->left;
                                        if(d==NULL)
                                                d=nlist->left;
                                        double complex z = t->value;
                                        //Add to voltage matrix
                                        A[t->n1][c->num]+=z;
                                        A[t->n1][d->num]+=-z;
                                        A[t->n2][c->num]+=-z;
                                        A[t->n2][d->num]+=z;
                                }
                        }
                        break;
                        case 'F' : case 'f' :
                        {
                                element * dep = searcheList(list,t->depname);
                                if(dep==NULL||
(dep->name[0]!='V' && dep->name[0]!='v'))
                                {
                                        printf("Non-existent dependency for \
CCCS %s\n",t->name);
                                        flag=0;
                                }
                                else
                                {
                                        double complex z = t->value;
                                        //Add to current matrix
                                        B[t->n1][dep->eqnid]+=z;
                                        B[t->n2][dep->eqnid]+=-z;
                                }
                        }
                        break;
                        case 'H' : case 'h' :
                        {
                                element * dep = searcheList(list,t->depname);
                                if(dep==NULL||
(dep->name[0]!='V' && dep->name[0]!='v'))
                                {
                                        printf("Non-existent dependency for \
CCVS %s\n",t->name);
                                        flag=0;
                                }
                                else
                                {
                                        double complex z = t->value;
                                        //Add to current matrix
                                        B[t->n1][t->eqnid]+=1.0;
                                        B[t->n2][t->eqnid]+=-1.0;
                                        B[nnodes+t->eqnid][dep->eqnid]+=-z;
                                        //Add to voltage matrix
                                        A[nnodes+t->eqnid][t->n1]+=1.0;
                                        A[nnodes+t->eqnid][t->n2]+=-1.0;
                                }
                        }
                        break;
                }               
                if(flag==0)
                {
                        int i=0;
                        for(i=0;i<nnodes+nvs;i++)
                        {
                                free(A[i]);
                                free(B[i]);
                        }
                        free(b);
                        return 0;
                }
                t=t->right;
        }while(t!=NULL);
        //Ground node equation
        A[0][0]+=1.0;
        //Generate MNA matrix
        double complex ** MAT =
malloc((nnodes+nvs)*sizeof(double complex*));
        for(i=0;i<nnodes+nvs;i++)
        {
                MAT[i]=malloc((nnodes+nvs)*sizeof(double complex));
                int j=0;
                for(;j<nnodes;j++)
                        MAT[i][j]=A[i][j];
                for(;j-nnodes<nvs;j++)
                        MAT[i][j]=B[i][j-nnodes];
        }
        //Try to solve equations
        int s = gaussianSolver(MAT,b,nnodes+nvs);
        if(s==0)
        {
                printf("Equations unsolvable. Inconsistent circuit.\n");
                *result = NULL;
                return 0;
        }
        //If successful , store pointer to solution vector in parameter
        *result=b;
        //Free other intermediate memory
        i=0;
        for(i=0;i<nnodes+nvs;i++)
        {
                free(A[i]);
                free(B[i]);
        }
        free(A);free(B);
        return 1;
}
//Simulation functions
int DCSim ( elemList * list , nameList * nlist , int nnodes , int nvs ,
 element * var , double freq , double vi , double vf , double step , FILE * fout)
{
        if(step<0 && vf>vi)
        {
                printf("Invalid limits and step\n");
                return 0;
        }
        var->value=vi;
        do
        {
                double complex * result;
                //call solver
                int r = circuitSolve ( nlist , list , nnodes , nvs , freq , &result);
                if (r==0)
                        return 0;
                //Write to file and free memory
                int i=0;
                fprintf(fout,"%.9f ",var->value);
                for(;i<nnodes+nvs;i++)
                        fprintf(fout,"%.9f ",creal(result[i]));
                fprintf(fout,"\n");
                free(result);
                //increment value
                var->value+=step;
        }while (fabs(step)>EPSILON && 
( (step<0 && var->value>=vf) || (step>0 && var->value<=vf) ) );
        return 1;
}
int ACSim ( elemList * list , nameList * nlist , int nnodes , int nvs ,
int n1 , int n2 , double fi , double ff , int ndec , FILE * fout)
 {
        if(fi>ff)
        {
                printf("Upper limit cannot be less than lower limit\n");
                return 0;
        }
        if(ndec<=0)
        {
                printf("Invalid values/decade\n");
                return 0;
        }
        int nfreqs = round(log10(ff/fi)*(double)ndec) + 1;
        double lgfreq = log10(fi);//logarithm of current frequency
        do
        {
                double complex * result;
                //call solver
                int r = circuitSolve ( nlist , list ,
nnodes , nvs , pow(10.0,lgfreq) , &result);
                if (r==0)
                        return 0;
                //Write to file and free memory
                int i=0;
                fprintf(fout,"%.9f ",pow(10.0,lgfreq)/(2*M_PI));
                for(;i<nnodes+nvs;i++)
                        fprintf(fout,"%.9f %.9f ",cabs(result[i]),carg(result[i]));
                fprintf(fout,"\n");
                free(result);
                //Logarithmically increment simulation frequency
                lgfreq=lgfreq+1.0/(double)ndec;
                nfreqs--;
        }while(nfreqs);
        //Open plotting thread and plot log-log plot of magnitude :
        FILE * gp = popen ("gnuplot -persistent", "w");
        fprintf(gp, "set logscale xy; plot \"spice.out\" \
        using 1:(abs($%d*exp({0,1}*$%d)-$%d*exp({0,1}*$%d))) \
        title 'Amplitude' linecolor rgb \"blue\" smooth unique\n"
        ,2*n1+2,2*n1+3,2*n2+2,2*n2+3);
        fclose(gp);
        return 1;
 }
int main( int argv ,char ** argc ) 
{
        if(argv != 2)
        {
                printf("Usage ./a.out <filename>\n");
                exit(1);
        }
        FILE *fp = fopen(argc[1], "r");
        if(fp == NULL)
        {
                printf("File could not be opened\n");
                exit(2);
        }
        elemList* list = initList();
        char buf[MAX_LENGTH];
        int found = 0 , line = 0;
        while(fgets(buf, MAX_LENGTH, fp))
        {
                line++;
                if(strcmp(buf,".circuit\n")==0)
                {
                        found=1;
                        break;
                }
        }
        if ( found==0 )
        {
                printf("No circuit description found. Exiting. \n");
                return 0;
        }
        while(fgets(buf, MAX_LENGTH, fp))
        {
                char * temp = strtok(buf,"#\n");
                if(temp != NULL)
                {
                strcpy(buf,temp);
                if(buf!=NULL && (strcmp(buf,".end\n")==0||strcmp(buf,".end")==0) )
                        break;
                element * t;
                int r = parseElement(buf,&t);
                if(r!=SUCCESS)
                {
                        switch(r)
                        {
                                case INVALID_NAME :
                                        printf("Invalid name on line %d\n",line );
                                break;
                                case INVALID_NODE :
                                        printf("Invalid node on line %d\n",line );
                                break;
                                case INVALID_NUMBER :
                                        printf("Invalid value on line %d\n",line );
                                break;
                                case INVALID_ELEMENT :
                                printf("Invalid element type on line %d\n",line );
                                break;
                                case EOL_JUNK :
                                        printf("Too many tokens on line %d\n",line );
                                break;
                        }
                        while(list->count)
                                deleteLeft(list);
                        free(list);
                        fclose(fp);
                        return 0;
                }
                else 
                        addRight(list,t);
                }
                line++;
        }
    int gnfound = 0;
    nameList * nlist = initnList();
    addRightn(nlist , createName("ground",0) );
    element* t = list->left;
    /*Allocate node numbers to names and store node numbers
      corresponding to an element's location.
    */
    do
        {
                if(getNode(t->loc[0])==0 || strcmp("ground",t->loc[0])==0)
                {
                    gnfound=1;
                    t->n1=0;
                }
                else
                {
                    nname * temp = searchnList(nlist,t->loc[0]);
                    if(temp==NULL)
                    {
                        t->n1=nlist->count;
                        addRightn(nlist,createName(t->loc[0],nlist->count));
                    }
                    else
                        t->n1=temp->num;
                }
                if(getNode(t->loc[1])==0 || strcmp("ground",t->loc[1])==0)
                {
                    gnfound=1;
                    t->n2=0;
                }
                else
                {
                    nname * temp = searchnList(nlist,t->loc[1]);
                    if(temp==NULL)
                    {
                        t->n2=nlist->count;
                        addRightn(nlist,createName(t->loc[1],nlist->count));
                    }
                    else
                        t->n2=temp->num;
                }
                t=t->right;
        }while(t!=NULL);
        if(gnfound==0)
        {
        printf("No reference node found. Exiting\n");
        while(list->count)
                    deleteLeft(list);
            free(list);
            while(nlist->count)
                    deleteLeftn(nlist);
            free(nlist);
            fclose(fp);
            return 0;
    }
        //See how many nodes are to be added for voltage sources.
    int nnodes=nlist->count,nvs=0;    
        t=list->left;
    do
        {
        if(t->name[0]=='V' || t->name[0]=='v' || t->name[0]=='E'
        ||t->name[0]=='e'||t->name[0]=='H'||t->name[0]=='h')
        {
            t->eqnid=nvs;
            nvs++;
        }
        t=t->right;
        }while(t!=NULL);
        //Print node names and numbers
        printf("Printing node names and numbers:\n");
        //Print node data
        nname * f = nlist->left;
        do
        {
                printf("Name:%s  Number:%d\n",f->name,f->num);
                f=f->right;
        }while(f!=NULL);
        //Print Voltage source equation numbers
        t=list->left;
    do
        {
        if(t->name[0]=='V' || t->name[0]=='v' || t->name[0]=='E'
        ||t->name[0]=='e'||t->name[0]=='H'||t->name[0]=='h')
            printf("VS Name:%s  Eqn. no. :%d\n",t->name,t->eqnid+nnodes);
        t=t->right;
        }while(t!=NULL); printf("\n");
        //File for output storage
        FILE * fout = fopen("spice.out","w");
        // Check for simulation instructions
        found=0;
        while(fgets(buf, MAX_LENGTH, fp))
        {
                line++;
                if(strcmp(buf,".command\n")==0)
                {
                        found=1;
                        break;
                }
        }
        if ( found==0 )
        {
                printf("No simulation instructions found. Exiting. \n");
                return 0;
        }
        //Accept simulation instruction from user
        while(fgets(buf, MAX_LENGTH, fp))
        {
                char * sim ;
                sim=strtok(buf,"#\n");
                if(sim!=NULL)
                {
                if(strcmp(sim,".end\n")==0||strcmp(sim,".end")==0)
                        break;
                char * tok;
                tok=strtok(sim," \n");
                if (strcmp(tok,"dc")==0)
                {
                        tok=strtok(NULL," \n");
                        element* var;
                        if(tok==NULL ||  (var=searcheList(list,tok))==NULL )
                        {
                                printf("Element to vary does not exist\n");
                                break;
                        }
                        else
                        {
                                //Accept and validate inputs
                                tok=strtok(NULL," \n");
                                double vinit,vfin,vstep; int flag=1;
                                if(tok==NULL || getValue(tok,&vinit)==0)
                                {
                                        printf("Invalid value\n");
                                        break;
                                }
                                tok=strtok(NULL," \n");
                                if(tok==NULL || getValue(tok,&vfin)==0)
                                {
                                        printf("Invalid value\n");
                                        break;
                                }
                                tok=strtok(NULL," \n");
                                if(tok==NULL || getValue(tok,&vstep)==0)
                                {
                                        printf("Invalid value\n");
                                        break;
                                }
                                //Call simulator
                                DCSim(list,nlist,nnodes,nvs,
var,1.0,vinit,vfin,vstep,fout);
                                break;
                        }
                }
                else if (strcmp(tok,"ac")==0)
                {
                        tok=strtok(NULL," \n");
                        char * tok1 = strtok (NULL," \n");
                        nname *n1=NULL , *n2=NULL ;
                        int node1=0 , node2=0;
                        if(tok && tok1)
                        {
                                n1=searchnList(nlist,tok);
                                n2=searchnList(nlist,tok1);
                        }
                        if((n1==NULL && getNode(tok)!=0)
||(n2==NULL && getNode(tok1)!=0))
                        {
                                printf("Invalid node\n");
                                break;
                        }
                        if(getNode(tok))
                                node1=n1->num;
                        if(getNode(tok1))
                                node2=n2->num;
                        double finit , ffin ; int ndec;
                        tok=strtok(NULL," \n");
                        tok1=strtok(NULL," \n");
                        if((tok==NULL || getValue(tok,&finit)==0) 
|| (tok1==NULL || getValue(tok1,&ffin)==0)
                        || ( (tok=strtok(NULL," \n"))==NULL 
|| (ndec=getNode(tok))<=0 ) )
                        {
                                printf("Invalid value\n");
                                break;
                        }
                        ACSim(list,nlist,nnodes,nvs,
node1,node2,2*M_PI*finit,2*M_PI*ffin,ndec,fout);
                        break;
                }
                else
                {
                        printf("Invalid command\n");
                        break;
                }
                }
                line++;
        }
        //Deallocate memory , close file and terminate
        while(list->count)
                deleteLeft(list);
        free(list);
        while(nlist->count)
                deleteLeftn(nlist);
        free(nlist);
        fclose(fp);
        fclose(fout);
        return 0;
}
