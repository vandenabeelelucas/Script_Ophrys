#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <sys/stat.h>

#define MAX(x,y) ( ( x > y ) ? x : y )
#define MIN(x,y) ( ( x > y ) ? y : x )
#define MAXFILEALI		100000
#define LENIDENT		80
#define LENIDENTPHYLIP	11   /* length of species names for PHYLIP */
#define LENSEQ			5000000
#define MAXSPECIES		5000
#define MAXCARSEQ		500000000
#define LENBUFFER		5000000
#define LENMOLDOM		100
#define LENFILENAME		300
#define MAXCLADE		1000
#define MAXSPECPERCLADE	200

/* list of input files */
char *buffer;
int nb_file_ali;
char file_ALI[MAXFILEALI][LENFILENAME];
char molecule[LENMOLDOM], domaine[LENMOLDOM], car_inconnu[100];
int type_seq;
int nb_spec_ali;
int nb_nuc_tot_ali;
long lenseq_tot;           /* Maximum length that sequences can have */
                           /* due to MAXCARSEQ and nb_spec_ali */
char seq[MAXCARSEQ];
char *ident;
int mask[MAXSPECIES], new[MAXSPECIES];

int nb_spec_uniq, nb_ident_uniq[MAXSPECIES];
char ident_uniq[MAXSPECIES][LENIDENT];
int min_spec, min_clade, nb_clade, nb_needed_clade, max_para, nb_file_OK;
char file_clade[100], file_needed_clade[100];
char name_clade[MAXCLADE][LENIDENT];
int nb_species_clade[MAXCLADE];
char species_clade[MAXCLADE][MAXSPECPERCLADE][LENIDENT];
int no_needed_clade[MAXCLADE], min_spec_needed_clade[MAXCLADE];
FILE *out_STAT;

static char Mois[12][10] = { "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"};
static char Jour[7][10] = { "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"};

void libere_memoire(void);
char read_line(FILE *in, int type, int *nb_line);
char interp_seq(int nb_line);
void date(char *version);
void date_fichier(char *nom, char *version);
FILE *ffopen(char *File, char *Mode);

/******************************************************************************/
/* Stores the current date in the variable 'version'                 */
/******************************************************************************/
void date(char *version)
{
struct tm *date;
time_t l_time;

time(&l_time);
date = localtime(&l_time);
sprintf(version,"%s %d %s %d at %d hours %d",Jour[date->tm_wday],date->tm_mday,Mois[date->tm_mon],1900 + date->tm_year,date->tm_hour,date->tm_min);
}

/******************************************************************************/
/* Stores the date of use of the file 'nom' in the variable 'version'         */
/******************************************************************************/
void date_fichier(char *nom, char *version)
{
struct tm *date_f;
struct stat buf;

stat(nom,&buf);
date_f = localtime(&buf.st_atime);
sprintf(version,"%d %s %d at %d hours %d",date_f->tm_mday,Mois[date_f->tm_mon],1900 + date_f->tm_year,date_f->tm_hour,date_f->tm_min);
}

/******************************************************************************/
/* Opens a file and verifies that it opened successfully                   */
/******************************************************************************/
FILE *ffopen(char *File, char *Mode)
{
FILE *in;

if ( (in = fopen(File,Mode)) == NULL )
	{
	fprintf(stderr,"Problem in opening file %s\n\n",File);
	exit(1);
	}
else
	return(in);
}

/******************************************************************************/
void strupr(char *chaine)
{
char *pc;

for ( pc = chaine ; *pc != 0 ; ++pc )
    *pc = toupper(*pc);
}

/******************************************************************************/
/* Reads the list of CARAC_FILE type files                    */
/******************************************************************************/
void read_name_file_ali()
{
FILE *in;

sprintf(buffer,"ls -1 *.ali >ali2paup.tmp");
system(buffer);

nb_file_ali = 0;
in = ffopen("ali2paup.tmp","r");
while ( fgets(buffer,LENBUFFER - 1,in) != NULL )
	{
	if ( nb_file_ali >= MAXFILEALI )
		{
		fprintf(stderr,"There are too many files for this program (MAX=%d)\n",MAXFILEALI);
		exit(1);
		}
	*strchr(buffer,'\n') = 0;
	strncpy(file_ALI[nb_file_ali++],buffer,LENFILENAME);
	}
fclose(in);
system("\"rm\" ali2paup.tmp");

if ( nb_file_ali == 0 )
	{
	fprintf(stderr,"There is no file for this program\n");
	exit(1);
	}
}

/******************************************************************************/
/* Frees dynamically allocated memory                       */
/******************************************************************************/
void libere_memoire()
{
if ( ident != NULL )
	{
	free(ident);
	ident = NULL;
	}
}

/******************************************************************************/
/* Reads a line in an .ALI file and discards comments  */
/******************************************************************************/
char read_line(FILE *in, int type, int *nb_line)
{
char *pc;

while ( fgets(buffer,LENBUFFER - 1,in) != NULL )
	{
	++*nb_line;
	if ( (pc = strchr(buffer,'\n')) == NULL )
		{
		sprintf(buffer,"Line %d is too long (MAXIMUM = %d)",*nb_line,LENSEQ - 1);
		return(-1);
		}
	if ( *buffer != '#' )
		{
		if ( (type == 0) && (*buffer != '>') )
			{
			sprintf(buffer,"Line %d is not properly formatted: need for >Identifier",*nb_line);
			return(-1);
			}
		*pc = 0;
		return(0);
		}
	}
return(-2);
}

/******************************************************************************/
/* Interprets a line containing a sequence                    */
/******************************************************************************/
char interp_seq(int nb_line)
{
char *pc;

strupr(buffer);
nb_nuc_tot_ali = MAX(nb_nuc_tot_ali,strlen(buffer));
if ( nb_spec_ali == 0 )
	{
	if ( strpbrk(buffer,"EFILPQ") != NULL )
		type_seq = 1;
	else
		type_seq = 0;
	}

if ( type_seq == 1 )
	{
	for ( pc = buffer ; *pc != 0 ; ++pc )
		{
		if ( *pc == 'B' )
			*pc = 'X';
		if ( *pc == 'Z' )
			*pc = 'X';
		if ( strchr("ACDEFGHIKLMNPQRSTVWY*-$ X?",*pc) == NULL )
			{
			sprintf(buffer,"Character %c is not an amino acid (on line %d)",*pc,nb_line);
			return(-1);
			}
		}
	}
else
	{
	for ( pc = buffer ; *pc != 0 ; ++pc )
		{
		if ( strchr("ACGUTN*-$ RYMWSKVHDBX?",*pc) == NULL )
			{
			sprintf(buffer,"Character %c is not a base (on line %d)",*pc,nb_line);
			return(-1);
			}
		}
	}
return(0);
}

/******************************************************************************/
/*  Reads aligned sequence files in ALI format        */
/******************************************************************************/
void read_file_ALI(int no_file_ali)
{
FILE *in_ali;
char *pc, *pb, *ps;
int i, j, k, nb_line, drap_erreur;

in_ali = ffopen(file_ALI[no_file_ali],"r");
fgets(buffer,LENBUFFER - 1,in_ali);
if ( *buffer != '#' )
	{
	fprintf(stderr,"First line does not begin with #\n");
	exit(1);
	}
*strchr(buffer,'\n') = 0;
strncpy(molecule,buffer + 1,LENMOLDOM - 1);
fgets(buffer,LENBUFFER - 1,in_ali);
if ( *buffer != '#' )
	{
	fprintf(stderr,"Second line does not begin with #\n");
	exit(1);
	}
*strchr(buffer,'\n') = 0;
strncpy(domaine,buffer + 1,LENMOLDOM - 1);
nb_line = 2;
for ( nb_line = 2 , nb_spec_ali = 0 , nb_nuc_tot_ali = 0 ; (drap_erreur = read_line(in_ali,0,&nb_line)) == 0 ; ++nb_spec_ali )
	{
	if ( nb_spec_ali >= MAXSPECIES )
		{
		fprintf(stderr,"Maximal number of species equals %d\n",MAXSPECIES);
		exit(1);
		}
	switch(read_line(in_ali,1,&nb_line))
		{
		case -2 :
			fprintf(stderr,"There is an identifier without a sequence at the end of the file\n");
			exit(1);
			break;
		case -1 :
			goto erreur;
			break;
		case 0 :
			if ( interp_seq(nb_line) == -1 )
				goto erreur;
			break;
		}
	}
if ( drap_erreur == -1 )
	goto erreur;
fclose(in_ali);

if ( nb_spec_ali == 0 )
	{
	fprintf(stderr,"The file %s is empty\n", file_ALI[no_file_ali]);
	exit(1);
	}
lenseq_tot = MAXCARSEQ / (long) nb_spec_ali;
if ( lenseq_tot < ((long) nb_nuc_tot_ali + 1) )
	{
	fprintf(stderr,"There are too many characters, %ld (MAX = %ld)\n",(long) nb_spec_ali * (long) (nb_nuc_tot_ali + 1),(long) MAXCARSEQ);
	exit(1);
	}

if ( type_seq == 0 )
	strcpy(car_inconnu," ?XRYMWSKVHDBN*-");
else
	strcpy(car_inconnu," ?X*-");

if ( (ident = (char *) malloc(nb_spec_ali * LENIDENT * sizeof(char))) == NULL )
	{
	fprintf(stderr,"Insufficient available memory (ident)\n");
	exit(-1);
	}

in_ali = ffopen(file_ALI[no_file_ali],"r");
for ( i = 0 ; fgets(buffer,LENBUFFER - 1,in_ali) != NULL ; )
	if ( *buffer == '>' )
		{
		if ( (pc = strstr(buffer,"#MASKED#")) == NULL )
			mask[i] = 0;
		else
			{
			mask[i] = -1;
			*pc = '\n';
			}
		if ( (pc = strstr(buffer,"#NEW#")) == NULL )
			new[i] = 0;
		else
			{
			new[i] = -1;
			*pc = '\n';
			}
		if ( strlen(buffer + 1) > LENIDENT )
			*(buffer + LENIDENT) = '\n';
		for ( pc = ident + i * LENIDENT , pb = buffer + 1 ; *pb != '\n' ; ++pc , ++pb )
			*pc = *pb;
		*pc = 0;
		}
	else
		if ( *buffer != '#' )
			{
			for ( ps = seq + (long) i * lenseq_tot , j = 0 , k = 0 , pc = buffer ; *pc != '\n' ; ++ps , ++pc , ++j )
				{
				if ( (*ps = *pc) == ' ' )
					++k;
				}
			for ( ; j < nb_nuc_tot_ali ; ++j , ++k , ++ps )
				*ps = ' ';
			*ps = 0;
			++i;
			}
fclose(in_ali);
return;

erreur:
	fclose(in_ali);
	fprintf(stderr,"%s\n",buffer);
	exit(-1);
}

/******************************************************************************/
void read_file_CLADE()
{
FILE *in_CLADE;
char *pc, *pf, stop;
int i, j;

if ( *file_clade == 0 )
	return;

nb_clade = 0;
in_CLADE = ffopen(file_clade,"r");
while ( fgets(buffer,LENBUFFER - 1,in_CLADE) != NULL )
	{
//fprintf(stdout,"%s",buffer); fflush(stdout);
	*(pc = strchr(buffer,':')) = 0;
	strcpy(name_clade[nb_clade],buffer);
	nb_species_clade[nb_clade] = 0;
	stop = 0;
	for ( ; stop == 0 ; pc = pf )
		{
		pf = strpbrk(pc + 1,",\n");
		if ( *pf == '\n' )
			stop = -1;
		*pf = 0;
		strcpy(species_clade[nb_clade][nb_species_clade[nb_clade]++],pc + 1);
		if ( nb_species_clade[nb_clade] >= MAXSPECPERCLADE )
			{
			fprintf(stderr,"Too many species for clade %s (MAX=%d)\n",name_clade[nb_clade],MAXSPECPERCLADE);
			exit(1);
			}
//fprintf(stdout,"pc + 1:%s|pf + 1:%s\n",pc + 1,pf + 1); fflush(stdout);
		}
	++nb_clade;
	}
fclose(in_CLADE);

/*for ( i = 0 ; i < nb_clade ; ++i )
	{
	fprintf(stdout,"[%d]%s:",i,name_clade[i]);
	for ( j = 0 ; j < nb_species_clade[i] ; ++j )
		fprintf(stdout,",%s",species_clade[i][j]);
	fprintf(stdout,"\n"); fflush(stdout);
	}
*/
}

/******************************************************************************/
void read_file_NEEDED_CLADE()
{
FILE *in_NEEDED_CLADE;
int i;
char *pc;

system("mkdir NOTOK");
if ( *file_needed_clade == 0 )
	{
	system("mkdir OK");
	return;
	}

if ( *file_clade == 0 )
	{
	fprintf(stderr,"You must specify a file_clade if you want to use file_needed_clade\n");
	exit(1);
	}

nb_needed_clade = 0;
in_NEEDED_CLADE = ffopen(file_needed_clade,"r");
while ( fgets(buffer,LENBUFFER - 1,in_NEEDED_CLADE) != NULL )
	{
//fprintf(stdout,"%s",buffer); fflush(stdout);
	if ( strchr(buffer,'[') == NULL  )
		{
		fprintf(stderr,"A '[' is missing on this line to indicate the minimum number of species (%s):%s",file_clade,buffer);
		exit(1);
		}
	*(pc = strchr(buffer,'[')) = 0;
	sscanf(pc + 1,"%d]",&min_spec_needed_clade[nb_needed_clade]);
	for ( i = 0 ; i < nb_clade ; ++i )
		if ( strcmp(name_clade[i],buffer) == 0 )
			break;
	if ( i == nb_clade )
		{
		fprintf(stderr,"The clade, %s, is missing in %s\n",buffer,file_clade);
		exit(1);
		}
	else
		no_needed_clade[nb_needed_clade++] = i;
//fprintf(stdout,"-> %s %d\n",name_clade[no_needed_clade[nb_needed_clade -1]],min_spec_needed_clade[nb_needed_clade - 1]); fflush(stdout);
	}
fclose(in_NEEDED_CLADE);

//for ( i = 0 ; i < nb_needed_clade ; ++i )
//	fprintf(stdout,"NEEDED:%s",name_clade[no_needed_clade[i]]);
for ( i = 0 ; i <= nb_needed_clade ; ++i )
	{
	sprintf(buffer,"mkdir OK-%d-needed-clades",i);
	system(buffer);
	}
}

/******************************************************************************/
int statistics_taxonomy(int no_file_ali)
{
int i, j, k, nb_clade_present, no_spec_clade[MAXCLADE], nb_para, nb_needed_clade_present;

for ( i = 0 , nb_spec_uniq = 0 ; i < nb_spec_ali ; ++i )
	{
	strcpy(buffer,ident + i * LENIDENT);
	*strchr(buffer,'@') = 0;
	for ( j = 0 ; j < nb_spec_uniq ; ++j )
		if ( strcmp(buffer,ident_uniq[j]) == 0 )
			break;
	if ( j == nb_spec_uniq )
		{
		strcpy(ident_uniq[nb_spec_uniq],buffer);
		nb_ident_uniq[nb_spec_uniq] = 1;
		++nb_spec_uniq;
		}
	else
		++nb_ident_uniq[j];
	}

nb_clade_present = 0;
for ( i = 0 ; i < nb_clade ; ++i )
	{
	no_spec_clade[i] = 0;
	for ( j = 0 ; j < nb_species_clade[i] ; ++j )
		{
		for ( k = 0 ; k < nb_spec_uniq ; ++k )
			if ( strcmp(species_clade[i][j],ident_uniq[k]) == 0 )
				break;
		if ( k < nb_spec_uniq )
			++no_spec_clade[i];
		}
	if ( no_spec_clade[i] > 0 )
		++nb_clade_present;
	}

fprintf(out_STAT,"%s",file_ALI[no_file_ali]);
for ( i = 0 ; i < nb_clade ; ++i )
	fprintf(out_STAT,"\t%d",no_spec_clade[i]);
fprintf(out_STAT,"\n");

nb_para = 0;
for ( j = 0 ; j < nb_spec_uniq ; ++j )
	nb_para += nb_ident_uniq[j] - 1;

fprintf(stdout,"%d species, %d paralogs",nb_spec_uniq,nb_para);
if ( *file_clade != 0 )
	fprintf(stdout,", %d clades\n",nb_clade_present);

//for ( j = 0 ; j < nb_spec_uniq ; ++j )
//	fprintf(stdout,"%s\t%d\n",ident_uniq[j],nb_ident_uniq[j]);
//for ( i = 0 ; i < nb_clade ; ++i )
//	fprintf(stdout,"%s:%d\n",name_clade[i],no_spec_clade[i]);

if ( (min_spec != 0) && (nb_spec_uniq < min_spec) )
	{
	fprintf(stdout," --> too few species (%d)\n",nb_spec_uniq);
	return(-1);
	}
if ( nb_para > max_para )
	{
	fprintf(stdout," --> too many paralogs (%d)\n",nb_para);
	return(-1);
	}
if ( *file_needed_clade != 0 )
	{
	nb_needed_clade_present = 0;
	for ( i = 0 ; i < nb_needed_clade ; ++i )
//		if ( no_spec_clade[no_needed_clade[i]] == 0 )

		if ( no_spec_clade[no_needed_clade[i]] < min_spec_needed_clade[i] )
			{
			fprintf(stdout," --> the clade %s is not sufficiently rich (%d species < %d)\n",name_clade[no_needed_clade[i]],no_spec_clade[no_needed_clade[i]],min_spec_needed_clade[i]);
//			return(-1);
			}
		else
			++nb_needed_clade_present;
	fprintf(stdout," --> %d needed clades OK\n",nb_needed_clade_present);
	sprintf(buffer,"mv %s OK-%d-needed-clades",file_ALI[no_file_ali],nb_needed_clade_present);
	system(buffer);
	return(0);
	}
else
	if ( (min_clade != 0) && (nb_clade_present < min_clade) )
		{
		fprintf(stdout," --> too few clades (%d)\n",nb_clade_present);
		return(-1);
		}

	/*
int min_spec, min_clade, nb_clade, nb_needed_clade, max_para, nb_file_OK;
char file_clade[100], file_needed_clade[100];
char name_clade[MAXCLADE];
int nb_species_clade[MAXCLADE];
char species_clade[MAXCLADE][MAXSPECPERCLADE];
int no_needed_clade[MAXCLADE];
*/
fprintf(stdout," --> moved to the directory OK\n");
++nb_file_OK;
sprintf(buffer,"mv %s OK",file_ALI[no_file_ali]);
system(buffer);
return(0);
}

/******************************************************************************/
void usage(void)
{
fprintf(stderr,"USAGE: select-aliMCL-taxonomy");
fprintf(stderr,"\n\tmin_spec=<minimum number of species>");
fprintf(stderr,"\n\tmin_clade=<minimum number of clades> (optional)");
fprintf(stderr,"\n\tmax_para=<maximum number of paralogs (i.e. >)> (optional)");
fprintf(stderr,"\n\tfile_clade=<file containing the definition of clades> (optional)");
fprintf(stderr,"\n\tfile_needed_clade=<file containing the list of clades that must be present> (optional)");
fprintf(stderr,"\n\tali=<name of file .ali> (optional)\n\n");
exit(1);
}

/****************************************************************************/
int main(int argc, char **argv)
{
int i;

if ( argc < 2 )
	usage();

if ( (buffer = (char *) malloc(LENBUFFER * sizeof(char))) == NULL )
	{
	fprintf(stderr,"Insufficient available memory\n\n");
	exit(1);
	}

ident = NULL;
nb_file_ali = 0;
min_spec = 0;
min_clade = 0;
max_para = 10000000;
*file_clade = 0;
nb_clade = 0;
*file_needed_clade = 0;
nb_needed_clade = 0;
nb_file_OK = 0;
for( i = 1 ; i < argc ; ++i )
	{
	if ( strncmp(argv[i],"ali=",4) == 0 )
		{
		strcpy(file_ALI[0],argv[i] + 4);
		nb_file_ali = 1;
		if ( strchr(file_ALI[0],'.') != NULL )
			*strchr(file_ALI[0],'.') = 0;
		strcat(file_ALI[0],".ali");
		}
	else
		if ( strncmp(argv[i],"min_spec=",9) == 0 )
			sscanf(argv[i],"min_spec=%d",&min_spec);
		else
			if ( strncmp(argv[i],"min_clade=",10) == 0 )
				sscanf(argv[i],"min_clade=%d",&min_clade);
			else
				if ( strncmp(argv[i],"file_clade=",11) == 0 )
					strcpy(file_clade,argv[i] + 11);
				else
					if ( strncmp(argv[i],"max_para=",9) == 0 )
						sscanf(argv[i],"max_para=%d",&max_para);
					else
						if ( strncmp(argv[i],"file_needed_clade=",18) == 0 )
							strcpy(file_needed_clade,argv[i] + 18);
						else
							{
							fprintf(stderr,"This argument is not correctly formatted: %s\n\n",argv[i]);
							usage();
							}
	}
if ( (min_spec < 0) ||  (min_clade < 0) || (max_para < 0) )
	{
	fprintf(stderr,"You must indicate a positive value for min_spec,min_clade and max_para\n");
	exit(1);
	}

if ( nb_file_ali == 0 )
	read_name_file_ali();
//printf("OK1"); fflush(stdout);
read_file_CLADE();
//printf("OK2"); fflush(stdout);
read_file_NEEDED_CLADE();

out_STAT = ffopen("statistics-ALI.xls","w");
fprintf(out_STAT,"FILE");
for ( i = 0 ; i < nb_clade ; ++i )
	fprintf(out_STAT,"\t%s(%d)",name_clade[i],nb_species_clade[i]);
fprintf(out_STAT,"\n");
for ( i = 0 ; i < nb_file_ali ; ++i )
	{
	read_file_ALI(i);
	fprintf(stdout,"%s (%d sequences and %d characters) -> ",file_ALI[i],nb_spec_ali,nb_nuc_tot_ali);
	if ( statistics_taxonomy(i) == -1 )
		{
		sprintf(buffer,"mv %s NOTOK",file_ALI[i]);
		system(buffer);
		}
	libere_memoire();
//	fprintf(stdout,"\n");
	}
fclose(out_STAT);
fprintf(stdout,"There are %d files that are compliant with the constraints\n",nb_file_OK);
exit(0);
}
