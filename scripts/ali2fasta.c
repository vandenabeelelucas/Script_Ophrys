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

#define LENIDENT		80
#define LENSEQ			1000000000
#define MAXFILEALI		25000
#define LENBUFFER		LENSEQ + 1

/* list of presentation files */
char *buffer;
char file_ALI[MAXFILEALI][100];
int nb_file_ALI, type_seq, min_length, nb_seq_removed;
static char seq[LENSEQ];
char opt_keep_ali, opt_ident;

static char Mois[12][10] = { "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"};
static char Jour[7][10] = { "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"};

void date(char *version);
void date_fichier(char *nom, char *version);
void erreur_fatale(char *message);
FILE *ffopen(char *File, char *Mode);

/******************************************************************************/
void strupr(char *chaine)
{
char *pc;

for ( pc = chaine ; *pc != 0 ; ++pc )
    *pc = toupper(*pc);
}

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
    struct tm *date;
    struct stat buf;

    stat(nom,&buf);
    date = localtime(&buf.st_atime);
    sprintf(version,"%d %s %d at %d hours %d",date->tm_mday,Mois[date->tm_mon],1900 + date->tm_year,date->tm_hour,date->tm_min);
}

/******************************************************************************/
/* Displays the fatal error message and exits the program       */
/******************************************************************************/
void erreur_fatale(char *message)
{
    fprintf(stderr,"The program has encountered the following FATAL ERROR :\n\n");
    fprintf(stderr,"%s\n\n\n",message);
    fprintf(stderr,"\n\n\n");
    exit(1);
}

/******************************************************************************/
/* Opens a file and verifies that it opened successfully                    */
/******************************************************************************/
FILE *ffopen(char *File, char *Mode)
{
    FILE *in;

    if ( (in = fopen(File,Mode)) == NULL ) {
        fprintf(stderr,"Problem in opening file %s\n\n",File);
        exit(1);
    }
    else
        return(in);
}

/******************************************************************************/
/* Reads the list of CARAC_FILE type files                   */
/******************************************************************************/
void read_name_file_ALI()
{
FILE *in;

sprintf(buffer,"ls -1 *.ali >ali2paup.tmp");
system(buffer);

nb_file_ALI = 0;
in = ffopen("ali2paup.tmp","r");
while ( fgets(buffer,1000,in) != NULL )
	{
	if ( nb_file_ALI >= MAXFILEALI )
		{
		fprintf(stderr,"There are too many files for this program (MAX=%d)\n",MAXFILEALI);
		exit(1);
		}
	*strchr(buffer,'\n') = 0;
	strcpy(file_ALI[nb_file_ALI++],buffer);
	}
fclose(in);
system("\"rm\" ali2paup.tmp");

if ( nb_file_ALI == 0 )
	{
	fprintf(stderr,"There is no file for this program\n");
	exit(1);
	}
}

/******************************************************************************/
void check_sequence(int no_file_ALI, int *nb_char)
{
char *pc;

//fprintf(stdout,"{%s}\n",seq); fflush(stdout);
strupr(seq);
*nb_char = 0;
if ( type_seq == 1 )
	{
	for ( pc = seq ; *pc != 0 ; ++pc )
		{
		if ( *pc == 'B' )
			*pc = 'X';
		if ( *pc == 'Z' )
			*pc = 'X';
		if ( strchr("ACDEFGHIKLMNPQRSTVWY*-$ X?",*pc) == NULL )
			{
			fprintf(stderr,"Character %c is not an amino acid (%s)\n",*pc,file_ALI[no_file_ALI]);
			*pc = 'X';
			}
		if ( strchr(" ?*-",*pc) == NULL )
			++*nb_char;
		}
	}
else
	{
	for ( pc = seq ; *pc != 0 ; ++pc )
		{
		if ( strchr("ACGUTN*-$ RYMWSKVHDBX?",*pc) == NULL )
			{
			fprintf(stderr,"Character %c is not a nucleotide (%s)\n",*pc,file_ALI[no_file_ALI]);
			*pc = 'N';
			}
		if ( strchr(" ?*-",*pc) == NULL )
			++*nb_char;
		}
	}
}

/******************************************************************************/
/* Reads aligned sequence files in ALI format         */
/******************************************************************************/
void ALI2FASTA(int no_file_ALI)
{
FILE *out_fasta, *in_ali;
int nb_char;
char *ps, ident[200];

nb_seq_removed = 0;
strcpy(buffer,file_ALI[no_file_ALI]);
if ( strchr(buffer,'.') != NULL )
	*strchr(buffer,'.') = 0;
strcat(buffer,".fasta");
fprintf(stdout,"Reformating %s into %s",file_ALI[no_file_ALI],buffer);
out_fasta = ffopen(buffer,"w");
in_ali = ffopen(file_ALI[no_file_ALI],"r");
type_seq = -1;
*seq = 0;
while ( fgets(buffer,LENBUFFER,in_ali) != NULL )
	{
//fprintf(stdout,"%s",buffer); fflush(stdout);
	if ( *buffer != '#' )
		{
		if ( *buffer == '>' )
			{
			if ( *seq != 0 )
				{
				if ( type_seq == -1 )
					{
					if ( strpbrk(seq,"EFILPQ") != NULL )
						type_seq = 1;
					else
						type_seq = 0;
					}
				check_sequence(no_file_ALI,&nb_char);
				if ( nb_char >= min_length )
					{
					fprintf(out_fasta,"%s\n",ident);
					if ( opt_keep_ali == 1 )
						fprintf(out_fasta,"%s\n",seq);
					else
						{
						for ( ps = seq ; *ps != 0 ; ++ps )
							if ( strchr(" ?*-",*ps) == NULL )
								fprintf(out_fasta,"%c",*ps);
						fprintf(out_fasta,"\n");
						}
					}
				else
					++nb_seq_removed;
				*seq = 0;
				}
			*strchr(buffer,'\n') = 0;
			if ( strchr(buffer,' ') == NULL )
				{
				if ( strchr(buffer,'_') != NULL )
					*strchr(buffer,'_') = ' ';
				else
					strcat(buffer," sp");
				}
			if ( strchr(buffer,'@') == NULL )
				strcat(buffer,"@auto");
			if ( opt_ident != 0 )
				for ( ps = buffer ; *ps != 0 ; ++ps )
					if ( *ps == ' ' )
						*ps = '_';
			strcpy(ident,buffer);
/*
			if ( opt_ident == 0 )
				fprintf(out_fasta,"%s\n",buffer);
			else
				{
				for ( ps = buffer ; *ps != 0 ; ++ps )
					if ( *ps == ' ' )
						fprintf(out_fasta,"_");
					else
						fprintf(out_fasta,"%c",*ps);
				fprintf(out_fasta,"\n");
				}
*/
			}
		else
			{
			*strchr(buffer,'\n') = 0;
			strcat(seq,buffer);
			}
		}
	}
check_sequence(no_file_ALI,&nb_char);
if ( nb_char >= min_length )
	{
	fprintf(out_fasta,"%s\n",ident);
	if ( opt_keep_ali == 1 )
		fprintf(out_fasta,"%s\n",seq);
	else
		{
		for ( ps = seq ; *ps != 0 ; ++ps )
			if ( strchr(" ?*-",*ps) == NULL )
				fprintf(out_fasta,"%c",*ps);
		}
	}
else
	++nb_seq_removed;
fclose(out_fasta);
fclose(in_ali);

if ( min_length > 0 )
	fprintf(stdout," --> %d sequences removed\n",nb_seq_removed);
else
	fprintf(stdout,"\n");

//sprintf(buffer,"\"mv\" %s %s",file_FASTA[no_file_fasta],file_FASTA[no_file_fasta]);
//*(buffer + strlen(buffer) - 3) = 0;
//strcat(buffer,"bak");
//system(buffer);
//sprintf(buffer,"\"mv\" file_FASTA %s",file_FASTA[no_file_fasta]);
//system(buffer);
}

/******************************************************************************/
void usage(void)
{
fprintf(stderr,"USAGE: ali2fasta\n");
fprintf(stderr,"keep_alignment=<yes or no> (default = yes)\n");
fprintf(stderr,"ident=<space or no> (default = space)\n");
fprintf(stderr,"min_length=<minimum number of characters to keep a sequence (default=0)> (optional)\n");
fprintf(stderr,"ali=<name of file under the MUST ED format (.ali)> (optional)\n\n");
exit(1);
}

/****************************************************************************/
int main(int argc, char **argv)
{
int i;

if ( (buffer = (char *) malloc(LENBUFFER * sizeof(char))) == NULL )
	erreur_fatale("Insufficient available memory");

nb_file_ALI = 0;
opt_keep_ali = 1;
opt_ident = 0;
min_length = 0;
for( i = 1 ; i < argc ; ++i )
	{
	if ( strncmp(argv[i],"ali=",4) == 0 )
		{
		strcpy(file_ALI[0],argv[i] + 4);
		nb_file_ALI = 1;
		if ( strchr(file_ALI[0],'.') != NULL )
			*strchr(file_ALI[0],'.') = 0;
		strcat(file_ALI[0],".ali");
		}
	else
		if ( strncmp(argv[i],"keep_alignment=yes",18) == 0 )
			opt_keep_ali = 1;
		else
			if ( strncmp(argv[i],"keep_alignment=no",17) == 0 )
				opt_keep_ali = 0;
			else
				if ( strncmp(argv[i],"ident=space",11) == 0 )
					opt_ident = 0;
				else
					if ( strncmp(argv[i],"ident=no",8) == 0 )
						opt_ident = 1;
					else
						if ( strncmp(argv[i],"min_length=",11) == 0 )
							sscanf(argv[i],"min_length=%d",&min_length);
						else
							{
							fprintf(stderr,"This argument is not correctly formatted: %s\n\n",argv[i]);
							usage();
							}
	}


if ( nb_file_ALI == 0 )
	read_name_file_ALI();

for ( i = 0 ; i < nb_file_ALI ; ++i )
	ALI2FASTA(i);
exit(0);
}
